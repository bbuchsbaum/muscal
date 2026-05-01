#' Multifer Component Tests for COVSTATIS
#'
#' Runs adapter-owned `multifer` component-significance tests for either the
#' COVSTATIS interstructure/RV axes across input tables or the compromise ROI
#' eigencomponents.
#'
#' @param x A fitted `covstatis` object.
#' @param data Original list of covariance/correlation matrices used to fit `x`.
#' @param axes Which component family to test: `"interstructure"` tests axes of
#'   the table-by-table RV/interstructure matrix, `"compromise"` tests ROI
#'   eigencomponents of the weighted compromise matrix, and `"both"` runs both
#'   tests and returns a named list.
#' @param B Number of Monte Carlo draws per ladder rung.
#' @param alpha Significance threshold passed to [multifer::infer()].
#' @param seed Optional random seed.
#' @param strict Logical; passed to [multifer::infer()].
#' @param ... Additional arguments passed to [multifer::infer()].
#'
#' @return A `multifer::infer_result` object, or a named list of two such
#'   objects when `axes = "both"`.
#' @export
infer_covstatis <- function(x,
                            data,
                            axes = c("compromise", "interstructure", "both"),
                            B = 1000L,
                            alpha = 0.05,
                            seed = NULL,
                            strict = TRUE,
                            ...) {
  if (!inherits(x, "covstatis")) {
    stop("`x` must be a fitted covstatis object.", call. = FALSE)
  }
  if (!requireNamespace("multifer", quietly = TRUE)) {
    stop("`infer_covstatis()` requires the multifer package.", call. = FALSE)
  }
  .covstatis_require_multifer_adapter_geometry()

  axes <- match.arg(axes)
  if (identical(axes, "both")) {
    out <- list(
      interstructure = infer_covstatis(
        x, data = data, axes = "interstructure", B = B, alpha = alpha,
        seed = seed, strict = strict, ...
      ),
      compromise = infer_covstatis(
        x, data = data, axes = "compromise", B = B, alpha = alpha,
        seed = seed, strict = strict, ...
      )
    )
    class(out) <- c("covstatis_multifer_result", "list")
    return(out)
  }

  payload <- .covstatis_multifer_payload(x, data)
  adapter <- .covstatis_multifer_adapter(axes)
  model <- .covstatis_adapter_fit(payload, axes = axes)

  multifer::infer(
    adapter = adapter,
    data = payload,
    geometry = "adapter",
    relation = "variance",
    targets = "component_significance",
    model = model,
    B = B,
    R = 0L,
    alpha = alpha,
    seed = seed,
    strict = strict,
    ...
  )
}

.covstatis_multifer_adapter <- function(axes = c("compromise", "interstructure")) {
  axes <- match.arg(axes)
  .covstatis_require_multifer_adapter_geometry()

  multifer::infer_adapter(
    adapter_id = paste0("muscal_covstatis_", axes),
    adapter_version = "0.1.0",
    shape_kinds = "adapter",
    capabilities = multifer::capability_matrix(
      list(
        geometry = "adapter",
        relation = "variance",
        targets = "component_significance"
      )
    ),
    roots = function(x, ...) {
      x$roots
    },
    refit = function(x, new_data, ...) {
      .covstatis_adapter_fit(new_data, axes = axes)
    },
    null_action = function(x, data, ...) {
      .covstatis_roi_permutation_null(data)
    },
    component_stat = function(x, data, k, ...) {
      if (!identical(as.integer(k), 1L)) {
        stop("COVSTATIS adapter statistics expect `k = 1` on the current residual.",
             call. = FALSE)
      }
      .covstatis_leading_root_ratio(data, axes = axes)
    },
    residualize = function(x, k, data, ...) {
      if (!identical(as.integer(k), 1L)) {
        stop("COVSTATIS adapter residualization expects `k = 1` on the current residual.",
             call. = FALSE)
      }
      .covstatis_residualize_payload(x, data, axes = axes)
    },
    validity_level = "conditional",
    declared_assumptions = c(
      "tables_exchangeable",
      "roi_labels_aligned_across_tables",
      "roi_permutation_null"
    ),
    checked_assumptions = list(
      list(
        name = "covstatis_payload",
        check = function(data, ...) {
          .covstatis_valid_payload(data)
        },
        detail = "COVSTATIS adapter data must be a list payload with finite symmetric square matrices of common dimension."
      )
    )
  )
}

.covstatis_multifer_payload <- function(x, data) {
  if (is.null(data)) {
    stop("`data` is required for COVSTATIS multifer inference.", call. = FALSE)
  }
  if (!is.list(data) || length(data) < 2L) {
    stop("`data` must be a list of at least two covariance matrices.",
         call. = FALSE)
  }
  tables <- lapply(data, function(mat) .pre_process_new_cov(x, as.matrix(mat)))
  names(tables) <- names(data) %||% x$block_labels %||% paste0("Table_", seq_along(tables))
  list(
    tables = tables,
    ncomp = ncol(stats::coef(x)),
    labels = x$labels
  )
}

.covstatis_adapter_fit <- function(data, axes = c("compromise", "interstructure")) {
  axes <- match.arg(axes)
  if (!isTRUE(.covstatis_valid_payload(data))) {
    stop("Invalid COVSTATIS adapter payload.", call. = FALSE)
  }

  if (identical(axes, "interstructure")) {
    M <- .covstatis_table_matrix(data$tables)
    C <- crossprod(M)
    eig <- eigen((C + t(C)) / 2, symmetric = TRUE)
    roots <- .covstatis_positive_roots(eig$values, data$ncomp)
    vectors <- eig$vectors[, seq_len(length(roots)), drop = FALSE]
    return(structure(
      list(
        roots = roots,
        v = vectors,
        eig = eig,
        axes = axes
      ),
      class = c("covstatis_multifer_fit", "list")
    ))
  }

  inter <- .covstatis_interstructure(data$tables)
  compromise <- Reduce("+", Map(`*`, data$tables, inter$alpha))
  compromise <- (compromise + t(compromise)) / 2
  eig <- eigen(compromise, symmetric = TRUE)
  roots <- .covstatis_positive_roots(eig$values, data$ncomp)
  vectors <- eig$vectors[, seq_len(length(roots)), drop = FALSE]
  structure(
    list(
      roots = roots,
      v = vectors,
      eig = eig,
      alpha = inter$alpha,
      compromise = compromise,
      axes = axes
    ),
    class = c("covstatis_multifer_fit", "list")
  )
}

.covstatis_interstructure <- function(tables) {
  M <- .covstatis_table_matrix(tables)
  C <- crossprod(M)
  eigC <- eigen((C + t(C)) / 2, symmetric = TRUE)
  alpha_vec <- eigC$vectors[, 1L]
  alpha_vec <- alpha_vec * sign(sum(alpha_vec))
  if (abs(sum(alpha_vec)) < 1e-12) {
    alpha_vec <- abs(alpha_vec)
  }
  alpha <- alpha_vec / sum(alpha_vec)
  list(C = C, eigC = eigC, alpha = alpha)
}

.covstatis_table_matrix <- function(tables) {
  vapply(tables, c, numeric(length(tables[[1L]])), USE.NAMES = FALSE)
}

.covstatis_positive_roots <- function(values, ncomp) {
  values <- as.numeric(values)
  values[!is.finite(values)] <- 0
  values <- pmax(values, 0)
  n <- min(as.integer(ncomp), length(values))
  if (n < 1L) {
    stop("COVSTATIS adapter fit has no component roots.", call. = FALSE)
  }
  values[seq_len(n)]
}

.covstatis_leading_root_ratio <- function(data, axes) {
  roots <- .covstatis_adapter_fit(data, axes = axes)$roots
  denom <- sum(roots)
  if (!is.finite(denom) || denom <= .Machine$double.eps) {
    return(0)
  }
  roots[[1L]] / denom
}

.covstatis_roi_permutation_null <- function(data) {
  out <- data
  out$tables <- lapply(data$tables, function(mat) {
    idx <- sample.int(nrow(mat))
    mat[idx, idx, drop = FALSE]
  })
  names(out$tables) <- names(data$tables)
  out
}

.covstatis_residualize_payload <- function(fit, data, axes) {
  out <- data
  if (identical(axes, "interstructure")) {
    g <- fit$v[, 1L, drop = FALSE]
    M <- .covstatis_table_matrix(data$tables)
    score <- M %*% g
    M_res <- M - score %*% t(g)
    p <- nrow(data$tables[[1L]])
    out$tables <- lapply(seq_len(ncol(M_res)), function(i) {
      mat <- matrix(M_res[, i], nrow = p, ncol = p)
      (mat + t(mat)) / 2
    })
    names(out$tables) <- names(data$tables)
    return(out)
  }

  v <- fit$v[, 1L, drop = FALSE]
  P <- diag(nrow(v)) - tcrossprod(v)
  out$tables <- lapply(data$tables, function(mat) {
    res <- P %*% mat %*% P
    (res + t(res)) / 2
  })
  names(out$tables) <- names(data$tables)
  out
}

.covstatis_valid_payload <- function(data) {
  if (!is.list(data) || is.null(data$tables) || !is.list(data$tables) ||
      length(data$tables) < 2L) {
    return(FALSE)
  }
  ok <- vapply(data$tables, function(mat) {
    is.matrix(mat) &&
      is.numeric(mat) &&
      nrow(mat) == ncol(mat) &&
      all(is.finite(mat)) &&
      isTRUE(isSymmetric(mat, tol = sqrt(.Machine$double.eps)))
  }, logical(1L))
  if (!all(ok)) {
    return(FALSE)
  }
  dims <- vapply(data$tables, nrow, integer(1L))
  if (length(unique(dims)) != 1L) {
    return(FALSE)
  }
  is.numeric(data$ncomp) && length(data$ncomp) == 1L && data$ncomp >= 1L
}

.covstatis_require_multifer_adapter_geometry <- function() {
  if (!.covstatis_multifer_adapter_geometry_available()) {
    stop(
      "`infer_covstatis()` requires a multifer build with `geometry = \"adapter\"` support.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.covstatis_multifer_adapter_geometry_available <- function() {
  if (!requireNamespace("multifer", quietly = TRUE)) {
    return(FALSE)
  }
  tryCatch({
    multifer::capability_matrix(
      list(
        geometry = "adapter",
        relation = "variance",
        targets = "component_significance"
      )
    )
    TRUE
  }, error = function(e) {
    FALSE
  })
}
