#' Test leading-root concentration in COVSTATIS
#'
#' Test whether the leading eigenvalue accounts for an unusually large share
#' of the spectrum under independent ROI relabeling across tables.
#'
#' @param x A fitted `covstatis` object.
#' @param data Original list of covariance/correlation matrices used to fit `x`.
#' @param axes Matrix to test. `"interstructure"` uses the table-by-table
#'   matrix of inner products, `"compromise"` uses the weighted compromise,
#'   and `"both"` returns one test for each.
#' @param B Maximum number of Monte Carlo draws per test.
#' @param alpha Significance threshold passed to [multifer::infer()].
#' @param seed Optional random seed.
#' @param strict Logical; passed to [multifer::infer()].
#' @param ... Additional arguments passed to [multifer::infer()].
#' @param max_steps Only the leading root can be tested (`1`, the default).
#' @param auto_subspace Group nearly tied leading roots in the reported unit.
#'   Roots tied to numerical precision always share a unit.
#' @param tie_threshold Relative eigenvalue gap used to group near ties.
#'
#' @details The statistic is the largest positive eigenvalue divided by the
#'   sum of **all** positive eigenvalues of the matrix being tested, or zero
#'   when there is no positive inertia. Both the observed statistic and every
#'   null draw use the full spectrum, regardless of the fitted `ncomp`.
#'
#'   Conditional on each preprocessed table up to ROI labeling, the null
#'   assumes the ROI labels are independently exchangeable within each table.
#'   Each draw applies the same permutation to rows and columns of a table and
#'   recomputes the interstructure and, for the compromise, its table weights.
#'   With unequal variable weights this assumption concerns the weighted tables.
#'
#'   The p-value tests concentration under this global correspondence null.
#'   The unit's `selected` flag records rejection of this global null. When
#'   leading roots are tied, the reported unit spans the tied eigenspace and
#'   the p-value still concerns the largest root alone. It supplies no estimate
#'   of how many dimensions to retain. Subsequent axes are descriptive.
#'   Independent ROI permutations of deflated tables can
#'   move their shared removed subspace. Conditional dimension selection
#'   would require a separate residual null.
#'
#' @return A `multifer::infer_result` object, or a named list of two such
#'   objects when `axes = "both"`.
#' @md
#' @export
infer_covstatis <- function(x,
                            data,
                            axes = c("compromise", "interstructure", "both"),
                            B = 1000L,
                            alpha = 0.05,
                            seed = NULL,
                            strict = TRUE,
                            ..., max_steps = 1L, auto_subspace = TRUE, tie_threshold = 0.01) {
  if (inherits(x, "distatis")) {
    stop("DISTATIS inference is not implemented; infer_covstatis() uses a different interstructure.",
         call. = FALSE)
  }
  if (!inherits(x, "covstatis")) {
    stop("`x` must be a fitted covstatis object.", call. = FALSE)
  }
  if (!requireNamespace("multifer", quietly = TRUE)) {
    stop("`infer_covstatis()` requires the multifer package.", call. = FALSE)
  }
  .covstatis_require_multifer_adapter_geometry()
  if (!is.null(max_steps) &&
      (!is.numeric(max_steps) || length(max_steps) != 1L || !isTRUE(max_steps == 1))) {
    stop("Only the leading root is supported. Use `max_steps = 1`; residual tests require a conditional null.",
         call. = FALSE)
  }

  axes <- match.arg(axes)
  if (identical(axes, "both")) {
    families <- c("interstructure", "compromise")
    out <- lapply(families, function(family) {
      infer_covstatis(x, data, axes = family, B = B, alpha = alpha, seed = seed,
                      strict = strict, ..., auto_subspace = auto_subspace,
                      tie_threshold = tie_threshold)
    })
    names(out) <- families
    class(out) <- c("covstatis_multifer_result", "list")
    return(out)
  }

  payload <- .covstatis_multifer_payload(x, data)
  adapter <- .covstatis_multifer_adapter(axes)
  model <- .covstatis_adapter_fit(payload, axes = axes)

  out <- multifer::infer(
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
    auto_subspace = auto_subspace,
    tie_threshold = tie_threshold,
    ...
  )

  # The engine sees one testable root. Use the full spectrum to describe ties
  # at that root, including a tie across the original fit's retention boundary.
  resolution <- length(model$roots) * .Machine$double.eps
  gap <- if (isTRUE(auto_subspace)) max(resolution, tie_threshold) else resolution
  units <- multifer::form_units(.covstatis_relative_roots(model$roots),
                                group_near_ties = TRUE, tie_threshold = gap)
  bundled <- inherits(out, "multifer_infer_bundle")
  result <- if (bundled) out$result else out
  unit <- units[1L, , drop = FALSE]
  unit$selected <- result$units$selected
  result$units <- unit
  labels <- list(statistic = "largest positive eigenvalue / total positive inertia",
                 null = "independent ROI relabeling within each preprocessed table",
                 estimand = paste("leading-root concentration of the", axes, "matrix"))
  result$mc$statistic_label <- labels$statistic
  result$mc$null_label <- result$component_tests$null_label <- labels$null
  result$mc$estimand_label <- labels$estimand
  if (!bundled) return(result)
  out$result <- result
  out$artifacts$engine$units <- unit
  out$artifacts$engine$roots_observed <- model$roots
  out$artifacts$engine$labels <- labels
  out
}

.covstatis_multifer_adapter <- function(axes = c("compromise", "interstructure")) {
  axes <- match.arg(axes)
  .covstatis_require_multifer_adapter_geometry()

  multifer::infer_adapter(
    adapter_id = paste0("muscal_covstatis_", axes),
    adapter_version = "0.2.0",
    shape_kinds = "adapter",
    capabilities = multifer::capability_matrix(
      list(
        geometry = "adapter",
        relation = "variance",
        targets = "component_significance"
      )
    ),
    roots = function(x, ...) {
      # multifer uses the number of exposed roots as the default ladder depth.
      x$roots[1L]
    },
    domains = function(x, data = NULL, ...) {
      "global"
    },
    loadings = function(x, domain = NULL, ...) {
      if (!is.null(domain) && !identical(domain, "global")) {
        stop("COVSTATIS multifer adapter only exposes the `global` domain.",
             call. = FALSE)
      }
      as.matrix(x$v)
    },
    variable_stat = function(x, data, domain = NULL, k, ...) {
      if (!is.null(domain) && !identical(domain, "global")) {
        stop("COVSTATIS multifer adapter only exposes the `global` domain.",
             call. = FALSE)
      }
      .muscal_squared_loading_stat(as.matrix(x$v), k)
    },
    refit = function(x, new_data, ...) {
      .covstatis_adapter_fit(new_data, axes = axes)
    },
    null_action = function(x, data, ...) {
      .covstatis_roi_permutation_null(data)
    },
    component_stat = function(x, data, k, ...) {
      if (!is.numeric(k) || length(k) != 1L || !isTRUE(k == 1)) {
        stop("The COVSTATIS correspondence test supports only the leading root.",
             call. = FALSE)
      }
      .covstatis_leading_root_ratio(data, axes = axes)
    },
    residualize = function(x, k, data, ...) {
      stop("COVSTATIS residual tests require a null that preserves the removed subspace.",
           call. = FALSE)
    },
    validity_level = "conditional",
    declared_assumptions = c(
      "roi_labels_aligned_across_tables",
      "roi_labels_independently_exchangeable_within_preprocessed_tables",
      "leading_root_only",
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
  list(tables = tables, labels = x$labels)
}

.covstatis_adapter_fit <- function(data, axes = c("compromise", "interstructure")) {
  axes <- match.arg(axes)
  if (!isTRUE(.covstatis_valid_payload(data))) {
    stop("Invalid COVSTATIS adapter payload.", call. = FALSE)
  }

  M <- .covstatis_table_matrix(data$tables)
  eig <- eigen(crossprod(M), symmetric = TRUE)
  fit <- list(axes = axes)
  if (axes == "compromise") {
    a <- eig$vectors[, 1L]
    if (abs(sum(a)) < 1e-12) a <- abs(a)
    fit$alpha <- a / sum(a)
    compromise <- Reduce("+", Map(`*`, data$tables, fit$alpha))
    fit$compromise <- (compromise + t(compromise)) / 2
    eig <- eigen(fit$compromise, symmetric = TRUE)
  }
  fit$eig <- eig
  fit$roots <- pmax(eig$values, 0)
  fit$v <- eig$vectors
  structure(fit, class = c("covstatis_multifer_fit", "list"))
}

.covstatis_relative_roots <- function(roots) {
  if (roots[[1L]] == 0) return(roots)
  roots / roots[[1L]]
}

.covstatis_leading_root_ratio <- function(data, axes) {
  fit <- .covstatis_adapter_fit(data, axes = axes)
  roots <- .covstatis_relative_roots(fit$roots)
  if (roots[[1L]] == 0) return(0)
  1 / sum(roots)
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
  length(unique(dims)) == 1L && dims[[1L]] > 0L
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
