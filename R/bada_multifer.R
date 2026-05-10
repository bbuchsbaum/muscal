#' Multifer Inference for BaDA
#'
#' Run `multifer` inference for a fitted BaDA model using subject-level
#' class-barycenter blocks as the multiblock data representation.
#'
#' @param x A fitted `bada` object.
#' @param data The original `multidesign` data used to fit `x`.
#' @param targets Inference targets passed to [multifer::infer()].
#' @param B Number of Monte Carlo draws per component-significance rung.
#' @param R Number of bootstrap replicates for stability targets.
#' @param alpha Significance level passed to [multifer::infer()].
#' @param seed Optional random seed.
#' @param strict Logical; passed to [multifer::infer()].
#' @param ... Additional arguments passed to [multifer::infer()].
#'
#' @return A `multifer::infer_result` object.
#' @export
infer_bada <- function(x,
                       data,
                       targets = "default",
                       B = 1000L,
                       R = 500L,
                       alpha = 0.05,
                       seed = NULL,
                       strict = TRUE,
                       ...) {
  if (!inherits(x, "bada")) {
    stop("`x` must be a fitted bada object.", call. = FALSE)
  }
  if (!requireNamespace("multifer", quietly = TRUE)) {
    stop("`infer_bada()` requires the multifer package.", call. = FALSE)
  }
  .bada_require_multifer_contract()
  if (isTRUE((x$rescomp %||% 0L) > 0L)) {
    stop(
      "`infer_bada()` currently supports BaDA fits with `rescomp = 0` only.",
      call. = FALSE
    )
  }

  blocks <- .bada_subject_barycenter_blocks(x, data)
  adapter <- .bada_multifer_adapter()

  multifer::infer(
    adapter = adapter,
    data = blocks,
    geometry = "multiblock",
    relation = "variance",
    targets = targets,
    model = x,
    B = B,
    R = R,
    alpha = alpha,
    seed = seed,
    strict = strict,
    ...
  )
}

.bada_multifer_adapter <- function() {
  .bada_require_multifer_contract()

  multifer::infer_adapter(
    adapter_id = "muscal_bada",
    adapter_version = "0.1.0",
    shape_kinds = "multiblock",
    capabilities = multifer::capability_matrix(
      list(
        geometry = "multiblock",
        relation = "variance",
        targets = c(
          "component_significance",
          "variable_stability",
          "score_stability",
          "subspace_stability"
        )
      )
    ),
    domains = function(x, data = NULL, ...) {
      "global"
    },
    roots = function(x, ...) {
      .bada_multifer_roots(x)
    },
    scores = function(x, domain = NULL, ...) {
      .bada_check_global_domain(domain)
      as.matrix(x$fscores)
    },
    loadings = function(x, domain = NULL, ...) {
      .bada_check_global_domain(domain)
      as.matrix(x$v)
    },
    variable_stat = function(x, data, domain = NULL, k, ...) {
      .bada_check_global_domain(domain)
      .muscal_squared_loading_stat(as.matrix(x$v), k)
    },
    project_scores = function(x, data, domain = NULL, ...) {
      .bada_check_global_domain(domain)
      Xc <- .bada_average_barycenter_blocks(data)
      Xc %*% as.matrix(x$v)
    },
    residualize = function(x, k, data, ...) {
      if (!identical(as.integer(k), 1L)) {
        stop("BaDA multiblock residualization expects `k = 1` on the current residual.",
             call. = FALSE)
      }
      V <- as.matrix(x$v)
      if (ncol(V) < 1L) {
        stop("BaDA fit has no loading columns to residualize.", call. = FALSE)
      }
      v1 <- V[, 1L, drop = FALSE]
      out <- lapply(data, function(block) {
        block - block %*% v1 %*% t(v1)
      })
      names(out) <- names(data)
      out
    },
    refit = function(x, new_data, ...) {
      .bada_fit_from_barycenter_blocks(new_data, reference = x)
    },
    bootstrap_action = function(x, data, design, replicate = NULL, ...) {
      idx <- sample.int(length(data), size = length(data), replace = TRUE)
      boot_data <- data[idx]
      fit <- .bada_fit_from_barycenter_blocks(boot_data, reference = x)
      list(
        fit = fit,
        resample_indices = idx,
        info = list(
          unit = "subject",
          sampled_blocks = names(data)[idx],
          replicate = replicate
        )
      )
    },
    null_action = function(x, data, ...) {
      out <- lapply(data, function(block) {
        block[sample.int(nrow(block)), , drop = FALSE]
      })
      names(out) <- names(data)
      out
    },
    component_stat = function(x, data, k, ...) {
      if (!identical(as.integer(k), 1L)) {
        stop("BaDA multiblock component statistics expect `k = 1` on the current residual.",
             call. = FALSE)
      }
      .bada_leading_root_ratio(data)
    },
    validity_level = "conditional",
    declared_assumptions = c(
      "subject_blocks_exchangeable",
      "class_rows_aligned_within_subject",
      "rescomp_zero"
    ),
    checked_assumptions = list(
      list(
        name = "bada_barycenter_blocks",
        check = function(data, ...) {
          .bada_valid_barycenter_blocks(data)
        },
        detail = "BaDA multifer data must be a list of finite numeric subject barycenter matrices with matching class rows and feature columns."
      )
    )
  )
}

.bada_require_multifer_contract <- function() {
  if (!.bada_multifer_contract_available()) {
    stop(
      paste0(
        "`infer_bada()` requires a multifer build with adapter-owned ",
        "`bootstrap_action()` and `project_scores()` hooks. Install the ",
        "updated multifer from ~/code/multifer or bbuchsbaum/multifer."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.bada_multifer_contract_available <- function() {
  if (!requireNamespace("multifer", quietly = TRUE)) {
    return(FALSE)
  }
  tryCatch({
    multifer::infer_adapter(
      adapter_id = "muscal_bada_contract_probe",
      adapter_version = "0.0.0",
      shape_kinds = "multiblock",
      capabilities = multifer::capability_matrix(
        list(
          geometry = "multiblock",
          relation = "variance",
          targets = "score_stability"
        )
      ),
      loadings = function(x, domain = NULL, ...) matrix(1, 1, 1),
      project_scores = function(x, data, domain = NULL, ...) matrix(1, 1, 1),
      bootstrap_action = function(x, data, design, replicate = NULL, ...) list(fit = x),
      validity_level = "heuristic"
    )
    TRUE
  }, error = function(e) {
    FALSE
  })
}

.bada_subject_barycenter_blocks <- function(x, data) {
  if (is.null(data)) {
    stop("`data` is required to reconstruct BaDA subject barycenter blocks.",
         call. = FALSE)
  }
  if (length(x$subjects) != nrow(data$x)) {
    stop("`data` must have the same number of rows as the fitted BaDA object.",
         call. = FALSE)
  }

  sdat <- split(data, x$subjects)
  subject_names <- names(x$proclist) %||% names(sdat)
  y_name <- rlang::as_name(rlang::quo_get_expr(x$y_var))

  blocks <- lapply(seq_along(sdat), function(i) {
    p <- x$proclist[[i]]
    Xi <- sdat[[i]]$x
    Xout <- multivarious::transform(p, Xi)
    md_i <- multidesign::multidesign(Xout, sdat[[i]]$design)
    grouped <- multidesign::summarize_by(md_i, !!x$y_var)

    labels_i <- if (!is.null(grouped$design[[y_name]])) {
      as.character(grouped$design[[y_name]])
    } else {
      rownames(grouped$x)
    }
    if (is.null(labels_i) || anyNA(match(x$label_set, labels_i))) {
      stop("Every BaDA subject block must contain all fitted class labels.",
           call. = FALSE)
    }

    out <- grouped$x[match(x$label_set, labels_i), , drop = FALSE]
    rownames(out) <- x$label_set
    out
  })

  names(blocks) <- subject_names
  blocks
}

.bada_fit_from_barycenter_blocks <- function(blocks, reference = NULL) {
  if (!isTRUE(.bada_valid_barycenter_blocks(blocks))) {
    stop("Invalid BaDA barycenter blocks.", call. = FALSE)
  }

  Xc <- .bada_average_barycenter_blocks(blocks)
  Xc[!is.finite(Xc)] <- 0
  if (max(abs(Xc)) < 1e-12) {
    Xc <- Xc + diag(1e-8, nrow(Xc), ncol(Xc))
  }

  k_ref <- if (!is.null(reference) && !is.null(reference$v)) {
    ncol(reference$v)
  } else {
    min(nrow(Xc), ncol(Xc))
  }
  ncomp <- min(k_ref, nrow(Xc), ncol(Xc))
  pca_group <- multivarious::pca(
    Xc,
    ncomp = ncomp,
    preproc = multivarious::pass(),
    method = "base"
  )
  v <- as.matrix(pca_group$v)
  s <- do.call(rbind, lapply(blocks, function(block) block %*% v))
  fscores <- Xc %*% v
  sdev <- if (nrow(s) > 1L) {
    apply(s, 2L, stats::sd)
  } else {
    rep(0, ncol(s))
  }
  roots <- .bada_roots_from_matrix(Xc, ncomp)

  structure(
    list(
      v = v,
      s = s,
      sdev = sdev,
      roots = roots,
      fscores = fscores,
      barycenters = Xc,
      label_set = rownames(Xc),
      rescomp = 0L
    ),
    class = c("bada_multifer_fit", "list")
  )
}

.bada_average_barycenter_blocks <- function(blocks) {
  Xc <- Reduce("+", blocks) / length(blocks)
  if (!is.null(rownames(blocks[[1L]]))) {
    rownames(Xc) <- rownames(blocks[[1L]])
  }
  Xc
}

.bada_multifer_roots <- function(x) {
  k <- ncol(as.matrix(x$v))
  if (!is.null(x$barycenters)) {
    return(.bada_roots_from_matrix(as.matrix(x$barycenters), k))
  }
  roots <- x$roots %||% (as.numeric(x$sdev)^2)
  roots <- as.numeric(roots)
  roots[seq_len(min(length(roots), k))]
}

.bada_roots_from_matrix <- function(X, k) {
  d <- svd(X, nu = 0L, nv = 0L)$d
  roots <- d^2
  roots[seq_len(min(length(roots), k))]
}

.bada_leading_root_ratio <- function(blocks) {
  Xc <- .bada_average_barycenter_blocks(blocks)
  roots <- .bada_roots_from_matrix(Xc, min(nrow(Xc), ncol(Xc)))
  denom <- sum(roots)
  if (!is.finite(denom) || denom <= .Machine$double.eps) {
    return(0)
  }
  roots[[1L]] / denom
}

.bada_valid_barycenter_blocks <- function(data) {
  if (!is.list(data) || length(data) < 2L) {
    return(FALSE)
  }
  ok_matrix <- vapply(data, function(block) {
    is.matrix(block) && is.numeric(block) && all(is.finite(block))
  }, logical(1L))
  if (!all(ok_matrix)) {
    return(FALSE)
  }
  nr <- vapply(data, nrow, integer(1L))
  nc <- vapply(data, ncol, integer(1L))
  if (any(nr < 2L) || any(nc < 1L)) {
    return(FALSE)
  }
  length(unique(nr)) == 1L && length(unique(nc)) == 1L
}

.bada_check_global_domain <- function(domain) {
  if (!is.null(domain) && !identical(domain, "global")) {
    stop("BaDA multifer adapter only exposes the `global` domain.",
         call. = FALSE)
  }
  invisible(TRUE)
}
