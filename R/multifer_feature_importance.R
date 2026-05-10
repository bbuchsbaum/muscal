#' Resolve a muscal fit to a multifer adapter
#'
#' Returns the `multifer` adapter used to expose feature-importance hooks for
#' fitted `muscal` methods. The returned adapter can be passed directly to
#' `multifer::feature_importance_from_fit()`,
#' `multifer::check_feature_importance_adapter()`, and, when the fit/data expose
#' refitting hooks, `multifer::feature_importance_pvalues()`.
#'
#' @param x A fitted `muscal` model.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `multifer_adapter`.
#' @export
as_multifer_adapter <- function(x, ...) {
  UseMethod("as_multifer_adapter")
}

#' @export
as_multifer_adapter.mcca <- function(x, ...) {
  .mcca_multifer_adapter()
}

#' @export
as_multifer_adapter.aligned_mcca <- function(x, ...) {
  .aligned_mcca_multifer_adapter()
}

#' @export
as_multifer_adapter.anchored_mcca <- function(x, ...) {
  .anchored_mcca_multifer_adapter()
}

#' @export
as_multifer_adapter.aligned_rrr <- function(x, ...) {
  .arrr_multifer_adapter()
}

#' @export
as_multifer_adapter.bada <- function(x, ...) {
  .bada_multifer_adapter()
}

#' @export
as_multifer_adapter.covstatis <- function(x, axes = c("compromise", "interstructure"), ...) {
  .covstatis_multifer_adapter(axes = axes)
}

.muscal_require_multifer_feature_importance <- function(caller = "as_multifer_adapter()") {
  if (!requireNamespace("multifer", quietly = TRUE)) {
    stop(sprintf("`%s` requires the multifer package.", caller), call. = FALSE)
  }
  needed <- c("infer_adapter", "capability_matrix", "feature_importance_from_fit")
  missing <- needed[!vapply(needed, exists, logical(1L), where = asNamespace("multifer"), inherits = FALSE)]
  if (length(missing) > 0L) {
    stop(
      sprintf(
        "`%s` requires an updated multifer with feature-importance support; missing: %s.",
        caller, paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.muscal_multifer_domains_from_blocks <- function(x) {
  nms <- names(x$scaled_loadings_by_block) %||% names(x$block_indices) %||% x$names
  if (is.null(nms) || length(nms) == 0L) {
    stop("This fit does not expose named loading domains.", call. = FALSE)
  }
  as.character(nms)
}

.muscal_multifer_block_loadings <- function(x, domain = NULL) {
  domains <- .muscal_multifer_domains_from_blocks(x)
  domain <- domain %||% domains[[1L]]
  if (!is.character(domain) || length(domain) != 1L || is.na(domain)) {
    stop("`domain` must be a single loading domain.", call. = FALSE)
  }
  if (!(domain %in% domains)) {
    stop(
      sprintf("Unknown loading domain '%s'. Available domains: %s.",
              domain, paste(domains, collapse = ", ")),
      call. = FALSE
    )
  }

  L <- x$scaled_loadings_by_block[[domain]]
  if (is.null(L)) {
    idx <- x$block_indices[[domain]]
    L <- as.matrix(x$scaled_loadings %||% x$cor_loadings %||% x$v)[idx, , drop = FALSE]
  }
  L <- as.matrix(L)
  if (is.null(rownames(L))) {
    rownames(L) <- paste0(domain, "_", seq_len(nrow(L)))
  }
  if (is.null(colnames(L))) {
    colnames(L) <- paste0("Comp", seq_len(ncol(L)))
  }
  L
}

.muscal_squared_loading_stat <- function(loadings, k) {
  L <- as.matrix(loadings)
  k <- as.integer(k)
  if (length(k) == 0L || anyNA(k) || any(k < 1L) || any(k > ncol(L))) {
    stop("`k` must identify available loading columns.", call. = FALSE)
  }
  out <- rowSums(L[, k, drop = FALSE]^2)
  names(out) <- rownames(L) %||% as.character(seq_along(out))
  out
}

.mcca_multifer_roots <- function(x) {
  roots <- x$lambda %||% (as.numeric(x$sdev)^2)
  roots <- as.numeric(roots)
  roots[!is.finite(roots)] <- 0
  pmax(roots, 0)
}

.mcca_multifer_component_stat <- function(data) {
  fit <- if (is.list(data) && !is.null(data$X)) {
    mcca(data$X, preproc = multivarious::pass(), ncomp = data$ncomp %||% 1L,
         ridge = data$ridge %||% 1e-6,
         block_weights = data$block_weights %||% NULL,
         use_future = FALSE)
  } else {
    mcca(data, preproc = multivarious::pass(), ncomp = 1L, ridge = 1e-6)
  }
  roots <- .mcca_multifer_roots(fit)
  denom <- sum(roots)
  if (!is.finite(denom) || denom <= .Machine$double.eps) 0 else roots[[1L]] / denom
}

.mcca_multifer_adapter <- function() {
  .muscal_require_multifer_feature_importance(".mcca_multifer_adapter()")

  multifer::infer_adapter(
    adapter_id = "muscal_mcca",
    adapter_version = "0.1.0",
    shape_kinds = "multiblock",
    capabilities = multifer::capability_matrix(
      list(
        geometry = "multiblock",
        relation = "correlation",
        targets = c("component_significance", "variable_stability", "subspace_stability")
      )
    ),
    domains = function(x, data = NULL, ...) {
      .muscal_multifer_domains_from_blocks(x)
    },
    roots = function(x, ...) {
      .mcca_multifer_roots(x)
    },
    scores = function(x, domain = NULL, ...) {
      multivarious::scores(x)
    },
    loadings = function(x, domain = NULL, ...) {
      .muscal_multifer_block_loadings(x, domain = domain)
    },
    variable_stat = function(x, data, domain = NULL, k, ...) {
      .muscal_squared_loading_stat(.muscal_multifer_block_loadings(x, domain = domain), k)
    },
    refit = function(x, new_data, ...) {
      if (!is.null(x$fit_spec$refit$fit_fn)) {
        return(x$fit_spec$refit$fit_fn(new_data))
      }
      mcca(new_data, preproc = multivarious::pass(), ncomp = ncol(x$v), ridge = x$ridge %||% 1e-6)
    },
    null_action = function(x, data, ...) {
      lapply(data, function(block) as.matrix(block)[sample.int(nrow(block)), , drop = FALSE])
    },
    component_stat = function(x, data, k, ...) {
      if (!identical(as.integer(k), 1L)) {
        stop("MCCA multifer statistics expect `k = 1` on the current residual.",
             call. = FALSE)
      }
      .mcca_multifer_component_stat(data)
    },
    residualize = function(x, k, data, ...) {
      if (!identical(as.integer(k), 1L)) {
        stop("MCCA multifer residualization expects `k = 1` on the current residual.",
             call. = FALSE)
      }
      domains <- .muscal_multifer_domains_from_blocks(x)
      out <- data
      for (d in domains) {
        L <- x$canonical_weights[[d]]
        if (is.null(L)) {
          L <- .muscal_multifer_block_loadings(x, d)
        }
        v <- as.matrix(L)[, 1L, drop = FALSE]
        denom <- as.numeric(crossprod(v))
        if (is.finite(denom) && denom > .Machine$double.eps) {
          v <- v / sqrt(denom)
          out[[d]] <- as.matrix(data[[d]]) - as.matrix(data[[d]]) %*% v %*% t(v)
        }
      }
      out
    },
    validity_level = "conditional",
    declared_assumptions = c(
      "rows_exchangeable_across_blocks",
      "permutation_null_breaks_cross_block_row_alignment"
    )
  )
}

.aligned_mcca_refit_data <- function(fit, null_payload, original_data) {
  if (is.list(null_payload) && !is.null(null_payload$X)) {
    null_payload
  } else {
    original_data
  }
}

.aligned_mcca_multifer_adapter <- function() {
  .muscal_require_multifer_feature_importance(".aligned_mcca_multifer_adapter()")

  multifer::infer_adapter(
    adapter_id = "muscal_aligned_mcca",
    adapter_version = "0.1.0",
    shape_kinds = "adapter",
    capabilities = multifer::capability_matrix(
      list(
        geometry = "adapter",
        relation = "correlation",
        targets = c("component_significance", "variable_stability", "subspace_stability")
      )
    ),
    domains = function(x, data = NULL, ...) {
      .muscal_multifer_domains_from_blocks(x)
    },
    roots = function(x, ...) {
      .mcca_multifer_roots(x)
    },
    scores = function(x, domain = NULL, ...) {
      multivarious::scores(x)
    },
    loadings = function(x, domain = NULL, ...) {
      .muscal_multifer_block_loadings(x, domain = domain)
    },
    variable_stat = function(x, data, domain = NULL, k, ...) {
      .muscal_squared_loading_stat(.muscal_multifer_block_loadings(x, domain = domain), k)
    },
    refit = function(x, new_data, ...) {
      if (!is.null(x$fit_spec$refit$fit_fn)) {
        return(x$fit_spec$refit$fit_fn(new_data))
      }
      aligned_mcca(
        X = new_data$X,
        row_index = new_data$row_index,
        N = new_data$N,
        preproc = multivarious::pass(),
        ncomp = ncol(x$v),
        normalization = x$normalization %||% "MFA",
        alpha = x$alpha_blocks %||% NULL,
        ridge = x$ridge %||% 1e-6,
        use_future = FALSE
      )
    },
    bootstrap_action = function(x, data, design, replicate = NULL, ...) {
      boot <- .muscal_bootstrap_aligned_data(data)
      fit <- if (!is.null(x$fit_spec$refit$fit_fn)) {
        x$fit_spec$refit$fit_fn(boot)
      } else {
        aligned_mcca(
          X = boot$X,
          row_index = boot$row_index,
          N = boot$N,
          preproc = multivarious::pass(),
          ncomp = ncol(x$v),
          normalization = x$normalization %||% "MFA",
          alpha = x$alpha_blocks %||% NULL,
          ridge = x$ridge %||% 1e-6,
          use_future = FALSE
        )
      }
      list(
        fit = fit,
        data = boot,
        info = list(unit = "reference_rows", replicate = replicate)
      )
    },
    null_action = function(x, data, ...) {
      .muscal_permute_aligned_data(data)
    },
    refit_data = .aligned_mcca_refit_data,
    component_stat = function(x, data, k, ...) {
      if (!identical(as.integer(k), 1L)) {
        stop("Aligned MCCA multifer statistics expect `k = 1` on the current residual.",
             call. = FALSE)
      }
      fit <- aligned_mcca(
        X = data$X,
        row_index = data$row_index,
        N = data$N,
        preproc = multivarious::pass(),
        ncomp = data$ncomp %||% 1L,
        normalization = "None",
        ridge = data$ridge %||% 1e-6,
        use_future = FALSE
      )
      roots <- .mcca_multifer_roots(fit)
      denom <- sum(roots)
      if (!is.finite(denom) || denom <= .Machine$double.eps) 0 else roots[[1L]] / denom
    },
    residualize = function(x, k, data, ...) {
      if (!identical(as.integer(k), 1L)) {
        stop("Aligned MCCA multifer residualization expects `k = 1` on the current residual.",
             call. = FALSE)
      }
      out <- data
      for (d in names(out$X)) {
        L <- x$canonical_weights[[d]]
        if (is.null(L)) next
        v <- as.matrix(L)[, 1L, drop = FALSE]
        denom <- as.numeric(crossprod(v))
        if (is.finite(denom) && denom > .Machine$double.eps) {
          v <- v / sqrt(denom)
          out$X[[d]] <- as.matrix(out$X[[d]]) - as.matrix(out$X[[d]]) %*% v %*% t(v)
        }
      }
      out
    },
    validity_level = "conditional",
    declared_assumptions = c(
      "reference_rows_exchangeable",
      "block_rows_linked_by_row_index",
      "permutation_null_breaks_row_index_alignment"
    )
  )
}

.anchored_mcca_multifer_adapter <- function() {
  .muscal_require_multifer_feature_importance(".anchored_mcca_multifer_adapter()")

  adapter <- .aligned_mcca_multifer_adapter()
  adapter$adapter_id <- "muscal_anchored_mcca"
  adapter$refit <- function(x, new_data, ...) {
    if (!is.null(x$fit_spec$refit$fit_fn)) {
      return(x$fit_spec$refit$fit_fn(new_data))
    }
    anchored_mcca(
      Y = new_data$Y,
      X = new_data$X,
      row_index = new_data$row_index,
      preproc = multivarious::pass(),
      ncomp = ncol(x$v),
      normalization = x$normalization %||% "MFA",
      alpha = x$alpha_blocks %||% NULL,
      ridge = x$ridge %||% 1e-6,
      use_future = FALSE
    )
  }
  adapter$bootstrap_action <- function(x, data, design, replicate = NULL, ...) {
    boot <- .muscal_bootstrap_anchored_data(data)
    fit <- if (!is.null(x$fit_spec$refit$fit_fn)) {
      x$fit_spec$refit$fit_fn(boot)
    } else {
      anchored_mcca(
        Y = boot$Y,
        X = boot$X,
        row_index = boot$row_index,
        preproc = multivarious::pass(),
        ncomp = ncol(x$v),
        normalization = x$normalization %||% "MFA",
        alpha = x$alpha_blocks %||% NULL,
        ridge = x$ridge %||% 1e-6,
        use_future = FALSE
      )
    }
    list(
      fit = fit,
      data = boot,
      info = list(unit = "anchor_rows", replicate = replicate)
    )
  }
  adapter$null_action <- function(x, data, ...) {
    .muscal_permute_anchored_data(data)
  }
  adapter$refit_data <- function(fit, null_payload, original_data) {
    if (is.list(null_payload) && !is.null(null_payload$Y) && !is.null(null_payload$X)) {
      null_payload
    } else {
      original_data
    }
  }
  adapter
}
