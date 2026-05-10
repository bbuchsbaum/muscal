#' Multifer Inference for Aligned RRR
#'
#' Run `multifer` inference for a fitted [aligned_rrr()] model. The adapter uses
#' multifer's opaque predictive geometry so the model can keep its native
#' multiblock predictor/response data shape.
#'
#' @param x A fitted `aligned_rrr` object.
#' @param data Optional original data list with elements `X` and `Y`. When
#'   `NULL`, the data stored in `x$fit_spec$refit$data` are used.
#' @param targets Inference targets passed to [multifer::infer()]. Defaults to
#'   all targets declared by the adapter: component significance plus loading,
#'   score, and subspace stability.
#' @param B Number of Monte Carlo draws per component-significance rung.
#' @param R Number of bootstrap replicates for stability targets.
#' @param alpha Significance level passed to [multifer::infer()].
#' @param seed Optional random seed.
#' @param strict Logical; passed to [multifer::infer()].
#' @param ... Additional arguments passed to [multifer::infer()].
#'
#' @return A `multifer::infer_result` object, or a `multifer` bundle when
#'   `return_bundle = TRUE` is passed through `...`.
#' @export
infer_aligned_rrr <- function(x,
                              data = NULL,
                              targets = "default",
                              B = 1000L,
                              R = 500L,
                              alpha = 0.05,
                              seed = NULL,
                              strict = TRUE,
                              ...) {
  if (!inherits(x, "aligned_rrr")) {
    stop("`x` must be a fitted aligned_rrr object.", call. = FALSE)
  }
  if (!requireNamespace("multifer", quietly = TRUE)) {
    stop("`infer_aligned_rrr()` requires the multifer package.", call. = FALSE)
  }
  .arrr_require_multifer_contract()

  payload <- .arrr_multifer_payload(x, data = data)
  adapter <- .arrr_multifer_adapter()
  model <- .arrr_multifer_prepare_fit(x, payload)

  multifer::infer(
    adapter = adapter,
    data = payload,
    geometry = "adapter",
    relation = "predictive",
    targets = targets,
    model = model,
    B = B,
    R = R,
    alpha = alpha,
    seed = seed,
    strict = strict,
    ...
  )
}

.arrr_multifer_adapter <- function() {
  .arrr_require_multifer_contract()

  multifer::infer_adapter(
    adapter_id = "muscal_aligned_rrr",
    adapter_version = "0.1.0",
    shape_kinds = "adapter",
    capabilities = multifer::capability_matrix(
      list(
        geometry = "adapter",
        relation = "predictive",
        targets = c(
          "component_significance",
          "variable_stability",
          "score_stability",
          "subspace_stability"
        )
      )
    ),
    domains = function(x, data = NULL, ...) {
      c("predictor", "response")
    },
    roots = function(x, ...) {
      .arrr_predictive_roots(x, x$multifer_payload)
    },
    scores = function(x, domain = NULL, ...) {
      domain <- .arrr_check_multifer_domain(domain)
      if (identical(domain, "predictor")) {
        return(.ramfa_stack_scores(x$Z_list)$s)
      }
      if (!is.null(x$multifer_response_scores)) {
        return(x$multifer_response_scores)
      }
      payload <- x$multifer_payload
      if (is.null(payload)) {
        stop("Response scores require the aligned_rrr multifer payload.",
             call. = FALSE)
      }
      .arrr_project_response_scores(x, payload)
    },
    loadings = function(x, domain = NULL, ...) {
      domain <- .arrr_check_multifer_domain(domain)
      if (identical(domain, "predictor")) {
        return(as.matrix(x$v))
      }
      as.matrix(x$B)
    },
    variable_stat = function(x, data, domain = NULL, k, ...) {
      domain <- .arrr_check_multifer_domain(domain)
      L <- if (identical(domain, "predictor")) as.matrix(x$v) else as.matrix(x$B)
      .muscal_squared_loading_stat(L, k)
    },
    project_scores = function(x, data, domain = NULL, ...) {
      domain <- .arrr_check_multifer_domain(domain)
      if (identical(domain, "predictor")) {
        return(.arrr_project_predictor_scores(x, data))
      }
      .arrr_project_response_scores(x, data)
    },
    refit = function(x, new_data, ...) {
      .arrr_fit_from_multifer_payload(new_data)
    },
    bootstrap_action = function(x, data, design, replicate = NULL, ...) {
      boot <- .arrr_bootstrap_multifer_payload(data)
      fit <- .arrr_fit_from_multifer_payload(boot$data)
      list(
        fit = fit,
        data = boot$data,
        resample_indices = boot$indices,
        info = list(
          unit = "block_rows",
          replicate = replicate
        )
      )
    },
    null_action = function(x, data, ...) {
      .arrr_permute_multifer_payload(data)
    },
    component_stat = function(x, data, k, split = NULL, ...) {
      if (!identical(as.integer(k), 1L)) {
        stop("Aligned RRR multifer statistics expect `k = 1` on the current residual.",
             call. = FALSE)
      }
      .arrr_leading_predictive_gain(x, data)
    },
    residualize = function(x, k, data, ...) {
      if (!identical(as.integer(k), 1L)) {
        stop("Aligned RRR multifer residualization expects `k = 1` on the current residual.",
             call. = FALSE)
      }
      .arrr_residualize_multifer_payload(x, data)
    },
    validity_level = "conditional",
    declared_assumptions = c(
      "block_rows_exchangeable_within_block",
      "predictor_response_rows_aligned_within_block",
      "permutation_null_breaks_within_block_x_y_association"
    ),
    checked_assumptions = list(
      list(
        name = "aligned_rrr_payload",
        check = function(data, ...) {
          .arrr_valid_multifer_payload(data)
        },
        detail = "Aligned RRR multifer data must be a finite numeric payload with named X/Y block lists, matching block row counts, common response dimension, and valid row weights."
      )
    )
  )
}

.arrr_require_multifer_contract <- function() {
  if (!.arrr_multifer_contract_available()) {
    stop(
      paste0(
        "`infer_aligned_rrr()` requires a multifer build with ",
        "`geometry = \"adapter\"`, `relation = \"predictive\"`, ",
        "and adapter-owned bootstrap/project-score hooks. Install the ",
        "updated multifer from ~/code/multifer or bbuchsbaum/multifer."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.arrr_multifer_contract_available <- function() {
  if (!requireNamespace("multifer", quietly = TRUE)) {
    return(FALSE)
  }
  tryCatch({
    multifer::infer_adapter(
      adapter_id = "muscal_aligned_rrr_contract_probe",
      adapter_version = "0.0.0",
      shape_kinds = "adapter",
      capabilities = multifer::capability_matrix(
        list(
          geometry = "adapter",
          relation = "predictive",
          targets = c("component_significance", "score_stability")
        )
      ),
      roots = function(x, ...) 1,
      loadings = function(x, domain = NULL, ...) matrix(1, 1, 1),
      project_scores = function(x, data, domain = NULL, ...) matrix(1, 1, 1),
      refit = function(x, new_data, ...) list(),
      bootstrap_action = function(x, data, design, replicate = NULL, ...) list(fit = x),
      null_action = function(x, data, ...) data,
      component_stat = function(x, data, k, split = NULL, ...) 1,
      residualize = function(x, k, data, ...) data,
      validity_level = "heuristic"
    )
    TRUE
  }, error = function(e) {
    FALSE
  })
}

.arrr_multifer_payload <- function(x, data = NULL) {
  data <- data %||% x$fit_spec$refit$data
  if (is.null(data) || is.null(data$X) || is.null(data$Y)) {
    stop("`data` must be a list with `X` and `Y`, or `x` must store refit data.",
         call. = FALSE)
  }

  X <- lapply(data$X, as.matrix)
  Y <- lapply(data$Y, as.matrix)
  block_names <- names(x$W_list) %||% names(X) %||% paste0("X", seq_along(X))
  if (is.null(names(X))) names(X) <- block_names
  if (is.null(names(Y))) names(Y) <- names(X)
  X <- .lmfa_align_named_index_list(X, block_names, what = "X")
  Y <- .lmfa_align_named_index_list(Y, block_names, what = "Y")

  Xp <- lapply(block_names, function(nm) {
    multivarious::transform(x$block_preproc[[nm]], X[[nm]])
  })
  Yp <- lapply(Y, function(y) multivarious::transform(x$response_preproc, y))
  names(Xp) <- block_names
  names(Yp) <- block_names

  response_weights <- data$response_weights %||% x$response_weights
  response_weights <- .muscal_expand_block_row_weights(
    weights = response_weights,
    lengths = vapply(Yp, nrow, integer(1)),
    block_names = block_names,
    what = "response_weights"
  )
  if (is.null(response_weights)) {
    response_weights <- lapply(Yp, function(y) rep(1, nrow(y)))
    names(response_weights) <- block_names
  }

  block_weight <- as.numeric(x$block_weight)
  names(block_weight) <- names(x$block_weight) %||% block_names
  block_weight <- block_weight[block_names]

  out <- list(
    X = Xp,
    Y = Yp,
    response_weights = response_weights,
    block_weight = block_weight,
    ncomp = ncol(x$B),
    ridge = x$ridge %||% 1e-8,
    max_iter = x$max_iter %||% max(50L, length(x$objective_trace) * 2L),
    tol = x$tol %||% 1e-8
  )

  if (!isTRUE(.arrr_valid_multifer_payload(out))) {
    stop("Invalid aligned_rrr multifer payload.", call. = FALSE)
  }
  out
}

.arrr_fit_from_multifer_payload <- function(data) {
  if (!isTRUE(.arrr_valid_multifer_payload(data))) {
    stop("Invalid aligned_rrr multifer payload.", call. = FALSE)
  }
  fit <- aligned_rrr(
    Y = data$Y,
    X = data$X,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    ncomp = data$ncomp,
    block_weight = data$block_weight,
    response_weights = data$response_weights,
    ridge = data$ridge,
    max_iter = data$max_iter,
    tol = data$tol
  )
  .arrr_multifer_prepare_fit(fit, data)
}

.arrr_multifer_prepare_fit <- function(fit, data) {
  fit$multifer_payload <- data
  fit$multifer_response_scores <- .arrr_project_response_scores(fit, data)
  fit
}

.arrr_valid_multifer_payload <- function(data) {
  if (!is.list(data) || is.null(data$X) || is.null(data$Y) ||
      !is.list(data$X) || !is.list(data$Y) ||
      length(data$X) < 1L || length(data$X) != length(data$Y)) {
    return(FALSE)
  }
  block_names <- names(data$X)
  if (is.null(block_names) || any(!nzchar(block_names)) ||
      !identical(names(data$Y), block_names)) {
    return(FALSE)
  }
  ok <- vapply(block_names, function(nm) {
    X <- data$X[[nm]]
    Y <- data$Y[[nm]]
    is.matrix(X) && is.numeric(X) &&
      is.matrix(Y) && is.numeric(Y) &&
      nrow(X) == nrow(Y) &&
      nrow(X) >= 2L &&
      ncol(X) >= 1L &&
      ncol(Y) >= 1L &&
      all(is.finite(X)) &&
      all(is.finite(Y))
  }, logical(1L))
  if (!all(ok)) {
    return(FALSE)
  }
  if (length(unique(vapply(data$Y, ncol, integer(1)))) != 1L) {
    return(FALSE)
  }
  if (!is.list(data$response_weights) ||
      !identical(names(data$response_weights), block_names)) {
    return(FALSE)
  }
  weights_ok <- vapply(block_names, function(nm) {
    w <- data$response_weights[[nm]]
    is.numeric(w) &&
      length(w) == nrow(data$Y[[nm]]) &&
      all(is.finite(w)) &&
      all(w >= 0)
  }, logical(1L))
  if (!all(weights_ok)) {
    return(FALSE)
  }
  is.numeric(data$block_weight) &&
    length(data$block_weight) == length(block_names) &&
    identical(names(data$block_weight), block_names) &&
    all(is.finite(data$block_weight)) &&
    all(data$block_weight >= 0) &&
    any(data$block_weight > 0) &&
    is.numeric(data$ncomp) &&
    length(data$ncomp) == 1L &&
    data$ncomp >= 1L &&
    data$ncomp <= ncol(data$Y[[1L]]) &&
    is.numeric(data$ridge) &&
    length(data$ridge) == 1L &&
    is.finite(data$ridge) &&
    data$ridge >= 0 &&
    is.numeric(data$max_iter) &&
    length(data$max_iter) == 1L &&
    is.finite(data$max_iter) &&
    data$max_iter >= 1L &&
    is.numeric(data$tol) &&
    length(data$tol) == 1L &&
    is.finite(data$tol) &&
    data$tol >= 0
}

.arrr_check_multifer_domain <- function(domain) {
  domain <- domain %||% "predictor"
  domain <- as.character(domain[[1L]])
  if (!(domain %in% c("predictor", "response"))) {
    stop("Aligned RRR multifer domains are `predictor` and `response`.",
         call. = FALSE)
  }
  domain
}

.arrr_project_predictor_scores <- function(fit, data) {
  block_names <- names(data$X)
  out <- lapply(block_names, function(nm) {
    data$X[[nm]] %*% fit$W_list[[nm]]
  })
  names(out) <- block_names
  do.call(rbind, out)
}

.arrr_project_response_scores <- function(fit, data) {
  block_names <- names(data$Y)
  out <- lapply(block_names, function(nm) {
    data$Y[[nm]] %*% fit$B
  })
  names(out) <- block_names
  do.call(rbind, out)
}

.arrr_component_prediction <- function(fit, k = 1L) {
  k <- min(as.integer(k), ncol(fit$B))
  if (k <= 0L) {
    return(lapply(fit$Z_list, function(Z) matrix(0, nrow(Z), nrow(fit$B))))
  }
  Bk <- fit$B[, seq_len(k), drop = FALSE]
  lapply(fit$Z_list, function(Z) {
    Z[, seq_len(k), drop = FALSE] %*% t(Bk)
  })
}

.arrr_weighted_sse <- function(data, prediction = NULL) {
  block_names <- names(data$Y)
  sum(vapply(block_names, function(nm) {
    pred <- if (is.null(prediction)) {
      matrix(0, nrow(data$Y[[nm]]), ncol(data$Y[[nm]]))
    } else {
      prediction[[nm]]
    }
    resid <- data$Y[[nm]] - pred
    as.numeric(data$block_weight[[nm]]) *
      sum(data$response_weights[[nm]] * rowSums(resid^2))
  }, numeric(1L)))
}

.arrr_predictive_roots <- function(fit, data) {
  denom <- .arrr_weighted_sse(data)
  if (!is.finite(denom) || denom <= .Machine$double.eps) {
    return(rep(0, ncol(fit$B)))
  }
  sse_prev <- denom
  out <- numeric(ncol(fit$B))
  for (k in seq_len(ncol(fit$B))) {
    pred <- .arrr_component_prediction(fit, k = k)
    sse_k <- .arrr_weighted_sse(data, pred)
    out[[k]] <- max(0, (sse_prev - sse_k) / denom)
    sse_prev <- sse_k
  }
  cummin(out)
}

.arrr_leading_predictive_gain <- function(fit, data) {
  denom <- .arrr_weighted_sse(data)
  if (!is.finite(denom) || denom <= .Machine$double.eps) {
    return(0)
  }
  pred <- .arrr_component_prediction(fit, k = 1L)
  max(0, (denom - .arrr_weighted_sse(data, pred)) / denom)
}

.arrr_residualize_multifer_payload <- function(fit, data) {
  out <- data
  pred <- .arrr_component_prediction(fit, k = 1L)
  block_names <- names(data$X)
  for (nm in block_names) {
    z <- fit$Z_list[[nm]][, 1L, drop = FALSE]
    w <- as.numeric(data$block_weight[[nm]]) * data$response_weights[[nm]]
    denom <- sum(w * z[, 1L]^2)
    if (is.finite(denom) && denom > .Machine$double.eps) {
      coef <- crossprod(z * w, data$X[[nm]]) / denom
      out$X[[nm]] <- data$X[[nm]] - z %*% coef
    }
    out$Y[[nm]] <- data$Y[[nm]] - pred[[nm]]
  }
  out
}

.arrr_permute_multifer_payload <- function(data) {
  out <- data
  for (nm in names(data$Y)) {
    perm <- sample.int(nrow(data$Y[[nm]]))
    out$Y[[nm]] <- data$Y[[nm]][perm, , drop = FALSE]
    out$response_weights[[nm]] <- data$response_weights[[nm]][perm]
  }
  out
}

.arrr_bootstrap_multifer_payload <- function(data) {
  out <- data
  indices <- vector("list", length(data$X))
  names(indices) <- names(data$X)
  for (nm in names(data$X)) {
    idx <- sample.int(nrow(data$X[[nm]]), size = nrow(data$X[[nm]]), replace = TRUE)
    indices[[nm]] <- idx
    out$X[[nm]] <- data$X[[nm]][idx, , drop = FALSE]
    out$Y[[nm]] <- data$Y[[nm]][idx, , drop = FALSE]
    out$response_weights[[nm]] <- data$response_weights[[nm]][idx]
  }
  list(data = out, indices = indices)
}
