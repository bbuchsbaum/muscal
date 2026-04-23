#' Aligned Reduced-Rank Regression
#'
#' @description
#' `aligned_rrr()` fits a shared reduced-rank response model for multiblock
#' predictors `X` paired with blockwise multivariate responses `Y`. It is the
#' clean supervised baseline against [response_aligned_mfa()]: each block `k`
#' learns a block-specific regression weight matrix `W_k`, all blocks share a
#' common response basis `B`, and latent scores are defined by `Z_k = X_k W_k`.
#'
#' @details
#' The fitted model minimizes
#' \deqn{
#' \sum_k \beta_k \|R_k^{1/2}(Y_k - X_k W_k B^\top)\|_F^2
#' }
#' plus ridge stabilization on the block regression weights `W_k`. In contrast
#' to [response_aligned_mfa()], there is no `X_k \approx Z_k V_k^\top`
#' reconstruction term and no anchor machinery. This makes `aligned_rrr()` the
#' appropriate baseline when the question is whether explicit predictor
#' reconstruction is buying anything over a shared reduced-rank regression fit.
#'
#' Identifiability is imposed on the score side rather than by constraining
#' `B`: after each iteration, the block score matrices `Z_k = X_k W_k` are
#' orthonormalized in the pooled weighted score space and the shared response
#' basis `B` is rotated accordingly. A deterministic sign convention is then
#' applied componentwise.
#'
#' Responses are preprocessed in a shared pooled response space just as in
#' [response_aligned_mfa()]. Out-of-sample prediction is pure:
#' `x_new -> z_hat -> y_hat`, with no conditional response completion path.
#'
#' @param Y A list of numeric matrices/data.frames. Each element `Y[[k]]` is
#'   `n_k x q` and all blocks must share the same response column dimension.
#' @param X A list of numeric matrices/data.frames. Each element `X[[k]]` is
#'   `n_k x p_k`.
#' @param preproc A `multivarious` preprocessing pipeline (a `pre_processor` or
#'   `prepper`) or a list of them for the `X` blocks.
#' @param response_preproc A single shared `multivarious` preprocessing pipeline
#'   used to define the common response space.
#' @param ncomp Integer latent rank. In v1, `aligned_rrr()` keeps the baseline
#'   honest and requires `ncomp <= ncol(Y[[1]])`.
#' @param block_weight Optional numeric scalar or vector of per-block weights.
#'   Defaults to `1` for every block.
#' @param response_weights Optional rowwise response weights. May be `NULL`, a
#'   single non-negative scalar, a length-`length(X)` vector of blockwise
#'   constants, or a list mirroring `X` with one non-negative weight vector per
#'   block row.
#' @param max_iter Maximum ALS iterations.
#' @param tol Relative convergence tolerance on the objective.
#' @param ridge Non-negative ridge stabilization applied to the block regression
#'   weight updates.
#' @param verbose Logical; if `TRUE`, prints iteration diagnostics.
#' @param ... Unused (reserved for future extensions).
#'
#' @return An object inheriting from `multivarious::multiblock_biprojector` with
#'   additional class `"aligned_rrr"`. The object stores block-specific
#'   regression weights in `W_list`, shared response basis `B`, block-specific
#'   latent scores `Z_list`, rowwise response weights in `response_weights`, and
#'   pooled score normalization weights in `score_weights`. For compatibility
#'   with generic score accessors, the `s` slot stores a concatenated stack of
#'   `Z_list`, `score_index` maps stacked rows back to blocks, and
#'   `score_representation = "stacked_block_scores"` records that this is not a
#'   shared observation-level score table.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(40 * 6), 40, 6)
#' X2 <- matrix(rnorm(36 * 5), 36, 5)
#' B <- qr.Q(qr(matrix(rnorm(3 * 2), 3, 2)))
#' W1 <- matrix(rnorm(6 * 2), 6, 2)
#' W2 <- matrix(rnorm(5 * 2), 5, 2)
#' Y1 <- X1 %*% W1 %*% t(B) + matrix(rnorm(40 * 3, sd = 0.05), 40, 3)
#' Y2 <- X2 %*% W2 %*% t(B) + matrix(rnorm(36 * 3, sd = 0.05), 36, 3)
#'
#' fit <- aligned_rrr(
#'   Y = list(X1 = Y1, X2 = Y2),
#'   X = list(X1 = X1, X2 = X2),
#'   ncomp = 2
#' )
#'
#' pred <- predict(fit, X1[1:5, , drop = FALSE], block = "X1", type = "response")
#' stopifnot(nrow(pred) == 5)
#' }
#' @export
aligned_rrr <- function(Y,
                        X,
                        preproc = multivarious::center(),
                        response_preproc = multivarious::center(),
                        ncomp = 2,
                        block_weight = 1,
                        response_weights = NULL,
                        max_iter = 50,
                        tol = 1e-6,
                        ridge = 1e-8,
                        verbose = FALSE,
                        ...) {
  fit_call <- match.call(expand.dots = FALSE)
  fit_dots <- list(...)

  chk::chk_list(X)
  chk::chk_true(length(X) >= 2L)
  chk::chk_list(Y)
  chk::chk_equal(length(X), length(Y))
  ncomp <- as.integer(ncomp)
  chk::chk_integer(ncomp)
  chk::chk_gte(ncomp, 1)
  max_iter <- as.integer(max_iter)
  chk::chk_integer(max_iter)
  chk::chk_gte(max_iter, 1)
  chk::chk_numeric(tol)
  chk::chk_gte(tol, 0)
  chk::chk_numeric(ridge)
  chk::chk_gte(ridge, 0)
  chk::chk_flag(verbose)

  X <- lapply(X, function(x) {
    chk::chk_true(is.matrix(x) || is.data.frame(x))
    as.matrix(x)
  })
  Y <- lapply(Y, function(y) {
    chk::chk_true(is.matrix(y) || is.data.frame(y))
    as.matrix(y)
  })

  if (is.null(names(X))) names(X) <- paste0("X", seq_along(X))
  if (is.null(names(Y))) {
    names(Y) <- names(X)
  } else {
    Y <- .lmfa_align_named_index_list(Y, names(X), what = "Y")
  }

  for (k in seq_along(X)) {
    chk::chk_equal(nrow(X[[k]]), nrow(Y[[k]]))
  }

  q_dims <- unique(vapply(Y, ncol, integer(1)))
  if (length(q_dims) != 1L) {
    stop("All response blocks in `Y` must have the same number of columns.", call. = FALSE)
  }
  q <- q_dims[[1]]
  if (ncomp > q) {
    stop("`aligned_rrr()` requires `ncomp <= ncol(Y[[1]])` because the shared response basis has orthonormal columns.", call. = FALSE)
  }

  if (is.list(response_preproc) &&
      !inherits(response_preproc, "pre_processor") &&
      !inherits(response_preproc, "prepper")) {
    stop("`response_preproc` must be a single shared preprocessing pipeline, not a list.", call. = FALSE)
  }

  if (is.list(preproc) && !inherits(preproc, "pre_processor") && !inherits(preproc, "prepper")) {
    chk::chk_equal(length(preproc), length(X))
  }
  prep_res <- prepare_block_preprocessors(X, preproc, check_consistent_ncol = FALSE)
  Xp <- prep_res$Xp
  proclist <- prep_res$proclist

  pooled_response_rows <- do.call(rbind, Y)
  response_preproc_fit <- multivarious::fit(
    multivarious::fresh(response_preproc),
    pooled_response_rows
  )
  Yp <- lapply(Y, function(y) multivarious::transform(response_preproc_fit, y))

  block_weight <- .ramfa_parse_block_weights(
    block_weight,
    lengths = vapply(Yp, nrow, integer(1)),
    block_names = names(X),
    what = "block_weight"
  )
  response_weights <- .ramfa_parse_response_weights(
    response_weights = response_weights,
    lengths = vapply(Yp, nrow, integer(1)),
    block_names = names(X)
  )

  if (all(block_weight <= 0) || all(vapply(response_weights, sum, numeric(1)) <= 0)) {
    stop("At least one block must contribute positive response weight.", call. = FALSE)
  }

  B <- .arrr_initialize_response_basis(Yp, K = ncomp)
  W_list <- lapply(Xp, function(Xk) matrix(0, ncol(Xk), ncomp))
  Z_list <- lapply(Xp, function(Xk) matrix(0, nrow(Xk), ncomp))
  names(W_list) <- names(Xp)
  names(Z_list) <- names(Xp)

  objective_trace <- numeric(0)
  prev_obj <- Inf

  for (iter in seq_len(max_iter)) {
    W_list <- lapply(seq_along(Xp), function(k) {
      .arrr_update_block_weights(
        Xk = Xp[[k]],
        Yk = Yp[[k]],
        B = B,
        block_weight = block_weight[[k]],
        response_weight = response_weights[[k]],
        ridge = ridge
      )
    })
    names(W_list) <- names(Xp)

    Z_list <- lapply(seq_along(Xp), function(k) Xp[[k]] %*% W_list[[k]])
    names(Z_list) <- names(Xp)

    ortho <- .arrr_orthonormalize_score_list(
      Z_list = Z_list,
      block_weight = block_weight,
      response_weights = response_weights,
      ridge = ridge
    )
    Z_list <- ortho$Z_list
    W_list <- lapply(W_list, function(W) W %*% ortho$transform)
    B <- B %*% ortho$loading_rot
    B <- .ramfa_update_shared_response_loadings(
      score_list = Z_list,
      response_list = Yp,
      response_alpha = block_weight,
      response_weights = response_weights,
      ridge = ridge
    )

    oriented <- .ramfa_orient_components(
      B = B,
      V_list = W_list,
      Z_list = Z_list
    )
    B <- oriented$B
    W_list <- oriented$V_list
    Z_list <- oriented$Z_list

    obj <- .arrr_objective(
      Y_list = Yp,
      X_list = Xp,
      W_list = W_list,
      B = B,
      block_weight = block_weight,
      response_weights = response_weights,
      ridge = ridge
    )
    objective_trace <- c(objective_trace, obj)
    rel_change <- abs(obj - prev_obj) / (abs(prev_obj) + ridge)
    if (isTRUE(verbose)) {
      message(sprintf("aligned_rrr iter %d: obj=%.6g, rel_change=%.3g", iter, obj, rel_change))
    }
    if (is.finite(prev_obj) && rel_change < tol) break
    prev_obj <- obj
  }

  score_weights <- .arrr_score_weights(block_weight, response_weights)
  score_stack <- .ramfa_stack_scores(Z_list)
  score_index <- score_stack$score_index

  block_indices <- list()
  current <- 1L
  for (k in seq_along(Xp)) {
    block_indices[[k]] <- current:(current + ncol(Xp[[k]]) - 1L)
    current <- current + ncol(Xp[[k]])
  }
  names(block_indices) <- names(Xp)

  proc <- multivarious::concat_pre_processors(proclist, block_indices)
  w_concat <- do.call(rbind, W_list)
  s_stack <- score_stack$s
  sdev <- apply(s_stack, 2, stats::sd)

  cor_x_list <- lapply(seq_along(Xp), function(k) {
    Ck <- tryCatch(stats::cor(Xp[[k]], Z_list[[k]]), error = function(e) NULL)
    if (!is.null(Ck)) Ck[!is.finite(Ck)] <- 0
    Ck
  })
  cor_x_list <- cor_x_list[!vapply(cor_x_list, is.null, logical(1))]
  cor_loadings <- if (length(cor_x_list) > 0) do.call(rbind, cor_x_list) else NULL

  response_fit <- do.call(rbind, lapply(seq_along(Yp), function(k) {
    Y_hat <- Xp[[k]] %*% W_list[[k]] %*% t(B)
    resid <- Yp[[k]] - Y_hat
    w <- response_weights[[k]]
    sse <- block_weight[[k]] * sum(w * rowSums(resid^2))
    tss <- block_weight[[k]] * sum(w * rowSums(Yp[[k]]^2))
    data.frame(
      block = names(Yp)[k],
      n = nrow(Yp[[k]]),
      q = ncol(Yp[[k]]),
      sse = sse,
      tss = tss,
      r2 = if (tss > 0) 1 - sse / tss else NA_real_,
      mean_weight = mean(w)
    )
  }))

  data_refit <- list(
    Y = Y,
    X = X,
    response_weights = response_weights
  )

  fit <- multivarious::multiblock_biprojector(
    v = w_concat,
    s = s_stack,
    sdev = sdev,
    preproc = proc,
    block_indices = block_indices,
    B = B,
    W_list = W_list,
    Z_list = Z_list,
    score_index = score_index,
    block_weight = block_weight,
    response_weights = response_weights,
    score_weights = score_weights,
    objective_trace = objective_trace,
    partial_scores = Z_list,
    cor_loadings = cor_loadings,
    response_fit = response_fit,
    block_preproc = setNames(proclist, names(Xp)),
    response_preproc = response_preproc_fit,
    ridge = ridge,
    score_representation = "stacked_block_scores",
    classes = "aligned_rrr"
  )

  .muscal_attach_fit_contract(
    fit,
    method = "aligned_rrr",
    task = "response_prediction",
    oos_types = c("response", "scores"),
    fit_call = fit_call,
    refit_supported = TRUE,
    prediction_target = "Y",
    refit = .muscal_make_refit_spec(
      data = data_refit,
      fit_fn = function(data) {
        do.call(
          aligned_rrr,
          c(
            list(
              Y = data$Y,
              X = data$X,
              preproc = preproc,
              response_preproc = response_preproc,
              ncomp = ncomp,
              block_weight = block_weight,
              response_weights = data$response_weights,
              max_iter = max_iter,
              tol = tol,
              ridge = ridge,
              verbose = FALSE
            ),
            fit_dots
          )
        )
      },
      bootstrap_fn = .muscal_bootstrap_response_aligned_data,
      permutation_fn = .muscal_permute_response_aligned_data,
      resample_unit = "block_rows"
    )
  )
}

.arrr_initialize_response_basis <- function(Y_list, K) {
  .ramfa_initialize_response_loadings(Y_list, K = K)
}

.arrr_update_block_weights <- function(Xk,
                                       Yk,
                                       B,
                                       block_weight,
                                       response_weight,
                                       ridge = 1e-8) {
  w <- sqrt(as.numeric(block_weight) * as.numeric(response_weight))
  keep <- which(w > 0)
  p <- ncol(Xk)
  K <- ncol(B)
  if (length(keep) == 0L) {
    return(matrix(0, nrow = p, ncol = K))
  }

  Xw <- Xk[keep, , drop = FALSE] * w[keep]
  Yw <- Yk[keep, , drop = FALSE] * w[keep]
  A <- crossprod(Xw)
  C <- crossprod(B)
  D <- crossprod(Xw, Yw) %*% B
  lhs <- kronecker(C, A) + ridge * diag(p * K)
  rhs <- as.vector(D)
  matrix(solve(lhs, rhs), nrow = p, ncol = K)
}

.arrr_orthonormalize_score_list <- function(Z_list,
                                            block_weight,
                                            response_weights,
                                            ridge = 1e-8) {
  K <- ncol(Z_list[[1]])
  Zw <- do.call(rbind, lapply(seq_along(Z_list), function(k) {
    sqrt(as.numeric(block_weight[[k]]) * as.numeric(response_weights[[k]])) * Z_list[[k]]
  }))
  sv <- svd(Zw, nu = K, nv = K)
  d <- sv$d[seq_len(K)]
  safe_d <- pmax(d, sqrt(ridge))
  transform <- sv$v[, seq_len(K), drop = FALSE] %*% diag(1 / safe_d, K, K)
  loading_rot <- sv$v[, seq_len(K), drop = FALSE] %*% diag(safe_d, K, K)
  list(
    Z_list = lapply(Z_list, function(Z) Z %*% transform),
    transform = transform,
    loading_rot = loading_rot
  )
}

.arrr_objective <- function(Y_list,
                            X_list,
                            W_list,
                            B,
                            block_weight,
                            response_weights,
                            ridge = 1e-8) {
  obj <- 0
  for (k in seq_along(X_list)) {
    resid <- Y_list[[k]] - X_list[[k]] %*% W_list[[k]] %*% t(B)
    obj <- obj + as.numeric(block_weight[[k]]) * sum(response_weights[[k]] * rowSums(resid^2))
    obj <- obj + ridge * sum(W_list[[k]]^2)
  }
  obj
}

.arrr_score_weights <- function(block_weight, response_weights) {
  out <- stats::setNames(
    as.numeric(block_weight) * vapply(response_weights, mean, numeric(1)),
    names(block_weight)
  )
  if (all(out <= 0)) {
    stop("At least one block must contribute positive weight to pooled score normalization.", call. = FALSE)
  }
  out
}

.arrr_project_scores <- function(object,
                                 new_data,
                                 block,
                                 preprocess = TRUE) {
  resolved <- .muscal_resolve_block_id(block, names(object$W_list))
  block_name <- resolved$name
  Xnew <- as.matrix(new_data)

  if (isTRUE(preprocess)) {
    Xnew <- multivarious::transform(object$block_preproc[[block_name]], Xnew)
  }

  Wk <- object$W_list[[block_name]]
  if (ncol(Xnew) != nrow(Wk)) {
    stop(
      sprintf("Block '%s' has %d columns, expected %d.", block_name, ncol(Xnew), nrow(Wk)),
      call. = FALSE
    )
  }

  Xnew %*% Wk
}

#' Project New Rows into an Aligned RRR Space
#'
#' @param x A fitted `aligned_rrr` object.
#' @param new_data Numeric matrix/data.frame of new rows for a known predictor
#'   block.
#' @param block Character or integer identifying the block used for projection.
#' @param preprocess Logical; if `TRUE` (default), applies the fitted predictor
#'   preprocessing pipeline before projection.
#' @param ... Unused.
#'
#' @return A numeric matrix of latent scores.
#' @export
project.aligned_rrr <- function(x,
                                new_data,
                                block,
                                preprocess = TRUE,
                                ...) {
  chk::chk_true(inherits(x, "aligned_rrr"))
  .arrr_project_scores(
    object = x,
    new_data = new_data,
    block = block,
    preprocess = preprocess
  )
}

#' Predict from an Aligned RRR Fit
#'
#' @param object A fitted `aligned_rrr` object.
#' @param new_data Numeric matrix/data.frame of new rows from a known predictor
#'   block.
#' @param block Character or integer identifying the block.
#' @param type One of `"response"` or `"scores"`.
#' @param preprocess Logical; if `TRUE` (default), applies the fitted predictor
#'   preprocessing pipeline before projection.
#' @param ... Unused.
#'
#' @return A numeric matrix of predicted responses or latent scores.
#' @export
predict.aligned_rrr <- function(object,
                                new_data,
                                block,
                                type = c("response", "scores"),
                                preprocess = TRUE,
                                ...) {
  type <- match.arg(type)
  scores_new <- .arrr_project_scores(
    object = object,
    new_data = new_data,
    block = block,
    preprocess = preprocess
  )

  if (type == "scores") {
    return(scores_new)
  }

  Yhat_p <- scores_new %*% t(object$B)
  multivarious::inverse_transform(object$response_preproc, Yhat_p)
}
