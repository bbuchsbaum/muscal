#' Aligned Multiple Factor Analysis (Aligned MFA)
#'
#' @description
#' Aligned MFA estimates a shared score matrix for a set of *latent reference rows*
#' when multiple blocks have different row sets. Each block `X[[k]]` is linked to
#' the shared rows via an integer index vector `row_index[[k]]`, but **no single
#' block is privileged** as an anchor.
#'
#' This is the symmetric counterpart to [anchored_mfa()], which uses an explicit
#' reference block `Y` and estimates scores in the row space of `Y`.
#'
#' @details
#' ## Model
#' The fitted model has the form:
#' \deqn{X_k \approx S[\mathrm{idx}_k,] V_k^\top}
#' where `S` is `N × ncomp`, each `V_k` is `p_k × ncomp`, and `idx_k` maps rows of
#' `X_k` to `1..N`. Repeated indices are allowed (e.g., repeated measures), and
#' contributions are aggregated in the score updates.
#'
#' ## Feature similarity prior
#' When `feature_lambda > 0` and `feature_groups` is supplied, Aligned MFA applies
#' the same group-shrinkage penalty as [anchored_mfa()], pulling loading vectors
#' of grouped features toward a shared group center.
#'
#' @param X A list of numeric matrices/data.frames. Each element `X[[k]]` is
#'   `n_k × p_k`. Alternatively a [block_source()], which is read one column
#'   chunk at a time so the blocks never need to be in memory together (see
#'   *Out-of-core fitting*).
#' @param row_index A list of integer vectors. `row_index[[k]]` has length `n_k`
#'   and maps rows of `X[[k]]` to shared rows in `1..N`.
#' @param N Optional integer specifying the number of shared rows. If `NULL`
#'   (default), `N` is inferred as `max(unlist(row_index))`.
#' @param preproc A `multivarious` preprocessing pipeline (a `pre_processor`/`prepper`)
#'   or a list of them. If a list, it must have length `length(X)` and will be
#'   applied to `X` in that order.
#' @param ncomp Integer number of components to extract.
#' @param normalization Block weighting scheme. `"MFA"` uses inverse squared first
#'   singular value per block; `"None"` uses uniform weights; `"custom"` uses `alpha`.
#' @param alpha Optional numeric vector of per-block weights (length `length(X)`),
#'   used when `normalization = "custom"`.
#' @param score_constraint Identification strategy for the shared score matrix.
#'   `"none"` uses the historical unconstrained update followed by QR
#'   normalization inside each ALS iteration. `"orthonormal"` treats
#'   `S transpose S = I` as part of the model and updates `S` with a constrained
#'   majorization/polar step.
#' @param feature_groups Feature prior specification. One of:
#'   * `NULL` (no feature prior),
#'   * `"colnames"` to group features with identical column names across blocks,
#'   * a `data.frame` with columns `block`, `feature`, `group` and optional `weight`.
#' @param feature_lambda Non-negative scalar controlling strength of the feature prior.
#' @param feature_lambda_scale How the feature prior is measured. `"raw"`
#'   (default, historical) penalizes `sum w ||V_k[j, ] - c_g||^2` on the
#'   loadings' own scale, so useful `feature_lambda` values depend on the size
#'   and scale of the data. `"block"` compares loadings in each block's
#'   MFA-normalized units, `sum w ||sqrt(alpha_k) V_k[j, ] - c_g||^2`, which
#'   makes `feature_lambda` unitless (relative to the block-weighted data term)
#'   and the fit invariant to rescaling any block (up to the absolute `ridge`,
#'   which is not scale-free).
#' @param max_iter Maximum number of alternating least-squares iterations.
#' @param tol Relative tolerance on the objective for convergence.
#' @param ridge Non-negative ridge stabilization added to normal equations.
#' @param verbose Logical; if `TRUE`, prints iteration progress.
#' @param use_future Logical; if `TRUE`, block-wise computations that do not
#'   depend on one another may be performed via `furrr::future_map()` when
#'   available. The main alternating-least-squares loop is intrinsically
#'   sequential and is unaffected. Accepted here primarily for interface parity
#'   with [aligned_mcca()] and [anchored_mfa()].
#' @param strategy How the data are accessed. `"auto"` (default) fits in memory
#'   when `X` is a list; for a [block_source()] it uses the Gram path, or the
#'   streaming path when a feature prior is supplied. `"in_memory"` loads every
#'   block (reading a block source in full). `"gram"` and `"stream"` force the
#'   corresponding out-of-core path for either input.
#' @param gram_tol Relative eigenvalue tolerance used to drop numerically null
#'   directions when compressing each block's Gram matrix (Gram path only).
#' @param gram_rank Optional maximum number of Gram eigenpairs kept per block
#'   (one value, or one per block) on the Gram path. `NULL` (default) keeps all
#'   of them, which is exact. A smaller rank makes the fit approximate but
#'   bounds each compressed block at `n_k × gram_rank`; the fraction of each
#'   block's variance left out is returned as `gram_discarded`.
#' @param ... Unused (reserved for future extensions).
#'
#' @section Out-of-core fitting:
#' Without a feature prior, the alternating least-squares fit depends on each
#' preprocessed block `X_k` only through its row Gram matrix
#' `G_k = X_k X_k^T` (`n_k × n_k`): the loadings are ridge regressions on the
#' scores and so lie in the row space of `X_k`, and the MFA weight is the
#' inverse of the largest eigenvalue of `G_k`. The Gram path therefore
#' (1) reads each block once, chunk by chunk, to accumulate its column-local
#' preprocessing statistics and `G_k`; (2) runs the ordinary fit on compressed
#' blocks `R_k` with `R_k R_k^T = G_k` (at most `n_k` columns each); and
#' (3) reads each block once more to compute the feature-space loadings
#' `V_k = X_k^T S_k (S_k^T S_k + (ridge / alpha_k) I)^{-1}` from the final
#' scores. Peak memory is one chunk, one block's Gram matrix (`n_k^2` values,
#' built and compressed one block at a time), the compressed blocks, and the
#' loadings; `gram_rank` caps the compressed blocks for subjects with many
#' rows. The scores match the in-memory fit; with
#' `score_constraint = "none"` they can differ by a rotation within the same
#' subspace because that initialization depends on an SVD sign. The returned
#' loadings are the ridge regressions on the final scores rather than the last
#' alternating update, which differ only by the convergence tolerance.
#'
#' The feature prior (`feature_lambda > 0` with `feature_groups`) couples
#' features across blocks, so the Gram shortcut does not apply. The streaming
#' path instead runs the same alternating least squares as the in-memory fit
#' with one read of the data per iteration: for each chunk it updates the
#' loadings (feature-wise, given the group centers) and accumulates
#' `X_k V_k` for the score step, evaluating the previous iteration's objective
#' in the same read. It holds the loadings (two copies) and the per-row normal
#' equations but never the data, and needs about `iterations + 3` passes.
#' Feature groups given by name need column names on the source
#' (`block_source(colnames = )`).
#'
#' Both out-of-core paths support column-local preprocessing only (`NULL`,
#' `multivarious::pass()`, `center()`, `standardize()`). Resampling works on
#' block-source fits: permutations shuffle `row_index` only, and each bootstrap
#' replicate is a derived source that reads the original one and keeps the
#' resampled rows, so every replicate refit re-reads the data.
#'
#' @return An object inheriting from `multivarious::multiblock_biprojector` with
#'   additional class `"aligned_mfa"`. The object contains global scores in `s`,
#'   concatenated loadings in `v`, and block mappings in `block_indices`.
#'   Additional fields include `V_list`, `row_index`, `alpha_blocks`, and
#'   `objective_trace`.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' N <- 30
#' X1 <- matrix(rnorm(20 * 10), 20, 10)
#' X2 <- matrix(rnorm(15 * 8), 15, 8)
#' idx1 <- sample.int(N, nrow(X1), replace = FALSE)
#' idx2 <- sample.int(N, nrow(X2), replace = FALSE)
#'
#' fit <- aligned_mfa(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
#' stopifnot(nrow(multivarious::scores(fit)) == N)
#' }
#' @export
aligned_mfa <- function(X,
                        row_index,
                        N = NULL,
                        preproc = multivarious::center(),
                        ncomp = 2,
                        normalization = c("MFA", "None", "custom"),
                        alpha = NULL,
                        score_constraint = c("none", "orthonormal"),
                        feature_groups = NULL,
                        feature_lambda = 0,
                        feature_lambda_scale = c("raw", "block"),
                        max_iter = 50,
                        tol = 1e-6,
                        ridge = 1e-8,
                        verbose = FALSE,
                        use_future = FALSE,
                        strategy = c("auto", "in_memory", "gram", "stream"),
                        gram_tol = 1e-10,
                        gram_rank = NULL,
                        ...) {
  fit_call <- match.call(expand.dots = FALSE)
  fit_dots <- list(...)
  normalization <- match.arg(normalization)
  score_constraint <- match.arg(score_constraint)
  strategy <- match.arg(strategy)
  feature_lambda_scale <- match.arg(feature_lambda_scale)
  is_source <- inherits(X, "muscal_block_source")
  has_prior <- !is.null(feature_groups) && isTRUE(feature_lambda > 0)
  if (identical(strategy, "auto")) {
    strategy <- if (!is_source) "in_memory" else if (has_prior) "stream" else "gram"
  }
  if (strategy %in% c("gram", "stream")) {
    refit_args <- list(
      preproc = preproc, ncomp = ncomp, normalization = normalization, alpha = alpha,
      score_constraint = score_constraint, feature_groups = feature_groups,
      feature_lambda = feature_lambda, feature_lambda_scale = feature_lambda_scale,
      max_iter = max_iter, tol = tol, ridge = ridge,
      verbose = FALSE, use_future = use_future, strategy = strategy, gram_tol = gram_tol,
      gram_rank = gram_rank
    )
    args <- list(
      src = as_block_source(X), row_index = row_index, N = N, preproc = preproc,
      ncomp = ncomp, normalization = normalization, alpha = alpha,
      score_constraint = score_constraint, feature_groups = feature_groups,
      feature_lambda = feature_lambda, max_iter = max_iter, tol = tol, ridge = ridge,
      verbose = verbose, fit_call = fit_call, fit_dots = fit_dots, refit_args = refit_args
    )
    if (identical(strategy, "gram")) {
      return(do.call(.amfa_fit_gram,
                     c(args, list(use_future = use_future, gram_tol = gram_tol,
                                  gram_rank = gram_rank)),
                     quote = TRUE))
    }
    return(do.call(.amfa_fit_stream, c(args, list(feature_lambda_scale = feature_lambda_scale)),
                   quote = TRUE))
  }
  if (is_source) {
    src_names <- names(X)
    X <- lapply(seq_along(src_names), function(k) read_block(X, k))
    names(X) <- src_names
  }
  chk::chk_flag(use_future)
  if (isTRUE(use_future) && !requireNamespace("furrr", quietly = TRUE)) {
    stop("use_future = TRUE requires the 'furrr' package.", call. = FALSE)
  }

  chk::chk_list(X)
  chk::chk_true(length(X) >= 2)
  chk::chk_list(row_index)
  chk::chk_equal(length(X), length(row_index))
  ncomp <- as.integer(ncomp)
  chk::chk_integer(ncomp)
  chk::chk_gte(ncomp, 1)
  chk::chk_numeric(feature_lambda)
  chk::chk_gte(feature_lambda, 0)
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
  data_refit <- lapply(X, function(x) x)

  if (is.null(names(X))) names(X) <- paste0("X", seq_along(X))
  if (is.null(names(row_index))) {
    names(row_index) <- names(X)
  } else {
    row_index <- .lmfa_align_named_index_list(row_index, names(X), what = "row_index")
  }

  # Validate row_index and infer N if needed
  infer_N <- is.null(N)
  if (!infer_N) {
    N <- as.integer(N)
    chk::chk_integer(N)
    chk::chk_gte(N, 1)
  } else {
    N <- 0L
  }
  for (k in seq_along(X)) {
    idx <- row_index[[k]]
    chk::chk_true(is.numeric(idx) || is.integer(idx))
    idx <- as.integer(idx)
    chk::chk_equal(length(idx), nrow(X[[k]]))
    if (anyNA(idx)) stop("row_index cannot contain NA.", call. = FALSE)
    if (infer_N) N <- max(c(N, idx), na.rm = TRUE)
    if (any(idx < 1L) || any(idx > N)) {
      stop("row_index values must be in 1..N.", call. = FALSE)
    }
    row_index[[k]] <- idx
  }

  # Preprocess blocks
  if (is.list(preproc) && !inherits(preproc, "pre_processor") && !inherits(preproc, "prepper")) {
    chk::chk_equal(length(preproc), length(X))
  }
  prep_res <- prepare_block_preprocessors(X, preproc, check_consistent_ncol = FALSE)
  Xp <- prep_res$Xp
  proclist <- prep_res$proclist
  fitted_proclist <- .muscal_materialize_block_preprocessors(X, proclist)

  # Block weights (alpha)
  alpha_blocks <- if (normalization == "custom") {
    if (is.null(alpha)) stop("When normalization = 'custom', 'alpha' must be provided.", call. = FALSE)
    alpha <- as.numeric(alpha)
    chk::chk_equal(length(alpha), length(X))
    if (any(!is.finite(alpha)) || any(alpha < 0)) stop("'alpha' must be finite and non-negative.", call. = FALSE)
    alpha
  } else if (normalization == "None") {
    rep(1, length(X))
  } else {
    first_sv <- function(mat) {
      method <- .muscal_svd_method(mat, ncomp = 1)
      tryCatch(
        multivarious::svd_wrapper(mat, ncomp = 1, method = method)$sdev[1],
        error = function(e) multivarious::svd_wrapper(mat, ncomp = 1, method = "base")$sdev[1]
      )
    }
    sapply(Xp, function(M) 1 / (first_sv(M)^2 + ridge))
  }
  names(alpha_blocks) <- names(X)

  # Feature groups (X-only; all blocks are X here)
  fg <- .lmfa_parse_feature_groups(feature_groups, Xp, feature_lambda)
  fg <- .lmfa_scale_feature_groups(fg, alpha_blocks, feature_lambda_scale)
  feature_lambda_eff <- fg$feature_lambda

  # Determine effective K from feature dimensions (K can exceed per-block row counts with ridge)
  K <- min(as.integer(ncomp), as.integer(N))
  if (K < 1L) stop("ncomp must be at least 1 and <= N.", call. = FALSE)

  if (identical(score_constraint, "orthonormal")) {
    S <- .amfa_initialize_scores(Xp, row_index, alpha_blocks, N = N, K = K)
  } else {
    # Historical initialization: informative block seed plus QR normalization.
    k0 <- which.max(vapply(Xp, function(M) min(nrow(M), ncol(M)), integer(1)))
    kdim0 <- min(K, min(dim(Xp[[k0]])))
    S_init <- matrix(rnorm(N * K), nrow = N, ncol = K)
    if (kdim0 >= 1L) {
      sv <- multivarious::svd_wrapper(Xp[[k0]], ncomp = kdim0, method = .muscal_svd_method(Xp[[k0]], ncomp = kdim0))
      scores0 <- sv$u[, seq_len(kdim0), drop = FALSE] %*% diag(sv$sdev[seq_len(kdim0)], kdim0, kdim0)
      idx0 <- row_index[[k0]]
      rs <- rowsum(scores0, idx0, reorder = FALSE)
      out <- matrix(0, nrow = N, ncol = kdim0)
      out[as.integer(rownames(rs)), ] <- rs
      cnt <- tabulate(idx0, nbins = N)
      nz <- cnt > 0
      out[nz, ] <- out[nz, , drop = FALSE] / cnt[nz]
      S_init[, seq_len(kdim0)] <- S_init[, seq_len(kdim0), drop = FALSE] + out
    }
    ortho <- .lmfa_orthonormalize_scores(S_init)
    S <- ortho$S
  }

  # Initialize V_k by ridge regression of X_k on S[idx_k,]
  V_list <- lapply(seq_along(Xp), function(k) {
    Sk <- S[row_index[[k]], , drop = FALSE]
    .lmfa_update_loadings(Sk, Xp[[k]], ridge = ridge)$V
  })
  names(V_list) <- names(Xp)

  objective_trace <- numeric(0)
  prev_obj <- Inf

  for (iter in seq_len(max_iter)) {
    centers <- .lmfa_group_centers(V_list, fg, block_weights = alpha_blocks)

    # Update V_k (loadings) given current S
    for (k in seq_along(Xp)) {
      Sk <- S[row_index[[k]], , drop = FALSE]
      V_list[[k]] <- .lmfa_update_V_block(
        Sk = Sk,
        Xk = Xp[[k]],
        alpha_block = alpha_blocks[k],
        fg_block = fg$by_block[[k]],
        weights_block = fg$weights_by_block[[k]],
        centers = centers,
        feature_lambda = feature_lambda_eff,
        ridge = ridge,
        block_scale = fg$block_scale[k]
      )
    }

    centers <- .lmfa_group_centers(V_list, fg, block_weights = alpha_blocks)

    local <- .amfa_score_system(
      X_list = Xp,
      V_list = V_list,
      row_index = row_index,
      alpha_blocks = alpha_blocks,
      N = N
    )
    step <- .amfa_score_step(S, V_list, local, row_index, alpha_blocks, ridge, N, score_constraint)
    S <- step$S
    if (!is.null(step$R)) {
      V_list <- lapply(V_list, function(V) V %*% t(step$R))
    }

    centers <- .lmfa_group_centers(V_list, fg, block_weights = alpha_blocks)

    obj <- .amfa_objective(
      X_list = Xp,
      S = S,
      V_list = V_list,
      row_index = row_index,
      alpha_blocks = alpha_blocks,
      fg = fg,
      centers = centers,
      feature_lambda = feature_lambda_eff
    )
    objective_trace <- c(objective_trace, obj)

    rel_change <- abs(obj - prev_obj) / (abs(prev_obj) + ridge)
    if (isTRUE(verbose)) {
      message(
        sprintf(
          "aligned_mfa iter %d: obj=%.6g, rel_change=%.3g",
          iter, obj, rel_change
        )
      )
    }
    if (is.finite(prev_obj) && rel_change < tol) break
    prev_obj <- obj
  }

  # Build return object
  block_indices <- list()
  current <- 1L
  for (k in seq_along(Xp)) {
    block_indices[[k]] <- current:(current + ncol(Xp[[k]]) - 1L)
    current <- current + ncol(Xp[[k]])
  }
  names(block_indices) <- names(Xp)
  proc <- multivarious::concat_pre_processors(fitted_proclist, block_indices)
  v_concat <- do.call(rbind, V_list)
  sdev <- apply(S, 2, stats::sd)

  partial_scores <- lapply(seq_along(Xp), function(k) {
    Vk <- V_list[[k]]
    Xt <- Xp[[k]]
    kdim <- ncol(Vk)
    XtV <- Xt %*% Vk
    G <- crossprod(Vk) + ridge * diag(kdim)
    XtV %*% solve(G)
  })
  names(partial_scores) <- names(Xp)

  cor_x_list <- lapply(seq_along(Xp), function(k) {
    idx <- row_index[[k]]
    Ck <- tryCatch(stats::cor(Xp[[k]], S[idx, , drop = FALSE]), error = function(e) NULL)
    if (!is.null(Ck)) Ck[!is.finite(Ck)] <- 0
    Ck
  })
  cor_x_list <- cor_x_list[!vapply(cor_x_list, is.null, logical(1))]
  cor_loadings <- if (length(cor_x_list) > 0) do.call(rbind, cor_x_list) else NULL

  block_fit <- {
    df <- NULL
    for (k in seq_along(Xp)) {
      idx <- row_index[[k]]
      X_hat <- S[idx, , drop = FALSE] %*% t(V_list[[k]])
      sse <- sum((Xp[[k]] - X_hat)^2)
      tss <- sum(Xp[[k]]^2)
      r2 <- if (tss > 0) 1 - sse / tss else NA_real_
      row <- data.frame(
        block = names(Xp)[k],
        n = nrow(Xp[[k]]),
        p = ncol(Xp[[k]]),
        sse = sse,
        tss = tss,
        r2 = r2
      )
      df <- if (is.null(df)) row else rbind(df, row)
    }
    df
  }

  fit <- multivarious::multiblock_biprojector(
    v = v_concat,
    s = S,
    sdev = sdev,
    preproc = proc,
    block_indices = block_indices,
    # aligned_mfa metadata
    V_list = V_list,
    row_index = row_index,
    alpha_blocks = alpha_blocks,
    normalization = normalization,
    feature_groups = feature_groups,
    feature_lambda = feature_lambda_eff,
    objective_trace = objective_trace,
    partial_scores = partial_scores,
    cor_loadings = cor_loadings,
    block_fit = block_fit,
    N = N,
    block_preproc = setNames(fitted_proclist, names(Xp)),
    ridge = ridge,
    classes = "aligned_mfa"
  )

  .muscal_attach_fit_contract(
    fit,
    method = "aligned_mfa",
    task = "row_alignment",
    oos_types = c("scores", "reconstruction"),
    fit_call = fit_call,
    refit_supported = TRUE,
    prediction_target = "blocks",
    refit = .muscal_make_refit_spec(
      data = list(
        X = data_refit,
        row_index = row_index,
        N = N
      ),
      fit_fn = function(data) {
        do.call(
          aligned_mfa,
          c(
            list(
              X = data$X,
              row_index = data$row_index,
              N = data$N,
              preproc = preproc,
              ncomp = ncomp,
              normalization = normalization,
              alpha = alpha,
              score_constraint = score_constraint,
              feature_groups = feature_groups,
              feature_lambda = feature_lambda,
              feature_lambda_scale = feature_lambda_scale,
              max_iter = max_iter,
              tol = tol,
              ridge = ridge,
              verbose = FALSE,
              use_future = use_future
            ),
            fit_dots
          )
        )
      },
      bootstrap_fn = .muscal_bootstrap_aligned_data,
      permutation_fn = .muscal_permute_aligned_data,
      resample_unit = "latent_rows"
    )
  )
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.amfa_initialize_scores <- function(X_list, row_index, alpha_blocks, N, K) {
  H <- matrix(0, nrow = N, ncol = N)
  for (k in seq_along(X_list)) {
    Xk <- X_list[[k]]
    H <- .amfa_init_accumulate(H, Xk %*% t(Xk), row_index[[k]], alpha_blocks[k])
  }
  .amfa_init_from_H(H, K)
}

# Add block k's weighted, row-aggregated Gram matrix to the N x N matrix H
# whose leading eigenvectors initialize the orthonormal score path.
.amfa_init_accumulate <- function(H, Gk, idx, alpha_k) {
  rs1 <- rowsum(alpha_k * Gk, idx, reorder = FALSE)
  rs2 <- rowsum(t(rs1), idx, reorder = FALSE)
  g_rows <- as.integer(rownames(rs2))
  g_cols <- as.integer(rownames(rs1))
  H[g_rows, g_cols] <- H[g_rows, g_cols] + rs2
  H
}

.amfa_init_from_H <- function(H, K) {
  H <- (H + t(H)) / 2
  eig <- eigen(H, symmetric = TRUE)
  S_init <- eig$vectors[, seq_len(K), drop = FALSE]
  .muscal_stiefel_retract(S_init)
}

.amfa_update_scores <- function(X_list,
                               V_list,
                               row_index,
                               alpha_blocks,
                               local = NULL,
                               ridge = 1e-8,
                               N = NULL) {
  if (is.null(local)) {
    local <- .amfa_score_system(
      X_list = X_list,
      V_list = V_list,
      row_index = row_index,
      alpha_blocks = alpha_blocks,
      N = N
    )
  }
  N <- nrow(local$rhs)
  K <- ncol(local$rhs)

  S_new <- matrix(0, nrow = N, ncol = K)
  I_K <- diag(K)
  for (i in seq_len(N)) {
    A <- local$A_list[[i]]
    b <- local$rhs[i, ]
    A <- A + ridge * I_K
    if (all(abs(b) < 1e-14)) {
      S_new[i, ] <- 0
    } else {
      S_new[i, ] <- solve(A, b)
    }
  }
  S_new
}

.amfa_score_system <- function(X_list,
                               V_list,
                               row_index,
                               alpha_blocks,
                               N = NULL,
                               XV_list = NULL) {
  if (is.null(N)) N <- max(unlist(row_index))
  K <- ncol(V_list[[1]])

  VtV_list <- lapply(V_list, crossprod)
  if (is.null(XV_list)) {
    XV_list <- lapply(seq_along(X_list), function(k) X_list[[k]] %*% V_list[[k]])
  }

  sums_by_i <- lapply(seq_along(V_list), function(k) {
    idx <- row_index[[k]]
    rs <- rowsum(XV_list[[k]], idx, reorder = FALSE)
    out <- matrix(0, nrow = N, ncol = K)
    out[as.integer(rownames(rs)), ] <- rs
    out
  })

  counts_by_i <- lapply(seq_along(V_list), function(k) {
    tabulate(row_index[[k]], nbins = N)
  })

  A_list <- vector("list", N)
  rhs <- matrix(0, nrow = N, ncol = K)
  for (i in seq_len(N)) {
    A_i <- matrix(0, nrow = K, ncol = K)
    b_i <- rep(0, K)
    for (k in seq_along(V_list)) {
      cnt <- counts_by_i[[k]][i]
      if (cnt == 0L) next
      a_k <- alpha_blocks[k]
      A_i <- A_i + a_k * cnt * VtV_list[[k]]
      b_i <- b_i + a_k * sums_by_i[[k]][i, ]
    }
    A_list[[i]] <- A_i
    rhs[i, ] <- b_i
  }

  list(A_list = A_list, rhs = rhs)
}

.amfa_score_gradient <- function(S,
                                 X_list,
                                 V_list,
                                 row_index,
                                 alpha_blocks,
                                 ridge = 0,
                                 N = NULL) {
  local <- .amfa_score_system(
    X_list = X_list,
    V_list = V_list,
    row_index = row_index,
    alpha_blocks = alpha_blocks,
    N = N
  )
  grad <- matrix(0, nrow = nrow(S), ncol = ncol(S))
  for (i in seq_len(nrow(S))) {
    grad[i, ] <- 2 * (local$A_list[[i]] %*% S[i, ] - local$rhs[i, ])
  }
  if (ridge > 0) {
    grad <- grad + 2 * ridge * S
  }
  grad
}

.amfa_objective <- function(X_list,
                            S,
                            V_list,
                            row_index,
                            alpha_blocks,
                            fg,
                            centers,
                            feature_lambda) {
  err_x <- 0
  for (k in seq_along(X_list)) {
    Sk <- S[row_index[[k]], , drop = FALSE]
    resid <- X_list[[k]] - Sk %*% t(V_list[[k]])
    err_x <- err_x + alpha_blocks[k] * sum(resid^2)
  }

  err_x + .amfa_prior_penalty(V_list, fg, centers, feature_lambda)
}

.amfa_prior_penalty <- function(V_list, fg, centers, feature_lambda) {
  .lmfa_feature_group_penalty(V_list, fg, centers, feature_lambda)
}

# ---------------------------------------------------------------------------
# Out-of-core paths (see the "Out-of-core fitting" section of ?aligned_mfa)
# ---------------------------------------------------------------------------

# Validate row_index against a block source and infer N, as the in-memory
# path does.
.amfa_source_row_index <- function(src, row_index, N) {
  n_blocks <- length(src$names)
  chk::chk_true(n_blocks >= 2)
  chk::chk_list(row_index)
  chk::chk_equal(length(row_index), n_blocks)
  if (is.null(names(row_index))) {
    names(row_index) <- src$names
  } else {
    row_index <- .lmfa_align_named_index_list(row_index, src$names, what = "row_index")
  }
  infer_N <- is.null(N)
  N <- if (infer_N) 0L else as.integer(N)
  for (k in seq_len(n_blocks)) {
    idx <- row_index[[k]]
    chk::chk_true(is.numeric(idx) || is.integer(idx))
    idx <- as.integer(idx)
    if (length(idx) != src$nrow[[k]]) {
      stop(sprintf("row_index[['%s']] has length %d but the block has %d rows.",
                   src$names[[k]], length(idx), src$nrow[[k]]), call. = FALSE)
    }
    if (anyNA(idx)) stop("row_index cannot contain NA.", call. = FALSE)
    if (infer_N) N <- max(c(N, idx))
    if (any(idx < 1L) || any(idx > N)) {
      stop("row_index values must be in 1..N.", call. = FALSE)
    }
    row_index[[k]] <- idx
  }
  list(row_index = row_index, N = N)
}

# Build an aligned_mfa result from per-block summaries computed out of core:
# loadings V_k, XV_k = Xp_k V_k, XSc_k = Xp_k' (S_k - colMeans(S_k)), and the
# unweighted residual sum of squares sse_k.
.amfa_assemble_source_fit <- function(src, kinds, prep, S, V_list, XV_list, XSc_list, sse,
                                      row_index, alpha_blocks, N, normalization,
                                      feature_groups, feature_lambda, objective_trace,
                                      ridge, strategy, extra, fit_call, refit_args,
                                      fit_dots) {
  n_blocks <- length(src$names)
  K <- ncol(S)
  for (k in seq_len(n_blocks)) {
    rownames(V_list[[k]]) <- src$colnames[[k]]
  }
  names(V_list) <- src$names

  partial_scores <- lapply(seq_len(n_blocks), function(k) {
    XV_list[[k]] %*% solve(crossprod(V_list[[k]]) + ridge * diag(K))
  })
  names(partial_scores) <- src$names

  cor_loadings <- do.call(rbind, lapply(seq_len(n_blocks), function(k) {
    Sk <- S[row_index[[k]], , drop = FALSE]
    s_norm <- sqrt(colSums(sweep(Sk, 2L, colMeans(Sk))^2))
    Ck <- XSc_list[[k]] / outer(sqrt(prep[[k]]$col_ss), s_norm)
    Ck[!is.finite(Ck)] <- 0
    Ck
  }))

  block_fit <- do.call(rbind, lapply(seq_len(n_blocks), function(k) {
    tss <- prep[[k]]$tss
    sse_k <- max(sse[[k]], 0)
    data.frame(
      block = src$names[[k]],
      n = src$nrow[[k]],
      p = src$ncol[[k]],
      sse = sse_k,
      tss = tss,
      r2 = if (tss > 0) 1 - sse_k / tss else NA_real_
    )
  }))

  block_indices <- list()
  current <- 1L
  for (k in seq_len(n_blocks)) {
    block_indices[[k]] <- current:(current + src$ncol[[k]] - 1L)
    current <- current + src$ncol[[k]]
  }
  names(block_indices) <- src$names
  fitted_proclist <- lapply(seq_len(n_blocks), function(k) {
    .muscal_stream_fitted_preproc(kinds[[k]], prep[[k]], src$ncol[[k]])
  })
  names(fitted_proclist) <- src$names

  fit <- do.call(multivarious::multiblock_biprojector, c(
    list(
      v = do.call(rbind, V_list),
      s = S,
      sdev = apply(S, 2, stats::sd),
      preproc = multivarious::concat_pre_processors(fitted_proclist, block_indices),
      block_indices = block_indices,
      V_list = V_list,
      row_index = row_index,
      alpha_blocks = alpha_blocks,
      normalization = normalization,
      feature_groups = feature_groups,
      feature_lambda = feature_lambda,
      objective_trace = objective_trace,
      partial_scores = partial_scores,
      cor_loadings = cor_loadings,
      block_fit = block_fit,
      N = N,
      block_preproc = fitted_proclist,
      ridge = ridge,
      strategy = strategy
    ),
    extra,
    list(classes = "aligned_mfa")
  ))

  .muscal_attach_fit_contract(
    fit,
    method = "aligned_mfa",
    task = "row_alignment",
    oos_types = c("scores", "reconstruction"),
    fit_call = fit_call,
    refit_supported = TRUE,
    prediction_target = "blocks",
    refit = .muscal_make_refit_spec(
      data = list(X = src, row_index = row_index, N = N),
      fit_fn = function(data) {
        do.call(
          aligned_mfa,
          c(list(X = data$X, row_index = data$row_index, N = data$N), refit_args, fit_dots)
        )
      },
      bootstrap_fn = .muscal_bootstrap_aligned_source,
      permutation_fn = .muscal_permute_aligned_data,
      resample_unit = "latent_rows"
    )
  )
}

# Gram path: without a feature prior the fit depends on each block only through
# G_k = Xp_k Xp_k', so fit on compressed blocks and read once more for loadings.
.amfa_fit_gram <- function(src, row_index, N, preproc, ncomp, normalization, alpha,
                           score_constraint, feature_groups, feature_lambda, max_iter,
                           tol, ridge, verbose, use_future, gram_tol, fit_call, fit_dots,
                           refit_args, gram_rank = NULL) {
  if (!is.null(feature_groups) && isTRUE(feature_lambda > 0)) {
    stop(
      "strategy = 'gram' cannot use the feature prior (feature_groups with ",
      "feature_lambda > 0); use strategy = 'stream' or 'auto'.",
      call. = FALSE
    )
  }
  chk::chk_number(gram_tol)
  chk::chk_gte(gram_tol, 0)
  ri <- .amfa_source_row_index(src, row_index, N)
  row_index <- ri$row_index
  n_blocks <- length(src$names)
  if (!is.null(gram_rank)) {
    gram_rank <- as.integer(gram_rank)
    if (length(gram_rank) == 1L) gram_rank <- rep(gram_rank, n_blocks)
    if (length(gram_rank) != n_blocks || anyNA(gram_rank) || any(gram_rank < 1L)) {
      stop("`gram_rank` must be NULL or positive integers (one, or one per block).",
           call. = FALSE)
    }
  }

  # Pass 1: preprocessing statistics, compressing each block as soon as its
  # Gram matrix is complete, so only one n_k x n_k matrix exists at a time.
  kinds <- .muscal_stream_preproc_kinds(preproc, n_blocks)
  R_list <- vector("list", n_blocks)
  lead <- numeric(n_blocks)
  gram <- .muscal_stream_gram(src, kinds, keep_gram = FALSE, on_gram = function(k, Gk) {
    R_list[[k]] <<- .muscal_gram_root(Gk, tol = gram_tol, max_rank = gram_rank[k])
    lead[[k]] <<- max(attr(R_list[[k]], "lead"), 0)
  })
  # MFA weights from the exact leading eigenvalue of each G_k; a randomized
  # first-singular-value solver on the small compressed blocks would make the
  # fit depend on the random seed.
  if (identical(normalization, "MFA")) {
    alpha <- 1 / (lead + ridge)
  }
  discarded <- vapply(R_list, attr, numeric(1L), which = "discarded")
  R_list <- lapply(R_list, function(R) { attr(R, "discarded") <- attr(R, "lead") <- NULL; R })
  names(R_list) <- names(discarded) <- src$names

  inner <- aligned_mfa(
    X = R_list,
    row_index = row_index,
    N = ri$N,
    preproc = NULL,
    ncomp = ncomp,
    normalization = if (identical(normalization, "MFA")) "custom" else normalization,
    alpha = alpha,
    score_constraint = score_constraint,
    max_iter = max_iter,
    tol = tol,
    ridge = ridge,
    verbose = verbose,
    use_future = use_future,
    strategy = "in_memory"
  )
  S <- multivarious::scores(inner)
  alpha_blocks <- inner$alpha_blocks
  K <- ncol(S)

  # Pass 2: loadings V_k = Xp_k' M_k from the final scores, plus Xp_k V_k, the
  # residual cross-term, and the score cross-products for correlation loadings,
  # all read from the data so they stay exact when the Gram matrices were
  # truncated.
  S_blocks <- lapply(row_index, function(idx) S[idx, , drop = FALSE])
  rhs_list <- lapply(seq_len(n_blocks), function(k) {
    Sk <- S_blocks[[k]]
    a <- alpha_blocks[[k]]
    Mk <- if (a > 0) {
      a * Sk %*% solve(a * crossprod(Sk) + ridge * diag(K))
    } else {
      matrix(0, nrow(Sk), K)
    }
    cbind(Mk, sweep(Sk, 2L, colMeans(Sk)), Sk)
  })
  V_list <- lapply(seq_len(n_blocks), function(k) matrix(0, src$ncol[[k]], K))
  XSc_list <- V_list
  XV_list <- lapply(seq_len(n_blocks), function(k) matrix(0, src$nrow[[k]], K))
  tr <- numeric(n_blocks)
  .muscal_stream_pass(src, kinds, gram$prep, function(k, cols, Xp) {
    P <- crossprod(Xp, rhs_list[[k]])
    Vc <- P[, seq_len(K), drop = FALSE]
    V_list[[k]][cols, ] <<- Vc
    XSc_list[[k]][cols, ] <<- P[, K + seq_len(K), drop = FALSE]
    XV_list[[k]] <<- XV_list[[k]] + Xp %*% Vc
    tr[[k]] <<- tr[[k]] + sum(Vc * P[, 2L * K + seq_len(K), drop = FALSE])
  })
  sse <- vapply(seq_len(n_blocks), function(k) {
    gram$prep[[k]]$tss - 2 * tr[[k]] + sum(crossprod(V_list[[k]]) * crossprod(S_blocks[[k]]))
  }, numeric(1L))

  .amfa_assemble_source_fit(
    src = src, kinds = kinds, prep = gram$prep, S = S, V_list = V_list,
    XV_list = XV_list, XSc_list = XSc_list,
    sse = sse, row_index = row_index, alpha_blocks = alpha_blocks, N = ri$N,
    normalization = normalization, feature_groups = NULL, feature_lambda = 0,
    objective_trace = inner$objective_trace, ridge = ridge, strategy = "gram",
    extra = list(
      gram_rank = vapply(R_list, ncol, integer(1L)),
      gram_discarded = discarded
    ),
    fit_call = fit_call, refit_args = refit_args, fit_dots = fit_dots
  )
}

# Streaming path: the in-memory alternating least squares with one read of the
# data per iteration. Each pass, for every chunk, (a) evaluates the previous
# iteration's objective, (b) updates that chunk's loadings (the update is
# feature-wise given the group centers), and (c) accumulates Xp_k V_k for the
# score step. Only loadings and per-row normal equations are kept in memory.
.amfa_fit_stream <- function(src, row_index, N, preproc, ncomp, normalization, alpha,
                             score_constraint, feature_groups, feature_lambda, max_iter,
                             tol, ridge, verbose, fit_call, fit_dots, refit_args,
                             feature_lambda_scale = "raw") {
  ri <- .amfa_source_row_index(src, row_index, N)
  row_index <- ri$row_index
  N <- ri$N
  n_blocks <- length(src$names)
  K <- min(as.integer(ncomp), as.integer(N))
  if (K < 1L) stop("ncomp must be at least 1 and <= N.", call. = FALSE)
  kinds <- .muscal_stream_preproc_kinds(preproc, n_blocks)

  if (normalization == "custom") {
    if (is.null(alpha)) stop("When normalization = 'custom', 'alpha' must be provided.", call. = FALSE)
    alpha <- as.numeric(alpha)
    chk::chk_equal(length(alpha), n_blocks)
    if (any(!is.finite(alpha)) || any(alpha < 0)) {
      stop("'alpha' must be finite and non-negative.", call. = FALSE)
    }
  }
  alpha_blocks <- switch(normalization,
    custom = alpha,
    None = rep(1, n_blocks),
    MFA = numeric(n_blocks)
  )
  names(alpha_blocks) <- src$names

  # Pass 1: preprocessing statistics, block weights, and the initialization,
  # holding one block's Gram matrix at a time.
  k0 <- which.max(pmin(src$nrow, src$ncol))
  kdim0 <- min(K, src$nrow[[k0]], src$ncol[[k0]])
  H <- if (identical(score_constraint, "orthonormal")) matrix(0, N, N) else NULL
  scores0 <- NULL
  gram <- .muscal_stream_gram(src, kinds, keep_gram = FALSE, on_gram = function(k, Gk) {
    if (identical(normalization, "MFA")) {
      lead <- max(eigen(Gk, symmetric = TRUE, only.values = TRUE)$values[[1L]], 0)
      alpha_blocks[[k]] <<- 1 / (lead + ridge)
    }
    if (!is.null(H)) {
      H <<- .amfa_init_accumulate(H, Gk, row_index[[k]], alpha_blocks[[k]])
    } else if (k == k0 && kdim0 >= 1L) {
      e <- eigen(Gk, symmetric = TRUE)
      keep <- seq_len(kdim0)
      scores0 <<- e$vectors[, keep, drop = FALSE] %*%
        diag(sqrt(pmax(e$values[keep], 0)), kdim0, kdim0)
    }
  })
  prep <- gram$prep

  dummy_blocks <- lapply(seq_len(n_blocks), function(k) {
    matrix(0, 0L, src$ncol[[k]], dimnames = list(NULL, src$colnames[[k]]))
  })
  names(dummy_blocks) <- src$names
  fg <- .lmfa_parse_feature_groups(feature_groups, dummy_blocks, feature_lambda)
  fg <- .lmfa_scale_feature_groups(fg, alpha_blocks, feature_lambda_scale)
  feature_lambda_eff <- fg$feature_lambda

  if (!is.null(H)) {
    S <- .amfa_init_from_H(H, K)
    rm(H)
  } else {
    S_init <- matrix(stats::rnorm(N * K), nrow = N, ncol = K)
    if (!is.null(scores0)) {
      idx0 <- row_index[[k0]]
      rs <- rowsum(scores0, idx0, reorder = FALSE)
      out <- matrix(0, nrow = N, ncol = kdim0)
      out[as.integer(rownames(rs)), ] <- rs
      cnt <- tabulate(idx0, nbins = N)
      nz <- cnt > 0
      out[nz, ] <- out[nz, , drop = FALSE] / cnt[nz]
      S_init[, seq_len(kdim0)] <- S_init[, seq_len(kdim0), drop = FALSE] + out
    }
    S <- .lmfa_orthonormalize_scores(S_init)$S
  }

  # Pass 2: initial loadings by ridge regression on the initial scores.
  V_list <- lapply(seq_len(n_blocks), function(k) matrix(0, src$ncol[[k]], K))
  S_blocks <- lapply(row_index, function(idx) S[idx, , drop = FALSE])
  .muscal_stream_pass(src, kinds, prep, function(k, cols, Xp) {
    V_list[[k]][cols, ] <<- .lmfa_update_loadings(S_blocks[[k]], Xp, ridge = ridge)$V
  })

  stream_iteration <- function(S, V_old, centers, update) {
    S_blocks <- lapply(row_index, function(idx) S[idx, , drop = FALSE])
    Sc_blocks <- lapply(S_blocks, function(Sk) sweep(Sk, 2L, colMeans(Sk)))
    tr <- numeric(n_blocks)
    XV_old <- lapply(seq_len(n_blocks), function(k) matrix(0, src$nrow[[k]], K))
    XSc <- lapply(seq_len(n_blocks), function(k) matrix(0, src$ncol[[k]], K))
    V_new <- if (update) V_old else NULL
    XV_new <- if (update) XV_old else NULL
    .muscal_stream_pass(src, kinds, prep, function(k, cols, Xp) {
      Sk <- S_blocks[[k]]
      Vc_old <- V_old[[k]][cols, , drop = FALSE]
      tr[[k]] <<- tr[[k]] + sum(Vc_old * crossprod(Xp, Sk))
      XV_old[[k]] <<- XV_old[[k]] + Xp %*% Vc_old
      XSc[[k]][cols, ] <<- crossprod(Xp, Sc_blocks[[k]])
      if (update) {
        Vc <- .lmfa_update_V_block(
          Sk = Sk,
          Xk = Xp,
          alpha_block = alpha_blocks[[k]],
          fg_block = fg$by_block[[k]][cols],
          weights_block = fg$weights_by_block[[k]][cols],
          centers = centers,
          feature_lambda = feature_lambda_eff,
          ridge = ridge,
          block_scale = fg$block_scale[k]
        )
        V_new[[k]][cols, ] <<- Vc
        XV_new[[k]] <<- XV_new[[k]] + Xp %*% Vc
      }
    })
    sse <- vapply(seq_len(n_blocks), function(k) {
      prep[[k]]$tss - 2 * tr[[k]] + sum(crossprod(V_old[[k]]) * crossprod(S_blocks[[k]]))
    }, numeric(1L))
    list(sse = sse, XV_old = XV_old, XSc = XSc, V_new = V_new, XV_new = XV_new)
  }

  objective_trace <- numeric(0)
  prev_obj <- Inf
  iter <- 0L
  repeat {
    centers <- .lmfa_group_centers(V_list, fg, block_weights = alpha_blocks)
    pass <- stream_iteration(S, V_list, centers, update = iter < max_iter)
    if (iter >= 1L) {
      obj <- unname(sum(alpha_blocks * pass$sse) +
        .amfa_prior_penalty(V_list, fg, centers, feature_lambda_eff))
      objective_trace <- c(objective_trace, obj)
      rel_change <- abs(obj - prev_obj) / (abs(prev_obj) + ridge)
      if (isTRUE(verbose)) {
        message(sprintf("aligned_mfa iter %d: obj=%.6g, rel_change=%.3g", iter, obj, rel_change))
      }
      if ((is.finite(prev_obj) && rel_change < tol) || iter >= max_iter) break
      prev_obj <- obj
    }
    V_list <- pass$V_new
    local <- .amfa_score_system(
      X_list = NULL,
      V_list = V_list,
      row_index = row_index,
      alpha_blocks = alpha_blocks,
      N = N,
      XV_list = pass$XV_new
    )
    step <- .amfa_score_step(S, V_list, local, row_index, alpha_blocks, ridge, N, score_constraint)
    S <- step$S
    if (!is.null(step$R)) {
      V_list <- lapply(V_list, function(V) V %*% t(step$R))
    }
    iter <- iter + 1L
  }

  .amfa_assemble_source_fit(
    src = src, kinds = kinds, prep = prep, S = S, V_list = V_list,
    XV_list = pass$XV_old, XSc_list = pass$XSc, sse = pass$sse,
    row_index = row_index, alpha_blocks = alpha_blocks, N = N,
    normalization = normalization, feature_groups = feature_groups,
    feature_lambda = feature_lambda_eff, objective_trace = objective_trace,
    ridge = ridge, strategy = "stream", extra = list(passes = iter + 3L),
    fit_call = fit_call, refit_args = refit_args, fit_dots = fit_dots
  )
}

# One score update given the per-row normal equations in `local`. With
# score_constraint = "none" the unconstrained solution is QR-normalized and the
# triangular factor is returned so callers can counter-rotate the loadings.
.amfa_score_step <- function(S, V_list, local, row_index, alpha_blocks, ridge, N,
                             score_constraint) {
  if (identical(score_constraint, "orthonormal")) {
    S_warm <- .amfa_update_scores(
      X_list = NULL,
      V_list = V_list,
      row_index = row_index,
      alpha_blocks = alpha_blocks,
      local = local,
      ridge = ridge,
      N = N
    )
    S_start <- .muscal_stiefel_retract(S)
    S_warm <- .muscal_stiefel_retract(S_warm)
    obj_start <- .muscal_score_objective_from_system(S_start, local$A_list, local$rhs, ridge = ridge)
    obj_warm <- .muscal_score_objective_from_system(S_warm, local$A_list, local$rhs, ridge = ridge)
    if (obj_warm <= obj_start) {
      S_start <- S_warm
    }
    s_opt <- .muscal_stiefel_mm(
      S = S_start,
      A_list = local$A_list,
      rhs = local$rhs,
      ridge = ridge
    )
    return(list(S = s_opt$S, R = NULL))
  }
  S_new <- .amfa_update_scores(
    X_list = NULL,
    V_list = V_list,
    row_index = row_index,
    alpha_blocks = alpha_blocks,
    local = local,
    ridge = ridge,
    N = N
  )
  ortho <- .lmfa_orthonormalize_scores(S_new)
  list(S = ortho$S, R = ortho$R)
}
