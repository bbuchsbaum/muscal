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
#'   `n_k × p_k`.
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
#' @param feature_groups Feature prior specification. One of:
#'   * `NULL` (no feature prior),
#'   * `"colnames"` to group features with identical column names across blocks,
#'   * a `data.frame` with columns `block`, `feature`, `group` and optional `weight`.
#' @param feature_lambda Non-negative scalar controlling strength of the feature prior.
#' @param max_iter Maximum number of alternating least-squares iterations.
#' @param tol Relative tolerance on the objective for convergence.
#' @param ridge Non-negative ridge stabilization added to normal equations.
#' @param verbose Logical; if `TRUE`, prints iteration progress.
#' @param ... Unused (reserved for future extensions).
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
                        feature_groups = NULL,
                        feature_lambda = 0,
                        max_iter = 50,
                        tol = 1e-6,
                        ridge = 1e-8,
                        verbose = FALSE,
                        ...) {
  normalization <- match.arg(normalization)

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

  if (is.null(names(X))) names(X) <- paste0("X", seq_along(X))
  if (is.null(names(row_index))) names(row_index) <- names(X)

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
      mdim <- min(dim(mat))
      method <- if (mdim < 3) "base" else "svds"
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
  feature_lambda_eff <- fg$feature_lambda

  # Determine effective K from feature dimensions (K can exceed per-block row counts with ridge)
  min_p <- min(vapply(Xp, ncol, integer(1)))
  K <- min(as.integer(ncomp), as.integer(N), min_p)
  if (K < 1L) stop("ncomp must be at least 1 and <= min(N, min ncol(X)).", call. = FALSE)

  # Initialization: seed S with mapped PCs from one informative block + small random jitter.
  k0 <- which.max(vapply(Xp, function(M) min(nrow(M), ncol(M)), integer(1)))
  kdim0 <- min(K, min(dim(Xp[[k0]])))
  S_init <- matrix(rnorm(N * K), nrow = N, ncol = K)
  if (kdim0 >= 1L) {
    sv <- multivarious::svd_wrapper(Xp[[k0]], ncomp = kdim0, method = if (min(dim(Xp[[k0]])) < 3) "base" else "svds")
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

  # Initialize V_k by ridge regression of X_k on S[idx_k,]
  V_list <- lapply(seq_along(Xp), function(k) {
    Sk <- S[row_index[[k]], , drop = FALSE]
    .lmfa_update_loadings(Sk, Xp[[k]], ridge = ridge)$V
  })
  names(V_list) <- names(Xp)

  objective_trace <- numeric(0)
  prev_obj <- Inf

  for (iter in seq_len(max_iter)) {
    centers <- .lmfa_group_centers(V_list, fg)

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
        ridge = ridge
      )
    }

    centers <- .lmfa_group_centers(V_list, fg)

    # Update S row-by-row (shared scores for the N reference rows)
    S <- .amfa_update_scores(
      X_list = Xp,
      V_list = V_list,
      row_index = row_index,
      alpha_blocks = alpha_blocks,
      ridge = ridge,
      N = N
    )

    # Orthonormalize S and rotate loadings to preserve fitted values
    ortho <- .lmfa_orthonormalize_scores(S)
    S <- ortho$S
    rot <- ortho$R
    V_list <- lapply(V_list, function(V) V %*% t(rot))

    centers <- .lmfa_group_centers(V_list, fg)

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
  proc <- multivarious::concat_pre_processors(proclist, block_indices)
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

  multivarious::multiblock_biprojector(
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
    classes = "aligned_mfa"
  )
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.amfa_update_scores <- function(X_list,
                               V_list,
                               row_index,
                               alpha_blocks,
                               ridge = 1e-8,
                               N = NULL) {
  if (is.null(N)) N <- max(unlist(row_index))
  K <- ncol(V_list[[1]])

  VtV_list <- lapply(V_list, crossprod) # each K x K
  XV_list <- lapply(seq_along(X_list), function(k) X_list[[k]] %*% V_list[[k]]) # n_k x K

  sums_by_i <- lapply(seq_along(X_list), function(k) {
    idx <- row_index[[k]]
    rs <- rowsum(XV_list[[k]], idx, reorder = FALSE)
    out <- matrix(0, nrow = N, ncol = K)
    out[as.integer(rownames(rs)), ] <- rs
    out
  })

  counts_by_i <- lapply(seq_along(X_list), function(k) {
    tabulate(row_index[[k]], nbins = N)
  })

  S_new <- matrix(0, nrow = N, ncol = K)
  I_K <- diag(K)
  for (i in seq_len(N)) {
    A <- matrix(0, nrow = K, ncol = K)
    b <- rep(0, K)
    for (k in seq_along(X_list)) {
      cnt <- counts_by_i[[k]][i]
      if (cnt == 0L) next
      a_k <- alpha_blocks[k]
      A <- A + a_k * cnt * VtV_list[[k]]
      b <- b + a_k * sums_by_i[[k]][i, ]
    }
    A <- A + ridge * I_K
    if (all(abs(b) < 1e-14)) {
      S_new[i, ] <- 0
    } else {
      S_new[i, ] <- solve(A, b)
    }
  }
  S_new
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

  pen <- 0
  if (isTRUE(feature_lambda > 0) && isTRUE(fg$enabled) && !is.null(centers) && length(fg$members) > 0) {
    groups <- names(fg$members)
    if (is.null(groups)) groups <- seq_along(fg$members)
    for (gi in seq_along(fg$members)) {
      gname <- groups[gi]
      cvec <- centers$centers[match(gname, rownames(centers$centers)), , drop = TRUE]
      for (m in fg$members[[gi]]) {
        v <- V_list[[m$block]][m$feature, , drop = TRUE]
        w <- m$weight
        pen <- pen + w * sum((v - cvec)^2)
      }
    }
    pen <- feature_lambda * pen
  }

  err_x + pen
}
