#' Aligned Multiblock Canonical Correlation Analysis (Aligned MCCA)
#'
#' @description
#' Aligned MCCA extends MAXVAR generalized CCA (GCCA) to the setting where each
#' block has its own row set, linked to a common set of latent reference rows via
#' an integer index vector. No block is privileged: the shared scores live in the
#' reference row space `1..N`, and each block contributes through its row mapping.
#'
#' This is the correlation-based counterpart to [aligned_mfa()], analogous to how
#' [mcca()] relates to [mfa()] when all blocks share rows.
#'
#' @details
#' ## Model (MAXVAR GCCA with row alignment)
#' Let `S` be the shared score matrix for the `N` reference rows. For each block
#' `X_k` with `n_k` rows, `row_index[[k]]` maps its rows to reference rows:
#' \deqn{X_k \; \mathrm{linked\ to}\; S[\mathrm{idx}_k, ]}
#'
#' Aligned MCCA computes `S` as the leading eigenvectors of a weighted sum of
#' block-wise ridge projection operators mapped into the reference space.
#'
#' ## High-dimensional stability
#' The implementation works in observation space and uses ridge-stabilized solves
#' on `n_k × n_k` systems, which remains well-posed when blocks have more variables
#' than rows (`p_k >> n_k`) or are rank deficient after preprocessing.
#'
#' @param X A list of numeric matrices/data.frames. Each element `X[[k]]` is
#'   `n_k × p_k`.
#' @param row_index A list of integer vectors. `row_index[[k]]` has length `n_k`
#'   and maps rows of `X[[k]]` to shared rows in `1..N`. Repeated indices are allowed.
#' @param N Optional integer specifying the number of shared rows. If `NULL`
#'   (default), `N` is inferred as `max(unlist(row_index))`.
#' @param preproc A `multivarious` preprocessing pipeline (a `pre_processor`/`prepper`)
#'   or a list of them. If a list, it must have length `length(X)` and will be
#'   applied to `X` in that order.
#' @param ncomp Integer; number of canonical components to compute.
#' @param ridge Non-negative numeric scalar (or vector of length equal to the
#'   number of blocks) controlling ridge stabilization. The effective ridge used
#'   per block is scaled by the block's overall energy so the default works
#'   across variable scalings.
#' @param block_weights Optional numeric vector of non-negative weights (length
#'   = number of blocks) controlling each block's influence on the compromise.
#' @param use_future Logical; if `TRUE`, block-wise computations are performed via
#'   `furrr::future_map()` when available.
#' @param ... Unused (reserved for future extensions).
#'
#' @return An object inheriting from `multivarious::multiblock_biprojector` with
#'   additional class `"aligned_mcca"`. The object contains reference-space scores
#'   in `s` (`N × ncomp`) and block-wise canonical weights in `canonical_weights`.
#'   Per-block row-level score estimates are available in `partial_scores`.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' N <- 30
#' X1 <- matrix(rnorm(20 * 10), 20, 10)
#' X2 <- matrix(rnorm(15 * 8), 15, 8)
#' idx1 <- sample.int(N, nrow(X1), replace = TRUE)
#' idx2 <- sample.int(N, nrow(X2), replace = TRUE)
#'
#' fit <- aligned_mcca(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), N = N, ncomp = 2)
#' stopifnot(nrow(multivarious::scores(fit)) == N)
#' }
#' @export
aligned_mcca <- function(X,
                         row_index,
                         N = NULL,
                         preproc = multivarious::center(),
                         ncomp = 2,
                         ridge = 1e-6,
                         block_weights = NULL,
                         use_future = FALSE,
                         ...) {
  chk::chk_list(X)
  chk::chk_true(length(X) > 1)
  chk::chk_list(row_index)
  chk::chk_equal(length(X), length(row_index))
  chk::chk_flag(use_future)

  X <- lapply(X, function(x) {
    chk::chk_true(is.matrix(x) || is.data.frame(x))
    as.matrix(x)
  })

  S_blocks <- length(X)
  if (is.null(names(X))) names(X) <- paste0("B", seq_len(S_blocks))
  if (is.null(names(row_index))) names(row_index) <- names(X)

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
    if (any(idx < 1L) || any(idx > N)) stop("row_index values must be in 1..N.", call. = FALSE)
    row_index[[k]] <- idx
  }

  # Preprocess blocks independently
  if (is.list(preproc) && !inherits(preproc, "pre_processor") && !inherits(preproc, "prepper")) {
    chk::chk_equal(length(preproc), length(X))
  }
  prep_res <- prepare_block_preprocessors(X, preproc, check_consistent_ncol = FALSE)
  Xp <- prep_res$Xp
  proclist <- prep_res$proclist

  # Block indices for concatenated feature space
  block_indices <- list()
  current <- 1L
  for (k in seq_along(Xp)) {
    block_indices[[k]] <- current:(current + ncol(Xp[[k]]) - 1L)
    current <- current + ncol(Xp[[k]])
  }
  names(block_indices) <- names(Xp)
  proc <- multivarious::concat_pre_processors(proclist, block_indices)

  # Weights and ridge vectors
  if (is.null(block_weights)) block_weights <- rep(1, S_blocks)
  block_weights <- as.numeric(block_weights)
  chk::chk_true(length(block_weights) == S_blocks)
  chk::chk_true(all(is.finite(block_weights)))
  chk::chk_true(all(block_weights >= 0))

  if (length(ridge) == 1) ridge <- rep(ridge, S_blocks)
  ridge <- as.numeric(ridge)
  chk::chk_true(length(ridge) == S_blocks)
  chk::chk_true(all(is.finite(ridge)))
  chk::chk_true(all(ridge >= 0))

  idx_blocks <- seq_along(Xp)
  map_fun <- if (isTRUE(use_future)) {
    function(.x, .f) {
      furrr::future_map(.x, .f, .options = furrr::furrr_options(seed = TRUE))
    }
  } else {
    function(.x, .f) lapply(.x, .f)
  }

  block_fits <- map_fun(idx_blocks, function(k) {
    .mcca_fit_one_block(
      X = Xp[[k]],
      ridge_mult = ridge[k],
      block_weight = block_weights[k],
      block_name = names(Xp)[k],
      n = nrow(Xp[[k]])
    )
  })

  # Map each block's projector into the reference row space and sum
  H <- matrix(0, nrow = N, ncol = N)
  for (k in seq_along(Xp)) {
    Pk <- block_fits[[k]]$P # n_k x n_k (already weighted)
    idx <- row_index[[k]]

    rs1 <- rowsum(Pk, idx, reorder = FALSE)     # groups x n_k
    rs2 <- rowsum(t(rs1), idx, reorder = FALSE) # groups x groups

    g_rows <- as.integer(rownames(rs2))
    g_cols <- as.integer(rownames(rs1))
    H[g_rows, g_cols] <- H[g_rows, g_cols] + rs2
  }
  H <- (H + t(H)) / 2

  # Cap ncomp
  chk::chk_numeric(ncomp)
  chk::chk_true(length(ncomp) == 1)
  chk::chk_true(is.finite(ncomp))
  chk::chk_gte(ncomp, 1)
  if (abs(ncomp - round(ncomp)) > 1e-8) stop("ncomp must be an integer-like value.", call. = FALSE)
  ncomp_eff <- min(as.integer(round(ncomp)), N)

  eig <- NULL
  if (ncomp_eff < N) {
    eig <- tryCatch(
      RSpectra::eigs_sym(H, k = ncomp_eff, which = "LA"),
      error = function(e) NULL
    )
  }
  if (is.null(eig)) {
    full <- eigen(H, symmetric = TRUE)
    eig <- list(
      values = full$values[seq_len(ncomp_eff)],
      vectors = full$vectors[, seq_len(ncomp_eff), drop = FALSE]
    )
  }

  ord <- order(eig$values, decreasing = TRUE)
  lambda <- eig$values[ord]
  G <- eig$vectors[, ord, drop = FALSE]

  lambda[!is.finite(lambda)] <- 0
  lambda[lambda < 0] <- 0
  sdev <- sqrt(lambda)

  S_scores <- G %*% diag(sdev, nrow = length(sdev), ncol = length(sdev))
  colnames(S_scores) <- paste0("Comp", seq_len(ncol(S_scores)))

  # Build concatenated canonical weights V (feature space) and block partial scores
  p_tot <- sum(vapply(Xp, ncol, integer(1)))
  v_concat <- matrix(0, p_tot, ncomp_eff)

  partial_scores <- vector("list", S_blocks)
  names(partial_scores) <- names(Xp)

  W_list <- vector("list", S_blocks)
  names(W_list) <- names(Xp)

  for (k in seq_along(Xp)) {
    Xk <- Xp[[k]]
    idx <- row_index[[k]]
    Gk <- G[idx, , drop = FALSE] # n_k x K (reference eigenvectors mapped to block rows)

    Rchol <- block_fits[[k]]$chol
    A <- .mcca_chol_solve(Rchol, Gk) # (K + kappa I)^{-1} Gk
    W_raw <- crossprod(Xk, A)        # p_k x K
    W_list[[k]] <- W_raw

    inv_sdev <- ifelse(sdev > 1e-12, 1 / sdev, 0)
    W_weighted <- block_fits[[k]]$block_weight * W_raw
    V_block <- sweep(W_weighted, 2, inv_sdev, `*`)
    v_concat[block_indices[[k]], ] <- V_block
    partial_scores[[k]] <- Xk %*% V_block
  }

  cor_list <- lapply(seq_along(Xp), function(k) {
    idx <- row_index[[k]]
    Ck <- tryCatch(
      suppressWarnings(stats::cor(Xp[[k]], S_scores[idx, , drop = FALSE])),
      error = function(e) NULL
    )
    if (!is.null(Ck)) Ck[!is.finite(Ck)] <- 0
    Ck
  })
  cor_list <- cor_list[!vapply(cor_list, is.null, logical(1))]
  cor_loadings <- if (length(cor_list) > 0) do.call(rbind, cor_list) else NULL

  multivarious::multiblock_biprojector(
    v = v_concat,
    s = S_scores,
    sdev = sdev,
    preproc = proc,
    block_indices = block_indices,
    # aligned_mcca metadata / diagnostics
    lambda = lambda,
    ridge = ridge,
    kappa = vapply(block_fits, `[[`, numeric(1), "kappa"),
    block_weights = block_weights,
    row_index = row_index,
    N = N,
    partial_scores = partial_scores,
    canonical_weights = W_list,
    cor_loadings = cor_loadings,
    names = names(Xp),
    classes = "aligned_mcca"
  )
}

#' Anchored Multiblock Canonical Correlation Analysis (Anchored MCCA)
#'
#' @description
#' Convenience wrapper around [aligned_mcca()] for the common case with a fully
#' observed reference block `Y` (e.g., stimulus features) and one or more linked
#' blocks `X_k` (e.g., brain data). Rows of each `X_k` are linked to rows of `Y`
#' via `row_index[[k]]`.
#'
#' Compared to the legacy name "linked", "anchored" makes explicit that the
#' shared scores live in the row space of `Y` (with `N = nrow(Y)`).
#'
#' @param Y Numeric matrix/data.frame (`N × q`) serving as the reference block.
#' @param X A list of numeric matrices/data.frames. Each element `X[[k]]` is `n_k × p_k`.
#' @param row_index A list of integer vectors. `row_index[[k]]` has length `n_k`
#'   and maps rows of `X[[k]]` to rows of `Y` (values in `1..N`).
#' @param preproc A `multivarious` preprocessing pipeline (a `pre_processor`/`prepper`)
#'   or a list of them. If a list, it must have length `1 + length(X)` and will be
#'   applied to `c(list(Y), X)` in that order.
#' @param ncomp Integer; number of canonical components to compute.
#' @param ridge Non-negative numeric scalar (or vector of length equal to the
#'   number of blocks) controlling ridge stabilization. The effective ridge used
#'   per block is scaled by the block's overall energy so the default works
#'   across variable scalings.
#' @param block_weights Optional numeric vector of non-negative weights. If unnamed
#'   and of length `length(X)`, it is interpreted as weights for `X` blocks and `Y`
#'   is given weight 1. If of length `1 + length(X)`, the first weight corresponds to `Y`.
#' @param use_future Logical; if `TRUE`, block-wise computations are performed via
#'   `furrr::future_map()` when available.
#' @param ... Additional arguments forwarded to [aligned_mcca()] (currently unused).
#'
#' @return An object of class `"anchored_mcca"` (and `"aligned_mcca"`).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' N <- 30
#' Y <- matrix(rnorm(N * 5), N, 5)
#' X1 <- matrix(rnorm(20 * 10), 20, 10)
#' X2 <- matrix(rnorm(15 * 8), 15, 8)
#' idx1 <- sample.int(N, nrow(X1), replace = TRUE)
#' idx2 <- sample.int(N, nrow(X2), replace = TRUE)
#'
#' fit <- anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
#' stopifnot(nrow(multivarious::scores(fit)) == N)
#' }
#' @export
anchored_mcca <- function(Y,
                          X,
                          row_index,
                          preproc = multivarious::center(),
                          ncomp = 2,
                          ridge = 1e-6,
                          block_weights = NULL,
                          use_future = FALSE,
                          ...) {
  chk::chk_true(is.matrix(Y) || is.data.frame(Y))
  Y <- as.matrix(Y)
  chk::chk_list(X)
  chk::chk_true(length(X) >= 1)
  chk::chk_list(row_index)
  chk::chk_equal(length(X), length(row_index))

  X <- lapply(X, function(x) {
    chk::chk_true(is.matrix(x) || is.data.frame(x))
    as.matrix(x)
  })
  if (is.null(names(X))) names(X) <- paste0("X", seq_along(X))
  if (is.null(names(row_index))) names(row_index) <- names(X)

  N <- nrow(Y)
  for (k in seq_along(X)) {
    idx <- as.integer(row_index[[k]])
    chk::chk_equal(length(idx), nrow(X[[k]]))
    if (anyNA(idx)) stop("row_index cannot contain NA.", call. = FALSE)
    if (any(idx < 1L) || any(idx > N)) stop("row_index values must be in 1..nrow(Y).", call. = FALSE)
    row_index[[k]] <- idx
  }

  X_all <- c(list(Y = Y), X)
  idx_all <- c(list(Y = seq_len(N)), row_index)

  # Adapt weights to include Y if user provided X-only weights
  if (!is.null(block_weights)) {
    bw <- as.numeric(block_weights)
    if (!is.null(names(block_weights))) names(bw) <- names(block_weights)
    if (is.null(names(block_weights))) {
      if (length(bw) == length(X)) {
        bw <- c(1, bw)
      } else if (length(bw) != (length(X) + 1L)) {
        stop("block_weights must have length length(X) or 1 + length(X).", call. = FALSE)
      }
      names(bw) <- c("Y", names(X))
    } else {
      nms <- names(block_weights)
      if (!("Y" %in% nms)) {
        # assume weights for X blocks only
        if (!all(names(X) %in% nms)) stop("Named block_weights must include all X block names.", call. = FALSE)
        bw <- c(Y = 1, bw)
      }
      bw <- bw[c("Y", names(X))]
    }
    block_weights <- bw
  }

  fit <- aligned_mcca(
    X = X_all,
    row_index = idx_all,
    N = N,
    preproc = preproc,
    ncomp = ncomp,
    ridge = ridge,
    block_weights = block_weights,
    use_future = use_future,
    ...
  )
  class(fit) <- unique(c("anchored_mcca", class(fit)))
  fit
}
