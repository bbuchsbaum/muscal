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
#' Aligned MCCA computes `S` as the leading eigenvectors of a weighted sum of
#' block-wise ridge projection operators lifted into the reference space:
#' \deqn{H = \sum_b alpha_b M_b}
#' with `M_b = lift(P_b)` and `P_b = X_b (X_b^T X_b + kappa_b I)^{-1} X_b^T`.
#'
#' ## Block-weighting schemes
#' The block-weight vector `alpha_blocks` is formed according to `normalization`:
#' * `"MFA"` (default) — `alpha_b = 1 / (sigma_1(X_b)^2)` (the CCA analogue of
#'   the MFA normalisation; prevents a block from dominating purely by scale).
#' * `"balanced"` — per-component projected gradient ascent on the geometric
#'   mean criterion `Σ_b β_b log(g^T M_b g)` so every block contributes
#'   meaningfully to every component (uses `target = alpha` if supplied,
#'   else uniform).
#' * `"None"` — uniform weights `alpha_b = 1`.
#' * `"custom"` — use `alpha` as supplied (length `length(X)`).
#'
#' ## High-dimensional stability (p > n)
#' The implementation works in observation space and uses ridge-stabilised
#' solves on `n_b × n_b` systems, so `p_b >> n_b` blocks are handled without
#' special casing. Ridge stabilisation (`kappa_b = ridge_b * mean(diag(K_b))`,
#' with an automatic ridge floor if `ridge = 0` yields a singular system) keeps
#' the Cholesky well-posed.
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
#' @param normalization Block weighting scheme. One of `"MFA"` (default),
#'   `"balanced"`, `"None"`, or `"custom"`. See "Block-weighting schemes".
#' @param alpha Optional numeric vector of non-negative weights of length
#'   `length(X)`. Required for `normalization = "custom"`; used as the IRLS
#'   target when `normalization = "balanced"`.
#' @param ridge Non-negative numeric scalar (or vector of length `length(X)`)
#'   controlling ridge stabilisation. The effective ridge used per block is
#'   scaled by the block's overall energy.
#' @param max_iter Maximum number of iterations for `normalization = "balanced"`.
#'   Ignored otherwise.
#' @param tol Relative tolerance on the gradient norm / objective change used
#'   to declare IRLS convergence for `normalization = "balanced"`.
#' @param verbose Logical; if `TRUE`, prints IRLS iteration progress.
#' @param use_future Logical; if `TRUE`, block-wise projector fitting is
#'   performed via `furrr::future_map()` when available.
#' @param block_weights Deprecated alias for `alpha`. When supplied, implicitly
#'   sets `normalization = "custom"` unless `normalization` is already explicit.
#' @param ... Unused (reserved for future extensions).
#'
#' @return An object inheriting from `multivarious::multiblock_biprojector` with
#'   additional class `"aligned_mcca"`. Relevant fields include `s`
#'   (`N × ncomp` reference-space scores), `v` (concatenated canonical
#'   directions), `sdev`, `block_indices`, `alpha_blocks` (block weights used),
#'   `block_weights` (alias), `normalization`, `canonical_weights`,
#'   `partial_scores`, `cor_loadings`, `block_contribs` (`ncomp × B` matrix of
#'   `s_k^T M_b s_k`), and (for `"balanced"`) `alpha_per_component`,
#'   `balance_trace`, `balance_converged`, `balance_iters`.
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
                         normalization = c("MFA", "balanced", "None", "custom"),
                         alpha = NULL,
                         ridge = 1e-6,
                         max_iter = 50,
                         tol = 1e-6,
                         verbose = FALSE,
                         use_future = FALSE,
                         block_weights = NULL,
                         ...) {
  normalization_missing <- missing(normalization)
  normalization <- match.arg(normalization)

  # ----- handle deprecated block_weights alias --------------------------------
  if (!is.null(block_weights)) {
    if (!is.null(alpha)) {
      stop("Supply either `alpha` or (deprecated) `block_weights`, not both.",
           call. = FALSE)
    }
    warning(
      "`block_weights` is deprecated in aligned_mcca(); use `alpha` with ",
      "`normalization = \"custom\"` instead.",
      call. = FALSE
    )
    alpha <- block_weights
    if (normalization_missing) normalization <- "custom"
  }

  # ----- validation / coercion ------------------------------------------------
  chk::chk_list(X)
  chk::chk_true(length(X) > 1)
  chk::chk_list(row_index)
  chk::chk_equal(length(X), length(row_index))
  chk::chk_flag(use_future)
  chk::chk_flag(verbose)
  max_iter <- as.integer(max_iter)
  chk::chk_integer(max_iter); chk::chk_gte(max_iter, 1)
  chk::chk_numeric(tol); chk::chk_gte(tol, 0)

  X <- lapply(X, function(x) {
    chk::chk_true(is.matrix(x) || is.data.frame(x))
    as.matrix(x)
  })

  S_blocks <- length(X)
  if (is.null(names(X))) names(X) <- paste0("B", seq_len(S_blocks))
  if (is.null(names(row_index))) {
    names(row_index) <- names(X)
  } else {
    row_index <- .lmfa_align_named_index_list(row_index, names(X), what = "row_index")
  }
  block_names <- names(X)

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

  # ----- ridge vector ---------------------------------------------------------
  if (length(ridge) == 1) ridge <- rep(ridge, S_blocks)
  ridge <- as.numeric(ridge)
  chk::chk_true(length(ridge) == S_blocks)
  chk::chk_true(all(is.finite(ridge))); chk::chk_true(all(ridge >= 0))

  # ----- alpha validation -----------------------------------------------------
  alpha_in <- NULL
  if (!is.null(alpha)) {
    alpha_in <- as.numeric(alpha)
    if (!is.null(names(alpha))) names(alpha_in) <- names(alpha)
    if (!is.null(names(alpha_in))) {
      if (!all(block_names %in% names(alpha_in))) {
        stop("Named `alpha` must include all X block names.", call. = FALSE)
      }
      alpha_in <- alpha_in[block_names]
    } else {
      if (length(alpha_in) != S_blocks) {
        stop("`alpha` must have length length(X).", call. = FALSE)
      }
      names(alpha_in) <- block_names
    }
    if (any(!is.finite(alpha_in)) || any(alpha_in < 0)) {
      stop("`alpha` must be finite and non-negative.", call. = FALSE)
    }
  }
  if (normalization == "custom" && is.null(alpha_in)) {
    stop("When normalization = 'custom', `alpha` must be provided.", call. = FALSE)
  }

  # ----- preprocess blocks ----------------------------------------------------
  if (is.list(preproc) && !inherits(preproc, "pre_processor") && !inherits(preproc, "prepper")) {
    chk::chk_equal(length(preproc), length(X))
  }
  prep_res <- prepare_block_preprocessors(X, preproc, check_consistent_ncol = FALSE)
  Xp <- prep_res$Xp
  proclist <- prep_res$proclist
  fitted_proclist <- .muscal_materialize_block_preprocessors(X, proclist)

  # Block indices for concatenated feature space
  block_indices <- list()
  current <- 1L
  for (k in seq_along(Xp)) {
    block_indices[[k]] <- current:(current + ncol(Xp[[k]]) - 1L)
    current <- current + ncol(Xp[[k]])
  }
  names(block_indices) <- block_names
  proc <- multivarious::concat_pre_processors(fitted_proclist, block_indices)

  # ----- fit per-block projectors once (unweighted; reused during IRLS) -------
  map_fun <- .mcca_map_fun(use_future)
  block_fits <- map_fun(seq_along(Xp), function(k) {
    .mcca_fit_one_block(
      X = Xp[[k]],
      ridge_mult = ridge[k],
      block_weight = 1,
      block_name = block_names[k],
      n = nrow(Xp[[k]])
    )
  })

  # Lifted N x N projector per block (block_weight baked in as 1 here)
  M_list <- lapply(seq_along(Xp), function(k) {
    .mcca_lift_to_N(block_fits[[k]]$P_raw, row_index[[k]], N)
  })
  names(M_list) <- block_names

  # ----- Resolve ncomp_eff ---------------------------------------------------
  chk::chk_numeric(ncomp)
  chk::chk_true(length(ncomp) == 1)
  chk::chk_true(is.finite(ncomp))
  chk::chk_gte(ncomp, 1)
  if (abs(ncomp - round(ncomp)) > 1e-8) stop("ncomp must be an integer-like value.", call. = FALSE)
  ncomp_eff <- min(as.integer(round(ncomp)), N)

  # ----- choose alpha_blocks via normalization -------------------------------
  balance_info <- list(trace = NULL, converged = NA, iters = 0L,
                       alpha_per_component = NULL)
  alpha_per_component <- NULL
  G <- NULL
  lambda <- NULL

  alpha_blocks <- switch(
    normalization,
    "None"     = rep(1, S_blocks),
    "MFA"      = unname(.mcca_mfa_weights(Xp)),
    "custom"   = unname(alpha_in),
    "balanced" = {
      balance_target <- if (!is.null(alpha_in)) unname(alpha_in) else rep(1, S_blocks)
      balance_w_init <- unname(.mcca_mfa_weights(Xp))
      balance_res <- .mcca_irls_balance(
        M_list = M_list,
        ncomp = ncomp_eff,
        target = balance_target,
        w_init = balance_w_init,
        max_iter = max_iter,
        tol = tol,
        verbose = verbose
      )
      G <- balance_res$G
      lambda <- balance_res$lambda
      alpha_per_component <- balance_res$weights_per_comp
      colnames(alpha_per_component) <- block_names
      rownames(alpha_per_component) <- paste0("Comp", seq_len(nrow(alpha_per_component)))
      balance_info <- list(
        trace = balance_res$trace,
        converged = balance_res$converged,
        converged_per_comp = balance_res$converged_per_comp,
        iters = balance_res$iters,
        iters_per_comp = balance_res$iters_per_comp,
        target = balance_target,
        w_init = balance_w_init,
        alpha_per_component = alpha_per_component
      )
      if (!isTRUE(balance_res$converged) && isTRUE(verbose)) {
        message("aligned_mcca: balanced IRLS did not converge for all components.")
      }
      balance_res$weights
    }
  )

  alpha_blocks <- as.numeric(alpha_blocks)
  names(alpha_blocks) <- block_names
  chk::chk_true(length(alpha_blocks) == S_blocks)
  if (any(!is.finite(alpha_blocks)) || any(alpha_blocks < 0)) {
    stop("Resolved block weights must be finite and non-negative.", call. = FALSE)
  }
  if (!any(alpha_blocks > 0)) {
    stop("At least one block weight must be strictly positive.", call. = FALSE)
  }

  # ----- assemble final H / eigenvectors -------------------------------------
  if (is.null(G)) {
    H <- matrix(0, nrow = N, ncol = N)
    for (b in seq_len(S_blocks)) H <- H + alpha_blocks[b] * M_list[[b]]
    H <- (H + t(H)) / 2

    eig <- .mcca_top_eigen(H, ncomp_eff, N)
    lambda <- eig$values
    G <- eig$vectors
  }
  lambda[!is.finite(lambda)] <- 0
  lambda[lambda < 0] <- 0
  sdev <- sqrt(lambda)

  S_scores <- G %*% diag(sdev, nrow = length(sdev), ncol = length(sdev))
  colnames(S_scores) <- paste0("Comp", seq_len(ncol(S_scores)))

  # ----- canonical weights / partial scores ----------------------------------
  p_tot <- sum(vapply(Xp, ncol, integer(1)))
  v_concat <- matrix(0, p_tot, ncomp_eff)

  partial_scores <- vector("list", S_blocks)
  names(partial_scores) <- block_names
  W_list <- vector("list", S_blocks)
  names(W_list) <- block_names

  inv_sdev <- ifelse(sdev > 1e-12, 1 / sdev, 0)
  for (k in seq_along(Xp)) {
    Xk <- Xp[[k]]
    Gk <- G[row_index[[k]], , drop = FALSE]
    Rchol <- block_fits[[k]]$chol
    A <- .mcca_chol_solve(Rchol, Gk)
    W_raw <- crossprod(Xk, A)
    W_list[[k]] <- W_raw

    alpha_col <- if (!is.null(alpha_per_component)) {
      alpha_per_component[, k]
    } else {
      rep(alpha_blocks[k], ncomp_eff)
    }
    W_weighted <- sweep(W_raw, 2, alpha_col, `*`)
    V_block <- sweep(W_weighted, 2, inv_sdev, `*`)
    v_concat[block_indices[[k]], ] <- V_block
    partial_scores[[k]] <- Xk %*% V_block
  }

  cor_list <- lapply(seq_along(Xp), function(k) {
    Ck <- tryCatch(
      suppressWarnings(stats::cor(Xp[[k]], S_scores[row_index[[k]], , drop = FALSE])),
      error = function(e) NULL
    )
    if (!is.null(Ck)) Ck[!is.finite(Ck)] <- 0
    Ck
  })
  cor_list <- cor_list[!vapply(cor_list, is.null, logical(1))]
  cor_loadings <- if (length(cor_list) > 0) do.call(rbind, cor_list) else NULL

  # Per-block per-component contribution matrix (ncomp x B)
  block_contribs <- vapply(M_list, function(Mb) {
    as.numeric(colSums(G * (Mb %*% G)))
  }, numeric(ncomp_eff))
  if (is.null(dim(block_contribs))) {
    block_contribs <- matrix(block_contribs, nrow = ncomp_eff, ncol = S_blocks)
  }
  rownames(block_contribs) <- paste0("Comp", seq_len(ncomp_eff))
  colnames(block_contribs) <- block_names

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
    block_weights = unname(alpha_blocks),
    alpha_blocks = alpha_blocks,
    normalization = normalization,
    row_index = row_index,
    N = N,
    partial_scores = partial_scores,
    canonical_weights = W_list,
    cor_loadings = cor_loadings,
    block_contribs = block_contribs,
    alpha_per_component = alpha_per_component,
    balance_trace = balance_info$trace,
    balance_converged = balance_info$converged,
    balance_converged_per_comp = balance_info$converged_per_comp,
    balance_iters = balance_info$iters,
    balance_iters_per_comp = balance_info$iters_per_comp,
    balance_target = if (!is.null(balance_info$target)) balance_info$target else NULL,
    names = block_names,
    classes = "aligned_mcca"
  )
}

#' Anchored Multiblock Canonical Correlation Analysis (Anchored MCCA)
#'
#' @description
#' Anchored MCCA estimates a shared canonical score matrix `S` in the row space
#' of a privileged anchor block `Y` (`N × q`), while requiring one or more linked
#' blocks `X_k` (with row-to-anchor maps `row_index[[k]]`) to agree with `S`.
#'
#' Compared with [aligned_mcca()] (which treats all blocks symmetrically) and
#' [anchored_mfa()] (reconstruction-driven, optimises summed explained variance),
#' Anchored MCCA is correlation-driven and anchor-privileged: it looks for
#' components that every block can *predict*, not just explain in isolation.
#'
#' @details
#' ## Interface parity with [anchored_mfa()]
#' The signature mirrors [anchored_mfa()] (`Y, X, row_index, preproc, ncomp,`
#' `normalization, alpha, ridge, max_iter, tol, verbose, use_future`) so that
#' the two anchored methods can be swapped with minimal code changes. Arguments
#' that are specific to the reconstruction-ALS machinery of anchored_mfa
#' (`score_constraint`, `feature_groups`, `feature_lambda`) are intentionally
#' absent here; MCCA uses closed-form eigendecomposition rather than ALS.
#'
#' ## Block-weighting schemes
#' The block-weight vector `alpha_blocks` (one weight per block, first entry = `Y`)
#' is formed according to `normalization`:
#' * `"MFA"` (default) — `alpha_b = 1 / (sigma_1(X_b)^2)`, where `sigma_1` is the
#'   first singular value of the preprocessed block. This is the CCA analogue of
#'   the MFA normalization in [anchored_mfa()] and prevents a block from
#'   dominating purely by scale.
#' * `"balanced"` — optimizes each component with a geometric-mean/log
#'   contribution objective, `sum_b beta_b log(g_j^T M_b g_j)`, under
#'   orthogonality to earlier components. This discourages block-specific
#'   directions because every positive-target block must contribute to the same
#'   component.
#' * `"None"` — uniform weights `alpha_b = 1`.
#' * `"custom"` — use `alpha` as supplied.
#'
#' ## Model
#' Let \eqn{P_b = X_b (X_b^\top X_b + \kappa_b I)^-1 X_b^\top} be the
#' ridge-stabilised hat matrix of block `b`, lifted to the anchor row space via
#' `row_index[[b]]` to produce `M_b` (`N × N`). Anchored MCCA computes the top
#' eigenpairs of
#' \deqn{H = \sum_b alpha_b M_b}
#' with `alpha_b` chosen by `normalization`. Y enters as the first block with
#' identity row-index, so its contribution is exact and the shared scores live
#' in its row space.
#'
#' For `normalization = "balanced"`, the method instead solves the per-component
#' block-balanced objective above and reports both raw contributions
#' `block_contribs` and criterion-weighted contributions
#' `weighted_block_contribs`.
#'
#' ## What balanced mode constrains
#' Balanced mode does **not** constrain the sum of squares or Euclidean norms of
#' the canonical weight vectors. The balance criterion is applied in the anchor
#' row space through each block's lifted ridge projector. For a unit score
#' direction `g_j`, block `b` contributes
#' \deqn{c_{bj} = g_j^\top M_b g_j.}
#' The balanced objective penalizes directions where any positive-target block
#' has near-zero \eqn{c_{bj}}, so the selected component must be supportable by all
#' participating blocks. Canonical weights are then computed afterward as the
#' ridge back-projection needed for each block to express that already-chosen
#' score direction. Their norms remain scale- and ridge-dependent diagnostics,
#' not constraints.
#'
#' ## High-dimensional stability (p > n)
#' The implementation works in observation space and uses ridge-stabilised
#' solves on `n_b × n_b` systems, so `p_b >> n_b` blocks are handled without
#' special casing. The kernel matrix `K_b = X_b X_b^T` has rank at most `n_b`;
#' ridge stabilisation (`kappa_b = ridge_b * mean(diag(K_b))`, with an
#' automatic ridge floor if `ridge = 0` yields a singular system) keeps the
#' Cholesky well-posed. All four `normalization` modes inherit this
#' stability — `"MFA"` only needs the first singular value per block (via
#' `multivarious::svd_wrapper()` which picks `base` or `irlba` automatically),
#' and `"balanced"` reuses the same cached projectors across IRLS iterations.
#'
#' @param Y Numeric matrix/data.frame (`N × q`) serving as the anchor block.
#' @param X A list of numeric matrices/data.frames. Each element `X[[k]]` is
#'   `n_k × p_k`.
#' @param row_index A list of integer vectors. `row_index[[k]]` has length `n_k`
#'   and maps rows of `X[[k]]` to rows of `Y` (values in `1..N`).
#' @param preproc A `multivarious` preprocessing pipeline (a `pre_processor`/`prepper`)
#'   or a list of them. If a list, it must have length `1 + length(X)` and will be
#'   applied to `c(list(Y), X)` in that order.
#' @param ncomp Integer; number of canonical components to compute.
#' @param normalization Block weighting scheme. One of `"MFA"` (default, scale-
#'   normalised via `1/sigma_1^2`), `"balanced"` (per-component geometric
#'   contribution balancing to prevent single-block domination), `"None"`
#'   (uniform weights), or `"custom"` (use `alpha`).
#' @param alpha Optional numeric vector of non-negative weights. If unnamed and of
#'   length `length(X)` it is interpreted as weights for the X blocks with `Y`
#'   receiving weight `1`; if of length `1 + length(X)` the first entry applies
#'   to `Y`. Required when `normalization = "custom"`; used as the IRLS prior
#'   when `normalization = "balanced"`.
#' @param ridge Non-negative numeric scalar (or vector of length `1 + length(X)`)
#'   controlling ridge stabilisation. The effective ridge used per block is
#'   scaled by the block's overall energy so the default works across variable
#'   scalings.
#' @param max_iter Maximum number of IRLS iterations for `normalization = "balanced"`.
#'   Ignored for other modes.
#' @param tol Relative tolerance on the projected-gradient norm and objective
#'   change used to declare convergence for `normalization = "balanced"`.
#' @param verbose Logical; if `TRUE`, prints IRLS iteration progress.
#' @param use_future Logical; if `TRUE`, block-wise projector fitting is performed
#'   via `furrr::future_map()` when available.
#' @param block_weights Deprecated alias for `alpha`. When supplied, implicitly
#'   sets `normalization = "custom"` unless `normalization` is already explicit.
#' @param ... Unused (reserved for future extensions).
#'
#' @return An object of class `"anchored_mcca"` (and `"aligned_mcca"`,
#'   `"multiblock_biprojector"`). Relevant fields include `s` (anchor-space
#'   scores, `N × ncomp`), `v` (concatenated canonical directions), `sdev`,
#'   `block_indices`, `alpha_blocks`, `block_weights` (alias), `normalization`,
#'   `canonical_weights`, `partial_scores`, `cor_loadings`, `block_contribs`
#'   (`ncomp × B` matrix of raw `g_j^T M_b g_j` contributions),
#'   `weighted_block_contribs`, `block_contrib_fraction`, and (for
#'   `"balanced"` only) `balance_trace`, `balance_converged`, `balance_iters`.
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
#' fit <- anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
#'                      ncomp = 2, normalization = "balanced")
#' stopifnot(nrow(multivarious::scores(fit)) == N)
#' }
#' @export
anchored_mcca <- function(Y,
                          X,
                          row_index,
                          preproc = multivarious::center(),
                          ncomp = 2,
                          normalization = c("MFA", "balanced", "None", "custom"),
                          alpha = NULL,
                          ridge = 1e-6,
                          max_iter = 50,
                          tol = 1e-6,
                          verbose = FALSE,
                          use_future = FALSE,
                          block_weights = NULL,
                          ...) {
  fit_call <- match.call(expand.dots = FALSE)
  fit_dots <- list(...)
  normalization_missing <- missing(normalization)
  normalization <- match.arg(normalization)

  # ----- handle deprecated block_weights alias --------------------------------
  if (!is.null(block_weights)) {
    if (!is.null(alpha)) {
      stop("Supply either `alpha` or (deprecated) `block_weights`, not both.",
           call. = FALSE)
    }
    warning(
      "`block_weights` is deprecated in anchored_mcca(); use `alpha` with ",
      "`normalization = \"custom\"` instead.",
      call. = FALSE
    )
    alpha <- block_weights
    if (normalization_missing) normalization <- "custom"
  }

  # ----- validation / coercion ------------------------------------------------
  chk::chk_true(is.matrix(Y) || is.data.frame(Y))
  Y <- as.matrix(Y)
  chk::chk_list(X)
  chk::chk_true(length(X) >= 1)
  chk::chk_list(row_index)
  chk::chk_equal(length(X), length(row_index))
  ncomp <- as.integer(ncomp)
  chk::chk_integer(ncomp); chk::chk_gte(ncomp, 1)
  max_iter <- as.integer(max_iter)
  chk::chk_integer(max_iter); chk::chk_gte(max_iter, 1)
  chk::chk_numeric(tol); chk::chk_gte(tol, 0)
  chk::chk_flag(verbose)
  chk::chk_flag(use_future)

  X <- lapply(X, function(x) {
    chk::chk_true(is.matrix(x) || is.data.frame(x))
    as.matrix(x)
  })
  if (is.null(names(X))) names(X) <- paste0("X", seq_along(X))
  if (is.null(names(row_index))) {
    names(row_index) <- names(X)
  } else {
    row_index <- .lmfa_align_named_index_list(row_index, names(X), what = "row_index")
  }

  N <- nrow(Y)
  for (k in seq_along(X)) {
    idx <- as.integer(row_index[[k]])
    chk::chk_equal(length(idx), nrow(X[[k]]))
    if (anyNA(idx)) stop("row_index cannot contain NA.", call. = FALSE)
    if (any(idx < 1L) || any(idx > N)) stop("row_index values must be in 1..nrow(Y).", call. = FALSE)
    row_index[[k]] <- idx
  }

  data_refit <- list(
    Y = Y,
    X = X,
    row_index = row_index
  )

  X_all <- c(list(Y = Y), X)
  idx_all <- c(list(Y = seq_len(N)), row_index)
  B_all <- length(X_all)
  block_names <- names(X_all)

  # ----- resolve ridge vector -------------------------------------------------
  if (length(ridge) == 1) ridge <- rep(ridge, B_all)
  ridge <- as.numeric(ridge)
  chk::chk_true(length(ridge) == B_all)
  chk::chk_true(all(is.finite(ridge))); chk::chk_true(all(ridge >= 0))

  # ----- expand alpha into a vector aligned with X_all ------------------------
  alpha_in <- NULL
  if (!is.null(alpha)) {
    alpha_in <- as.numeric(alpha)
    if (!is.null(names(alpha))) names(alpha_in) <- names(alpha)

    if (is.null(names(alpha_in))) {
      if (length(alpha_in) == length(X)) {
        alpha_in <- c(1, alpha_in)
      } else if (length(alpha_in) != B_all) {
        stop("`alpha` must have length length(X) or 1 + length(X).", call. = FALSE)
      }
      names(alpha_in) <- block_names
    } else {
      nms <- names(alpha_in)
      if (!("Y" %in% nms)) {
        if (!all(names(X) %in% nms)) {
          stop("Named `alpha` must include all X block names.", call. = FALSE)
        }
        alpha_in <- c(Y = 1, alpha_in[names(X)])
      } else {
        alpha_in <- alpha_in[block_names]
      }
    }
    if (any(!is.finite(alpha_in)) || any(alpha_in < 0)) {
      stop("`alpha` must be finite and non-negative.", call. = FALSE)
    }
  }

  if (normalization == "custom" && is.null(alpha_in)) {
    stop("When normalization = 'custom', `alpha` must be provided.", call. = FALSE)
  }

  # ----- preprocess blocks ----------------------------------------------------
  if (is.list(preproc) && !inherits(preproc, "pre_processor") && !inherits(preproc, "prepper")) {
    chk::chk_equal(length(preproc), B_all)
  }
  prep_res <- prepare_block_preprocessors(X_all, preproc, check_consistent_ncol = FALSE)
  Xp <- prep_res$Xp
  proclist <- prep_res$proclist
  fitted_proclist <- .muscal_materialize_block_preprocessors(X_all, proclist)

  # Concatenated feature indices
  block_indices <- list()
  current <- 1L
  for (k in seq_along(Xp)) {
    block_indices[[k]] <- current:(current + ncol(Xp[[k]]) - 1L)
    current <- current + ncol(Xp[[k]])
  }
  names(block_indices) <- block_names
  proc <- multivarious::concat_pre_processors(fitted_proclist, block_indices)

  # ----- fit per-block projectors once (unweighted, reused during IRLS) -------
  map_fun <- .mcca_map_fun(use_future)
  block_fits <- map_fun(seq_along(Xp), function(k) {
    .mcca_fit_one_block(
      X = Xp[[k]],
      ridge_mult = ridge[k],
      block_weight = 1,
      block_name = block_names[k],
      n = nrow(Xp[[k]])
    )
  })

  # Lifted N x N projector per block (block_weight baked in as 1 here)
  M_list <- lapply(seq_along(Xp), function(k) {
    .mcca_lift_to_N(block_fits[[k]]$P_raw, idx_all[[k]], N)
  })
  names(M_list) <- block_names

  # ----- choose alpha_blocks according to normalization -----------------------
  balance_info <- list(trace = NULL, converged = NA, iters = 0L,
                       alpha_per_component = NULL)
  alpha_per_component <- NULL  # (ncomp × B) when balanced; NULL otherwise

  ncomp_eff <- min(ncomp, N)

  G <- NULL
  lambda <- NULL

  alpha_blocks <- switch(
    normalization,
    "None"     = rep(1, B_all),
    "MFA"      = unname(.mcca_mfa_weights(Xp)),
    "custom"   = unname(alpha_in),
    "balanced" = {
      # Target contribution ratio β_b: user-supplied alpha (if any) or uniform.
      # Starting weights: MFA-normalised for fast convergence on mixed-scale blocks.
      balance_target <- if (!is.null(alpha_in)) unname(alpha_in) else rep(1, B_all)
      balance_w_init <- unname(.mcca_mfa_weights(Xp))
      balance_res <- .mcca_irls_balance(
        M_list = M_list,
        ncomp = ncomp_eff,
        target = balance_target,
        w_init = balance_w_init,
        max_iter = max_iter,
        tol = tol,
        verbose = verbose
      )
      G <- balance_res$G
      lambda <- balance_res$lambda
      alpha_per_component <- balance_res$weights_per_comp
      colnames(alpha_per_component) <- block_names
      rownames(alpha_per_component) <- paste0("Comp", seq_len(nrow(alpha_per_component)))
      balance_info <- list(
        trace = balance_res$trace,
        converged = balance_res$converged,
        converged_per_comp = balance_res$converged_per_comp,
        iters = balance_res$iters,
        iters_per_comp = balance_res$iters_per_comp,
        target = balance_target,
        w_init = balance_w_init,
        alpha_per_component = alpha_per_component
      )
      if (!isTRUE(balance_res$converged) && isTRUE(verbose)) {
        message("anchored_mcca: balanced IRLS did not converge for all components.")
      }
      balance_res$weights
    }
  )

  alpha_blocks <- as.numeric(alpha_blocks)
  names(alpha_blocks) <- block_names
  chk::chk_true(length(alpha_blocks) == B_all)
  if (any(!is.finite(alpha_blocks)) || any(alpha_blocks < 0)) {
    stop("Resolved block weights must be finite and non-negative.", call. = FALSE)
  }
  if (!any(alpha_blocks > 0)) {
    stop("At least one block weight must be strictly positive.", call. = FALSE)
  }

  # ----- assemble final H / eigenvectors ---------------------------------------
  # For non-balanced modes we solve a single eigenproblem on H = Σ w_b M_b.
  # For "balanced" mode we reuse G / λ from the per-component IRLS.
  if (is.null(G)) {
    H <- matrix(0, N, N)
    for (b in seq_len(B_all)) H <- H + alpha_blocks[b] * M_list[[b]]
    H <- (H + t(H)) / 2

    eig <- .mcca_top_eigen(H, ncomp_eff, N)
    lambda <- eig$values
    G <- eig$vectors
  }
  lambda[!is.finite(lambda)] <- 0
  lambda[lambda < 0] <- 0
  sdev <- sqrt(lambda)

  S_scores <- G %*% diag(sdev, nrow = length(sdev), ncol = length(sdev))
  colnames(S_scores) <- paste0("Comp", seq_len(ncol(S_scores)))

  # ----- loadings / canonical weights / partial scores ------------------------
  # alpha_col[k, b] is the per-component weight used to scale block b's
  # contribution to component k. In "balanced" mode this varies by component
  # (taken from alpha_per_component); otherwise it is constant across components.
  p_tot <- sum(vapply(Xp, ncol, integer(1)))
  v_concat <- matrix(0, p_tot, ncomp_eff)
  partial_scores <- vector("list", B_all)
  W_list <- vector("list", B_all)
  names(partial_scores) <- block_names
  names(W_list) <- block_names

  inv_sdev <- ifelse(sdev > 1e-12, 1 / sdev, 0)
  for (k in seq_len(B_all)) {
    Xk <- Xp[[k]]
    Gk <- G[idx_all[[k]], , drop = FALSE]
    Rchol <- block_fits[[k]]$chol
    A <- .mcca_chol_solve(Rchol, Gk)
    W_raw <- crossprod(Xk, A)
    W_list[[k]] <- W_raw

    alpha_col <- if (!is.null(alpha_per_component)) {
      alpha_per_component[, k]
    } else {
      rep(alpha_blocks[k], ncomp_eff)
    }
    W_weighted <- sweep(W_raw, 2, alpha_col, `*`)
    V_block <- sweep(W_weighted, 2, inv_sdev, `*`)
    v_concat[block_indices[[k]], ] <- V_block
    partial_scores[[k]] <- Xk %*% V_block
  }

  # cor_loadings (feature-vs-score correlations)
  cor_list <- lapply(seq_len(B_all), function(k) {
    Ck <- tryCatch(
      suppressWarnings(stats::cor(Xp[[k]], S_scores[idx_all[[k]], , drop = FALSE])),
      error = function(e) NULL
    )
    if (!is.null(Ck)) Ck[!is.finite(Ck)] <- 0
    Ck
  })
  cor_list <- cor_list[!vapply(cor_list, is.null, logical(1))]
  cor_loadings <- if (length(cor_list) > 0) do.call(rbind, cor_list) else NULL

  # Per-block per-component contribution matrix (ncomp x B). Raw contributions
  # are g_j^T M_b g_j; weighted contributions include the block/component weight
  # used by the fitted criterion.
  block_contribs <- vapply(M_list, function(Mb) {
    as.numeric(colSums(G * (Mb %*% G)))
  }, numeric(ncomp_eff))
  if (is.null(dim(block_contribs))) {
    block_contribs <- matrix(block_contribs, nrow = ncomp_eff, ncol = B_all)
  }
  rownames(block_contribs) <- paste0("Comp", seq_len(ncomp_eff))
  colnames(block_contribs) <- block_names

  alpha_mat <- if (!is.null(alpha_per_component)) {
    alpha_per_component
  } else {
    matrix(
      rep(alpha_blocks, each = ncomp_eff),
      nrow = ncomp_eff,
      ncol = B_all,
      dimnames = dimnames(block_contribs)
    )
  }
  weighted_block_contribs <- block_contribs * alpha_mat
  block_contrib_fraction <- weighted_block_contribs
  row_totals <- rowSums(block_contrib_fraction)
  nz <- is.finite(row_totals) & row_totals > 0
  block_contrib_fraction[nz, ] <- sweep(block_contrib_fraction[nz, , drop = FALSE], 1, row_totals[nz], `/`)
  block_contrib_fraction[!nz, ] <- 0

  out <- multivarious::multiblock_biprojector(
    v = v_concat,
    s = S_scores,
    sdev = sdev,
    preproc = proc,
    block_indices = block_indices,
    # anchored_mcca / aligned_mcca metadata
    lambda = lambda,
    ridge = ridge,
    kappa = vapply(block_fits, `[[`, numeric(1), "kappa"),
    block_weights = unname(alpha_blocks),
    alpha_blocks = alpha_blocks,
    normalization = normalization,
    row_index = idx_all,
    N = N,
    partial_scores = partial_scores,
    canonical_weights = W_list,
    cor_loadings = cor_loadings,
    block_contribs = block_contribs,
    weighted_block_contribs = weighted_block_contribs,
    block_contrib_fraction = block_contrib_fraction,
    alpha_per_component = alpha_per_component,
    balance_trace = balance_info$trace,
    balance_converged = balance_info$converged,
    balance_converged_per_comp = balance_info$converged_per_comp,
    balance_iters = balance_info$iters,
    balance_iters_per_comp = balance_info$iters_per_comp,
    balance_target = if (!is.null(balance_info$target)) balance_info$target else NULL,
    names = block_names,
    classes = c("anchored_mcca", "aligned_mcca")
  )

  refit_fn <- .mcca_make_anchored_mcca_refit_fn(
    preproc = preproc,
    ncomp = ncomp,
    normalization = normalization,
    alpha = alpha_in,
    ridge = ridge,
    max_iter = max_iter,
    tol = tol,
    use_future = use_future,
    fit_dots = fit_dots
  )

  .muscal_attach_fit_contract(
    out,
    method = "anchored_mcca",
    task = "association",
    oos_types = "scores",
    fit_call = fit_call,
    refit_supported = TRUE,
    prediction_target = "shared_scores",
    refit = .muscal_make_refit_spec(
      data = data_refit,
      fit_fn = refit_fn,
      bootstrap_fn = .muscal_bootstrap_anchored_data,
      permutation_fn = .muscal_permute_anchored_data,
      resample_unit = "anchor_rows"
    )
  )
}
