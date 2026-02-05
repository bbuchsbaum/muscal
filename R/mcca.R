#' @md
#' @rdname mcca
#' @details
#' The `mcca.list` method applies MCCA to a list of data matrices or data frames.
#' This method first converts the list to a `multiblock` object and then calls
#' `mcca.multiblock()`.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- replicate(3, matrix(rnorm(50 * 20), 50, 20), simplify = FALSE)
#' fit <- mcca(X, ncomp = 3)
#' }
#' @export
mcca.list <- function(data, preproc = multivarious::center(), ncomp = 2,
                      ridge = 1e-6, block_weights = NULL, use_future = FALSE, ...) {
  data <- multidesign::multiblock(data)
  mcca.multiblock(
    data = data,
    preproc = preproc,
    ncomp = ncomp,
    ridge = ridge,
    block_weights = block_weights,
    use_future = use_future,
    ...
  )
}

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.mcca_as_dense <- function(x) {
  if (inherits(x, "Matrix")) return(as.matrix(x))
  x
}

.mcca_chol_solve <- function(R, B) {
  # Solve (t(R) %*% R) X = B where R is upper-triangular chol factor.
  y <- forwardsolve(t(R), B, upper.tri = FALSE, transpose = FALSE)
  backsolve(R, y, upper.tri = TRUE, transpose = FALSE)
}

.mcca_fit_one_block <- function(X, ridge_mult, block_weight, block_name, n) {
  K <- .mcca_as_dense(Matrix::tcrossprod(X))
  if (!all(is.finite(K))) stop("Non-finite values encountered in block Gram matrix.")

  scale_k <- mean(diag(K))
  if (!is.finite(scale_k) || scale_k <= 0) scale_k <- 1

  kappa <- ridge_mult * scale_k
  if (!is.finite(kappa) || kappa < 0) {
    stop("ridge must be non-negative and finite.")
  }

  # If ridge is zero, try the unregularized system first (may be ill-posed after
  # centering); fall back to a tiny ridge for numerical invertibility.
  K_reg <- NULL
  R <- NULL
  if (kappa == 0) {
    R0 <- tryCatch(chol(K), error = function(e) NULL)
    if (!is.null(R0)) {
      K_reg <- K
      R <- R0
    } else {
      kappa <- 1e-8 * scale_k
      warning(
        sprintf(
          "mcca: ridge=0 yields a singular system for block '%s'; using ridge=1e-8 (scaled).",
          block_name
        ),
        call. = FALSE
      )
      K_reg <- K + kappa * diag(n)
    }
  } else {
    K_reg <- K + kappa * diag(n)
  }

  if (is.null(R)) {
    R <- tryCatch(
      chol(K_reg),
      error = function(e) {
        # Retry with increased ridge if needed
        kappa2 <- max(kappa, 1e-6 * scale_k) * 100
        K_reg2 <- K + kappa2 * diag(n)
        R2 <- chol(K_reg2)
        warning(
          sprintf(
            "mcca: increased ridge for block '%s' from %.3g to %.3g for numerical stability.",
            block_name, kappa, kappa2
          ),
          call. = FALSE
        )
        kappa <<- kappa2
        R2
      }
    )
  }

  # P = (K + kappa I)^{-1} K; K commutes with K + kappa I, so this equals K (K + kappa I)^{-1}.
  P <- .mcca_chol_solve(R, K)
  P <- (P + t(P)) / 2

  list(
    P = block_weight * P,
    chol = R,
    kappa = kappa,
    ridge_mult = ridge_mult,
    scale_k = scale_k,
    block_weight = block_weight
  )
}

#' @md
#' @rdname mcca
#' @details
#' The `mcca.multiblock` method implements a ridge-stabilized MAXVAR generalized
#' CCA. Internally it computes a compromise score space by eigendecomposition of
#' a weighted sum of block projection operators:
#' \deqn{H = \sum_s w_s \, X_s (X_s^T X_s + \kappa_s I)^{-1} X_s^T.}
#'
#' The ridge parameter is applied per block as \eqn{\kappa_s = \text{ridge}_s \cdot
#' \bar{d}_s}, where \eqn{\bar{d}_s} is the mean diagonal of the block Gram matrix
#' \eqn{X_s X_s^T}. This scaling makes the default ridge behave reasonably across
#' different variable scalings.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' blocks <- replicate(3, matrix(rnorm(20 * 50), 20, 50), simplify = FALSE)
#' fit <- mcca(blocks, ncomp = 2)
#' stopifnot(ncol(multivarious::scores(fit)) == 2)
#' }
#' @export
mcca.multiblock <- function(data, preproc = multivarious::center(), ncomp = 2,
                            ridge = 1e-6, block_weights = NULL,
                            use_future = FALSE, ...) {
  chk::chk_true(length(data) > 1)
  for (i in seq_along(data)) {
    chk::chkor_vld(chk::vld_matrix(data[[i]]), chk::vld_s4_class(data[[i]], "Matrix"))
  }

  nrs <- sapply(data, nrow)
  chk::chk_true(all(nrs == nrs[1]))
  n <- nrs[1]

  S <- length(data)
  if (is.null(names(data))) names(data) <- paste0("B", seq_len(S))

  if (is.null(block_weights)) block_weights <- rep(1, S)
  chk::chk_numeric(block_weights)
  chk::chk_true(length(block_weights) == S)
  chk::chk_true(all(is.finite(block_weights)))
  chk::chk_true(all(block_weights >= 0))

  # Preprocess blocks independently
  prep <- prepare_block_preprocessors(data, preproc, check_consistent_ncol = FALSE)
  strata <- prep$Xp
  proclist <- prep$proclist

  # Block indices for concatenated loadings space
  block_indices <- list()
  ind <- 1
  for (i in seq_along(strata)) {
    block_indices[[i]] <- seq(ind, ind + ncol(strata[[i]]) - 1)
    ind <- ind + ncol(strata[[i]])
  }
  names(block_indices) <- names(data)

  proc <- multivarious::concat_pre_processors(proclist, block_indices)

  # Ridge vector (scaled per block)
  if (length(ridge) == 1) ridge <- rep(ridge, S)
  chk::chk_numeric(ridge)
  chk::chk_true(length(ridge) == S)
  chk::chk_true(all(is.finite(ridge)))
  chk::chk_true(all(ridge >= 0))

  idx <- seq_along(strata)
  map_fun <- if (isTRUE(use_future)) {
    function(.x, .f) {
      furrr::future_map(.x, .f, .options = furrr::furrr_options(seed = TRUE))
    }
  } else {
    function(.x, .f) lapply(.x, .f)
  }

  block_fits <- map_fun(idx, function(i) {
    .mcca_fit_one_block(
      X = strata[[i]],
      ridge_mult = ridge[i],
      block_weight = block_weights[i],
      block_name = names(data)[i],
      n = n
    )
  })

  H <- Reduce(`+`, lapply(block_fits, `[[`, "P"))
  H <- (H + t(H)) / 2

  # Cap ncomp
  chk::chk_numeric(ncomp)
  chk::chk_true(length(ncomp) == 1)
  chk::chk_true(is.finite(ncomp))
  chk::chk_gte(ncomp, 1)
  if (abs(ncomp - round(ncomp)) > 1e-8) {
    stop("ncomp must be an integer-like value.", call. = FALSE)
  }
  ncomp_eff <- min(as.integer(round(ncomp)), n)

  eig <- NULL
  if (ncomp_eff < n) {
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

  # Scale scores like PCA: S = U * sdev so that crossprod(S) = diag(lambda).
  S_scores <- G %*% diag(sdev, nrow = length(sdev), ncol = length(sdev))
  colnames(S_scores) <- paste0("Comp", seq_len(ncol(S_scores)))

  # Build concatenated loading matrix v (so that X_concat %*% v == scores)
  p_tot <- sum(vapply(strata, ncol, integer(1)))
  v_concat <- matrix(0, p_tot, ncomp_eff)

  partial_scores <- vector("list", S)
  names(partial_scores) <- names(data)

  W_list <- vector("list", S)
  names(W_list) <- names(data)

  for (i in seq_along(strata)) {
    X <- strata[[i]]
    R <- block_fits[[i]]$chol
    A <- .mcca_chol_solve(R, G) # (K + kappa I)^{-1} G
    W_raw <- crossprod(X, A)    # p_i x k (canonical weights in feature space)
    W_list[[i]] <- W_raw

    # Scale so that concatenated scores match S_scores.
    inv_sdev <- ifelse(sdev > 1e-12, 1 / sdev, 0)
    W_weighted <- block_fits[[i]]$block_weight * W_raw
    V_block <- sweep(W_weighted, 2, inv_sdev, `*`)
    v_concat[block_indices[[i]], ] <- V_block
    partial_scores[[i]] <- X %*% V_block
  }

  cor_loadings <- tryCatch(
    {
      X_concat <- do.call(cbind, strata)
      C <- suppressWarnings(stats::cor(X_concat, S_scores))
      C[!is.finite(C)] <- 0
      C
    },
    error = function(e) NULL
  )

  rv <- tryCatch(
    {
      M <- compute_sim_mat(strata, function(x1, x2) MatrixCorrelation::RV(x1, x2))
      diag(M) <- 1
      rownames(M) <- colnames(M) <- names(data)
      M
    },
    error = function(e) NULL
  )

  rv2 <- tryCatch(
    {
      M <- compute_sim_mat(strata, function(x1, x2) MatrixCorrelation::RV2(x1, x2))
      diag(M) <- 1
      rownames(M) <- colnames(M) <- names(data)
      M
    },
    error = function(e) NULL
  )

  multivarious::multiblock_biprojector(
    v = v_concat,
    s = S_scores,
    sdev = sdev,
    preproc = proc,
    block_indices = block_indices,
    # MCCA metadata / diagnostics
    lambda = lambda,
    ridge = ridge,
    kappa = vapply(block_fits, `[[`, numeric(1), "kappa"),
    block_weights = block_weights,
    partial_scores = partial_scores,
    canonical_weights = W_list,
    cor_loadings = cor_loadings,
    rv = rv,
    rv2 = rv2,
    names = names(data),
    classes = "mcca"
  )
}
