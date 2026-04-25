library(testthat)
library(muscal)

.expand_rowsum_to_N <- function(x, idx, N) {
  rs <- rowsum(x, idx, reorder = FALSE)
  out <- matrix(0, nrow = N, ncol = ncol(x))
  out[as.integer(rownames(rs)), ] <- rs
  out
}

.permute_reference_labels <- function(idx, perm) {
  perm[idx]
}

.sim_aligned_mcca <- function(N, k, p_vec, n_vec, noise_sd = 0.05) {
  Z <- matrix(rnorm(N * k), nrow = N, ncol = k)
  idx_list <- lapply(n_vec, function(nk) sample.int(N, nk, replace = TRUE))
  blocks <- Map(function(p, idx) {
    A <- matrix(rnorm(p * k), nrow = p, ncol = k)
    signal <- Z[idx, , drop = FALSE] %*% t(A)
    signal + matrix(rnorm(length(idx) * p, sd = noise_sd), nrow = length(idx), ncol = p)
  }, p_vec, idx_list)
  list(blocks = blocks, idx = idx_list, Z = Z)
}

test_that("aligned_mcca validates row_index and returns expected structure", {
  set.seed(1)
  N <- 30
  X1 <- matrix(rnorm(10 * 6), 10, 6)
  X2 <- matrix(rnorm(12 * 5), 12, 5)
  idx1 <- sample.int(N, nrow(X1), replace = FALSE)
  idx2 <- sample.int(N, nrow(X2), replace = FALSE)

  expect_error(aligned_mcca(list(X1, X2), list(idx1)))
  expect_error(aligned_mcca(list(X1, X2), list(idx1, c(idx2, 1L))))
  expect_error(aligned_mcca(list(X1, X2), list(idx1, rep(999L, nrow(X2))), N = N))
  expect_error(aligned_mcca(list(X1, X2), list(idx1, c(idx2[1], NA_integer_, idx2[-c(1, 2)]))))
  expect_error(aligned_mcca(list(X1, X2), list(idx1, idx2), N = 5))

  fit <- aligned_mcca(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), N = N, ncomp = 2)
  expect_s3_class(fit, "aligned_mcca")
  expect_true(inherits(fit, "multiblock_biprojector"))
  expect_equal(length(fit$block_indices), 2)
  expect_equal(names(fit$block_indices), c("X1", "X2"))
  expect_equal(nrow(multivarious::scores(fit)), N)
  expect_equal(multivarious::ncomp(fit), 2)
  expect_equal(fit$N, N)
  expect_equal(length(fit$partial_scores), 2)
  expect_equal(length(fit$canonical_weights), 2)
})

test_that("aligned_mcca infers N when not supplied", {
  set.seed(1)
  X1 <- matrix(rnorm(10 * 6), 10, 6)
  X2 <- matrix(rnorm(12 * 5), 12, 5)
  idx1 <- c(1:9, 15)
  idx2 <- sample.int(15, nrow(X2), replace = TRUE)

  fit <- aligned_mcca(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
  expect_equal(fit$N, 15)
  expect_equal(nrow(multivarious::scores(fit)), 15)

  idx_bad <- idx2
  idx_bad[1] <- 0L
  expect_error(aligned_mcca(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx_bad)))
})

test_that("aligned_mcca supports per-block preproc lists and validates ncomp", {
  set.seed(2)
  N <- 20
  X1 <- matrix(rnorm(N * 6), N, 6)
  X2 <- matrix(rnorm(N * 5), N, 5)
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  fit <- aligned_mcca(
    list(X1 = X1, X2 = X2),
    idx,
    N = N,
    ncomp = 2,
    preproc = list(multivarious::center(), multivarious::center())
  )
  expect_s3_class(fit, "aligned_mcca")

  expect_error(aligned_mcca(list(X1 = X1, X2 = X2), idx, N = N, ncomp = 1.5))
})

test_that("aligned_mcca with normalization='None' agrees with mcca when all blocks share rows", {
  set.seed(2)
  N <- 40
  X1 <- matrix(rnorm(N * 8), N, 8)
  X2 <- matrix(rnorm(N * 5), N, 5)
  blocks <- list(X1 = X1, X2 = X2)
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  fit_aligned <- aligned_mcca(blocks, idx, N = N, ncomp = 2, ridge = 1e-6,
                              normalization = "None")
  fit_mcca <- mcca(blocks, ncomp = 2, ridge = 1e-6)

  S1 <- multivarious::scores(fit_aligned)
  S2 <- multivarious::scores(fit_mcca)
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))

  rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
  expect_lt(rel, 1e-10)
})

test_that("aligned_mcca custom alpha can drop a block cleanly", {
  set.seed(7)
  N <- 35
  X1 <- matrix(rnorm(N * 15), N, 15)
  X2 <- matrix(rnorm(N * 10), N, 10)
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  res_drop <- aligned_mcca(list(X1 = X1, X2 = X2), idx, N = N, ncomp = 2,
                           normalization = "custom", alpha = c(1, 0))
  expect_lt(max(abs(res_drop$partial_scores[[2]])), 1e-10)
  idx2 <- res_drop$block_indices[[2]]
  expect_lt(max(abs(res_drop$v[idx2, , drop = FALSE])), 1e-10)
})

test_that("aligned_mcca block_weights argument is deprecated but still works", {
  set.seed(7)
  N <- 35
  X1 <- matrix(rnorm(N * 15), N, 15)
  X2 <- matrix(rnorm(N * 10), N, 10)
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  res <- expect_warning(
    aligned_mcca(list(X1 = X1, X2 = X2), idx, N = N, ncomp = 2,
                 block_weights = c(1, 0)),
    regexp = "deprecated"
  )
  expect_equal(res$normalization, "custom")
  expect_equal(res$block_weights, c(1, 0))
  # Supplying both alpha and block_weights is an error.
  expect_error(
    aligned_mcca(list(X1 = X1, X2 = X2), idx, N = N, ncomp = 2,
                 alpha = c(1, 1), block_weights = c(1, 1))
  )
})

test_that("aligned_mcca uses eigen() fallback when ncomp == N", {
  set.seed(8)
  N <- 6
  X1 <- matrix(rnorm(N * 8), N, 8)
  X2 <- matrix(rnorm(N * 7), N, 7)
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  fit <- aligned_mcca(list(X1 = X1, X2 = X2), idx, N = N, ncomp = N, preproc = multivarious::pass())
  expect_equal(ncol(multivarious::scores(fit)), N)
})

test_that("aligned_mcca passes mathematical sanity checks on synthetic aligned low-rank data", {
  set.seed(10)
  N <- 80
  k <- 3
  Z <- matrix(rnorm(N * k), nrow = N, ncol = k)
  idx <- list(
    B1 = seq_len(N),                       # full coverage of reference rows
    B2 = sample.int(N, N, replace = FALSE), # permutation (row order differs)
    B3 = sample(rep(seq_len(N), 2))         # balanced duplicates (each row twice)
  )

  p_vec <- c(B1 = 25, B2 = 30, B3 = 20)
  blocks <- Map(function(p, idxk) {
    A <- matrix(rnorm(p * k), nrow = p, ncol = k)
    signal <- Z[idxk, , drop = FALSE] %*% t(A)
    signal + matrix(rnorm(length(idxk) * p, sd = 0.02), nrow = length(idxk), ncol = p)
  }, p_vec, idx)

  res <- aligned_mcca(blocks, idx, N = N, ncomp = k, ridge = 1e-6,
                      normalization = "None")

  S <- multivarious::scores(res)
  expect_equal(dim(S), c(N, k))

  XtX <- crossprod(S)
  expect_lt(max(abs(XtX - diag(diag(XtX)))), 1e-6)
  expect_equal(unname(diag(XtX)), unname(res$lambda[seq_len(k)]), tolerance = 1e-6)
  expect_equal(unname(res$sdev[seq_len(k)]^2), unname(res$lambda[seq_len(k)]), tolerance = 1e-10)

  # Sum of block partial scores mapped into reference space reproduces compromise scores
  S_sum <- Reduce(`+`, lapply(seq_along(res$partial_scores), function(k) {
    .expand_rowsum_to_N(res$partial_scores[[k]], res$row_index[[k]], res$N)
  }))
  expect_equal(unname(S_sum), unname(S), tolerance = 1e-6)

  # Recover latent space (up to rotation/sign) via canonical correlations
  cc <- cancor(scale(Z, center = TRUE, scale = FALSE),
               scale(S, center = TRUE, scale = FALSE))$cor
  expect_gt(min(cc[seq_len(k)]), 0.98)

  # Internal consistency between v and canonical_weights
  for (k in seq_along(res$block_indices)) {
    idx_feat <- res$block_indices[[k]]
    V_block <- res$v[idx_feat, , drop = FALSE]
    W_raw <- res$canonical_weights[[k]]
    expect_equal(
      V_block %*% diag(res$sdev),
      res$block_weights[k] * W_raw,
      tolerance = 1e-6
    )
  }
})

test_that("aligned_mcca handles p > n blocks without error", {
  set.seed(3)
  N <- 25
  X1 <- matrix(rnorm(10 * 60), 10, 60)
  X2 <- matrix(rnorm(12 * 55), 12, 55)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit <- aligned_mcca(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), N = N, ncomp = 3)
  expect_s3_class(fit, "aligned_mcca")
  expect_equal(dim(multivarious::scores(fit)), c(N, 3))
})

test_that("aligned_mcca ridge=0 warns and still returns a fit for rank-deficient blocks", {
  set.seed(5)
  N <- 20
  X_singular <- matrix(rnorm(N * 5), nrow = N, ncol = 5)   # K is singular (rank <= 5)
  X_fullrank <- matrix(rnorm(N * 30), nrow = N, ncol = 30) # typically full row rank
  idx <- list(B1 = seq_len(N), B2 = seq_len(N))

  expect_warning(
    res <- aligned_mcca(list(B1 = X_singular, B2 = X_fullrank), idx, N = N, ncomp = 2, ridge = 0),
    regexp = "ridge=0 yields a singular system"
  )
  expect_s3_class(res, "aligned_mcca")
})

test_that("anchored_mcca with normalization='None' matches aligned_mcca with Y included", {
  set.seed(4)
  N <- 30
  Y <- matrix(rnorm(N * 5), N, 5)
  X1 <- matrix(rnorm(20 * 10), 20, 10)
  X2 <- matrix(rnorm(15 * 8), 15, 8)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit_anchor <- anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                              ncomp = 2, normalization = "None")
  expect_s3_class(fit_anchor, "anchored_mcca")
  expect_s3_class(fit_anchor, "aligned_mcca")
  expect_equal(fit_anchor$normalization, "None")

  fit_aligned <- aligned_mcca(
    X = c(list(Y = Y), list(X1 = X1, X2 = X2)),
    row_index = c(list(Y = seq_len(N)), list(X1 = idx1, X2 = idx2)),
    N = N,
    ncomp = 2,
    normalization = "None"
  )

  S1 <- multivarious::scores(fit_anchor)
  S2 <- multivarious::scores(fit_aligned)
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
  expect_lt(rel, 1e-10)
})

test_that("anchored_mcca default normalization is 'MFA' and matches anchored_mfa interface", {
  set.seed(4)
  N <- 30
  Y <- matrix(rnorm(N * 5), N, 5)
  X1 <- matrix(rnorm(20 * 10), 20, 10)
  X2 <- matrix(rnorm(15 * 8), 15, 8)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit <- anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
  expect_equal(fit$normalization, "MFA")
  # Interface parity: same core args should be accepted as anchored_mfa
  shared_args <- c("Y", "X", "row_index", "preproc", "ncomp", "normalization",
                   "alpha", "ridge", "max_iter", "tol", "verbose", "use_future")
  fn_mcca <- names(formals(anchored_mcca))
  fn_mfa <- names(formals(anchored_mfa))
  expect_true(all(shared_args %in% fn_mcca))
  expect_true(all(shared_args %in% fn_mfa))
})

test_that("anchored_mcca adapts alpha for Y and validates names/lengths", {
  set.seed(6)
  N <- 30
  Y <- matrix(rnorm(N * 5), N, 5)
  X1 <- matrix(rnorm(20 * 10), 20, 10)
  X2 <- matrix(rnorm(15 * 8), 15, 8)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  # X-only alpha -> prepend Y weight 1
  fit_xonly <- anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                             ncomp = 2, normalization = "custom", alpha = c(0.5, 2))
  expect_equal(fit_xonly$block_weights, c(1, 0.5, 2))
  expect_equal(unname(fit_xonly$alpha_blocks), c(1, 0.5, 2))

  # Named X-only alpha -> still prepend Y=1 and order as (Y, X1, X2)
  fit_named <- anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                             ncomp = 2, normalization = "custom",
                             alpha = c(X1 = 0.5, X2 = 2))
  expect_equal(fit_named$block_weights, c(1, 0.5, 2))

  # Wrong length
  expect_error(
    anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                  ncomp = 2, normalization = "custom", alpha = 1)
  )

  # Missing X block name
  expect_error(
    anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                  ncomp = 2, normalization = "custom", alpha = c(X1 = 1))
  )

  # 'custom' without alpha is an error
  expect_error(
    anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                  ncomp = 2, normalization = "custom")
  )
})

test_that("anchored_mcca block_weights argument is deprecated but still works", {
  set.seed(6)
  N <- 30
  Y <- matrix(rnorm(N * 5), N, 5)
  X1 <- matrix(rnorm(20 * 10), 20, 10)
  X2 <- matrix(rnorm(15 * 8), 15, 8)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit <- expect_warning(
    anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                  ncomp = 2, block_weights = c(0.5, 2)),
    regexp = "deprecated"
  )
  expect_equal(fit$block_weights, c(1, 0.5, 2))
  expect_equal(fit$normalization, "custom")

  # Supplying both alpha and block_weights is an error.
  expect_error(
    anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                  ncomp = 2, alpha = c(0.5, 2), block_weights = c(0.5, 2))
  )
})

test_that("anchored_mcca assigns default X names and validates row_index values", {
  set.seed(9)
  N <- 25
  Y <- matrix(rnorm(N * 4), N, 4)
  X1 <- matrix(rnorm(10 * 6), 10, 6)
  X2 <- matrix(rnorm(12 * 5), 12, 5)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit <- anchored_mcca(Y, list(X1, X2), list(idx1, idx2), ncomp = 2)
  expect_equal(fit$names, c("Y", "X1", "X2"))

  idx_bad <- idx2
  idx_bad[1] <- NA_integer_
  expect_error(anchored_mcca(Y, list(X1, X2), list(idx1, idx_bad), ncomp = 2))

  idx_oob <- idx2
  idx_oob[1] <- N + 1L
  expect_error(anchored_mcca(Y, list(X1, X2), list(idx1, idx_oob), ncomp = 2))
})

test_that("anchored_mcca 'balanced' mode reduces single-block domination", {
  set.seed(33)
  N <- 40
  Y <- matrix(rnorm(N * 5), N, 5)
  X1 <- matrix(rnorm(25 * 10), 25, 10)
  # X2 has 10x scale so the unweighted eigenvalue landscape is dominated by X2.
  X2 <- 10 * matrix(rnorm(20 * 8), 20, 8)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit_none <- anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                            ncomp = 2, normalization = "None")
  fit_bal  <- anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                            ncomp = 2, normalization = "balanced", max_iter = 80)

  imbal <- function(fit) {
    apply(fit$block_contribs, 1, function(r) {
      pos <- r[r > 0]; if (length(pos) > 1) max(pos) / min(pos) else 1
    })
  }
  imb_none <- imbal(fit_none)
  imb_bal  <- imbal(fit_bal)

  # Balanced mode should give a strictly smaller per-component imbalance ratio.
  expect_true(all(imb_bal <= imb_none + 1e-8))
  # And a substantially more balanced worst-case.
  expect_lt(max(imb_bal), max(imb_none))

  # Diagnostics populated.
  expect_equal(fit_bal$normalization, "balanced")
  expect_true(!is.null(fit_bal$alpha_per_component))
  expect_equal(dim(fit_bal$alpha_per_component), c(2L, 3L))
  expect_true(all(is.finite(fit_bal$alpha_blocks)))
  expect_length(fit_bal$balance_converged_per_comp, 2L)
})

test_that("anchored_mcca balanced diagnostics satisfy contribution algebra", {
  set.seed(34)
  N <- 36
  Y <- matrix(rnorm(N * 5), N, 5)
  X1 <- matrix(rnorm(28 * 7), 28, 7)
  X2 <- matrix(rnorm(30 * 6), 30, 6)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit <- anchored_mcca(
    Y,
    list(X1 = X1, X2 = X2),
    list(X1 = idx1, X2 = idx2),
    ncomp = 2,
    normalization = "balanced",
    max_iter = 80
  )

  expect_equal(dim(fit$weighted_block_contribs), dim(fit$block_contribs))
  expect_equal(dim(fit$block_contrib_fraction), dim(fit$block_contribs))
  expect_equal(
    unname(fit$weighted_block_contribs),
    unname(fit$block_contribs * fit$alpha_per_component),
    tolerance = 1e-10
  )
  expect_equal(unname(rowSums(fit$block_contrib_fraction)), rep(1, nrow(fit$block_contrib_fraction)), tolerance = 1e-10)
  expect_true(all(is.finite(fit$weighted_block_contribs)))
})

test_that("anchored_mcca balanced optimization trace is monotone per component", {
  set.seed(35)
  N <- 45
  Z <- matrix(rnorm(N * 2), N, 2)
  Y <- Z %*% matrix(rnorm(2 * 5), 2, 5) + matrix(rnorm(N * 5, sd = 0.05), N, 5)
  X1 <- Z %*% matrix(rnorm(2 * 8), 2, 8) + matrix(rnorm(N * 8, sd = 0.05), N, 8)
  X2 <- Z %*% matrix(rnorm(2 * 7), 2, 7) + matrix(rnorm(N * 7, sd = 0.05), N, 7)

  fit <- anchored_mcca(
    Y,
    list(X1 = X1, X2 = X2),
    list(X1 = seq_len(N), X2 = seq_len(N)),
    ncomp = 2,
    normalization = "balanced",
    ridge = 1e-6,
    max_iter = 80
  )

  for (trace_j in fit$balance_trace) {
    objectives <- vapply(trace_j, `[[`, numeric(1), "objective")
    expect_true(all(is.finite(objectives)))
    expect_true(all(diff(objectives) >= -1e-10))
  }
})

test_that("anchored_mcca exposes a refit contract without retaining fit-time work objects", {
  set.seed(37)
  N <- 24
  Y <- matrix(rnorm(N * 4), N, 4)
  X1 <- matrix(rnorm(18 * 5), 18, 5)
  X2 <- matrix(rnorm(20 * 6), 20, 6)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit <- anchored_mcca(
    Y,
    list(X1 = X1, X2 = X2),
    list(X1 = idx1, X2 = idx2),
    ncomp = 2,
    normalization = "balanced",
    max_iter = 5
  )

  expect_equal(fit$task, "association")
  expect_equal(fit$fit_spec$method, "anchored_mcca")
  expect_true(fit$fit_spec$refit_supported)
  expect_setequal(fit$oos_types, "scores")

  refit <- fit$fit_spec$refit$fit_fn(fit$fit_spec$refit$data)
  expect_s3_class(refit, "anchored_mcca")
  expect_equal(dim(refit$s), dim(fit$s))
  expect_true(all(is.finite(refit$sdev)))

  callback_env_names <- ls(environment(fit$fit_spec$refit$fit_fn), all.names = TRUE)
  expect_false(any(c("X", "Xp", "M_list", "block_fits", "fit", "out") %in% callback_env_names))
  expect_true(all(c("preproc", "ncomp", "normalization", "fit_dots") %in% callback_env_names))
})

test_that("anchored_mcca 'balanced' prefers shared components over block-specific ones", {
  # Simulate: 1 shared latent component, 2 block-specific latents per X block,
  # scaled/amplified so the unweighted MAXVAR eigenproblem is dominated by X2's
  # block-specific signal. The anchor Y is deliberately noisy so it cannot
  # single-handedly drag the top eigenvectors toward the shared direction.
  make_hard <- function(seed, N = 100, noise_sd = 0.5) {
    set.seed(seed)
    Z_shared <- matrix(rnorm(N), N, 1)
    Z_X1_spec <- matrix(rnorm(N * 2), N, 2)
    Z_X2_spec <- matrix(rnorm(N * 2), N, 2)

    BY <- matrix(rnorm(1 * 4), 1, 4)
    Y <- 0.5 * Z_shared %*% BY + matrix(rnorm(N * 4, sd = noise_sd), N, 4)

    A1  <- matrix(rnorm(1 * 10), 1, 10)
    A1s <- matrix(rnorm(2 * 10), 2, 10)
    X1 <- Z_shared %*% A1 + Z_X1_spec %*% A1s +
          matrix(rnorm(N * 10, sd = noise_sd), N, 10)

    A2  <- matrix(rnorm(1 * 8), 1, 8)
    A2s <- matrix(rnorm(2 * 8), 2, 8)
    X2 <- 20 * (Z_shared %*% A2 + 5 * (Z_X2_spec %*% A2s) +
                matrix(rnorm(N * 8, sd = noise_sd), N, 8))

    list(Y = Y, X1 = X1, X2 = X2,
         Z_shared = Z_shared, Z_X1_spec = Z_X1_spec, Z_X2_spec = Z_X2_spec,
         idx1 = seq_len(N), idx2 = seq_len(N))
  }

  # "Breadth" per component: min/max block contribution.
  # ~1 → every block contributes (shared-like). ~0 → dominated by one block.
  breadth <- function(bc) {
    apply(bc, 1, function(r) {
      pos <- r[r > 0]
      if (length(pos) > 1) min(pos) / max(pos) else 1
    })
  }

  # Run across several seeds; claim the averages over seeds, not a single draw.
  n_seeds <- 5
  res <- lapply(seq_len(n_seeds) + 100L, function(seed) {
    sim <- make_hard(seed)
    f_mfa <- anchored_mcca(sim$Y, list(X1 = sim$X1, X2 = sim$X2),
                           list(X1 = sim$idx1, X2 = sim$idx2),
                           ncomp = 3, normalization = "MFA", max_iter = 80)
    f_bal <- anchored_mcca(sim$Y, list(X1 = sim$X1, X2 = sim$X2),
                           list(X1 = sim$idx1, X2 = sim$idx2),
                           ncomp = 3, normalization = "balanced", max_iter = 80)
    list(
      b_mfa = breadth(f_mfa$block_contribs),
      b_bal = breadth(f_bal$block_contribs),
      cc_mfa = abs(as.numeric(cor(multivarious::scores(f_mfa)[, 1], sim$Z_shared))),
      cc_bal = abs(as.numeric(cor(multivarious::scores(f_bal)[, 1], sim$Z_shared)))
    )
  })

  # Core claim 1: balanced's WORST per-component breadth is substantially higher
  # than MFA's worst — i.e., balanced never extracts a block-specific component,
  # while MFA does. Under this sim MFA's later comps collapse onto Y alone.
  min_breadth_mfa <- vapply(res, function(r) min(r$b_mfa), numeric(1))
  min_breadth_bal <- vapply(res, function(r) min(r$b_bal), numeric(1))
  expect_true(all(min_breadth_bal > min_breadth_mfa))
  expect_gt(mean(min_breadth_bal), 0.5)
  expect_lt(mean(min_breadth_mfa), 0.15)

  # Core claim 2: leading-component alignment with the shared latent is at
  # least as good under balanced as under MFA.
  cc_mfa <- vapply(res, `[[`, numeric(1), "cc_mfa")
  cc_bal <- vapply(res, `[[`, numeric(1), "cc_bal")
  expect_gte(mean(cc_bal), mean(cc_mfa) - 1e-6)
  expect_gt(mean(cc_bal), 0.9)

  # Core claim 3: balanced's leading shared-correlation is strictly higher on
  # average (the net result of rebalancing when Y is noisy).
  expect_gt(mean(cc_bal - cc_mfa), 0)
})

test_that("anchored_mcca handles p > n blocks under every normalization mode", {
  set.seed(7)
  N <- 25
  Y <- matrix(rnorm(N * 5), N, 5)
  # Both X blocks have p >> n.
  X1 <- matrix(rnorm(10 * 200), 10, 200)
  X2 <- matrix(rnorm(12 * 150), 12, 150)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  for (mode in c("None", "MFA", "balanced", "custom")) {
    args <- list(Y = Y, X = list(X1 = X1, X2 = X2),
                 row_index = list(X1 = idx1, X2 = idx2),
                 ncomp = 3, normalization = mode, max_iter = 30)
    if (mode == "custom") args$alpha <- c(1, 1, 1)
    fit <- do.call(anchored_mcca, args)

    expect_s3_class(fit, "anchored_mcca")
    expect_equal(dim(multivarious::scores(fit)), c(N, 3L))
    # Concatenated loadings match total feature count.
    expect_equal(nrow(fit$v), ncol(Y) + ncol(X1) + ncol(X2))
    # Per-block partial scores are n_k × ncomp, not p_k × ncomp.
    expect_equal(dim(fit$partial_scores$Y),  c(nrow(Y),  3L))
    expect_equal(dim(fit$partial_scores$X1), c(nrow(X1), 3L))
    expect_equal(dim(fit$partial_scores$X2), c(nrow(X2), 3L))
    # All eigenvalues are finite and non-negative.
    expect_true(all(is.finite(fit$lambda)) && all(fit$lambda >= 0))
  }
})

test_that("anchored_mcca handles extreme p >> n and rank-deficient blocks with ridge", {
  set.seed(42)
  N <- 20
  Y <- matrix(rnorm(N * 3), N, 3)
  # Extreme ratios: p/n = 125 and 80.
  X1 <- matrix(rnorm(8 * 1000), 8, 1000)
  X2 <- matrix(rnorm(10 * 800), 10, 800)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  # All three "production" modes should complete without error.
  for (mode in c("MFA", "balanced")) {
    fit <- anchored_mcca(Y, list(X1 = X1, X2 = X2),
                         list(X1 = idx1, X2 = idx2),
                         ncomp = 2, normalization = mode, max_iter = 30)
    expect_s3_class(fit, "anchored_mcca")
    expect_equal(dim(multivarious::scores(fit)), c(N, 2L))
  }

  # ridge=0 on a rank-deficient block (p < n, so K is low-rank) should warn
  # and still return a valid fit via the automatic ridge floor.
  X_singular <- matrix(rnorm(N * 5), N, 5)    # rank ≤ 5 < N after centering
  X_wide     <- matrix(rnorm(N * 100), N, 100)
  expect_warning(
    fit <- anchored_mcca(Y, list(Xs = X_singular, Xw = X_wide),
                         list(Xs = seq_len(N), Xw = seq_len(N)),
                         ncomp = 2, normalization = "None", ridge = 0),
    regexp = "ridge=0 yields a singular system"
  )
  expect_s3_class(fit, "anchored_mcca")
  expect_equal(dim(multivarious::scores(fit)), c(N, 2L))
})

test_that("anchored_mcca accepts max_iter/tol/verbose/use_future for interface parity", {
  set.seed(12)
  N <- 20
  Y <- matrix(rnorm(N * 3), N, 3)
  X1 <- matrix(rnorm(15 * 6), 15, 6)
  X2 <- matrix(rnorm(18 * 5), 18, 5)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  # These arguments should all be accepted without error.
  fit <- anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                       ncomp = 2, normalization = "balanced",
                       max_iter = 10L, tol = 1e-3, verbose = FALSE,
                       use_future = FALSE)
  expect_s3_class(fit, "anchored_mcca")
})

test_that("aligned_mcca use_future path matches serial", {
  skip_if_not_installed("future")
  skip_if_not_installed("furrr")

  op <- future::plan()
  on.exit(future::plan(op), add = TRUE)
  future::plan(future::sequential)

  set.seed(22)
  N <- 40
  X1 <- matrix(rnorm(N * 10), N, 10)
  X2 <- matrix(rnorm(N * 12), N, 12)
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  res_serial <- aligned_mcca(list(X1 = X1, X2 = X2), idx, N = N, ncomp = 2, use_future = FALSE)
  res_future <- aligned_mcca(list(X1 = X1, X2 = X2), idx, N = N, ncomp = 2, use_future = TRUE)

  S1 <- multivarious::scores(res_serial)
  S2 <- multivarious::scores(res_future)
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
  expect_lt(rel, 1e-12)
})

test_that("aligned_mcca is invariant to relabeling the reference row space", {
  set.seed(32)
  sim <- .sim_aligned_mcca(
    N = 70,
    k = 2,
    p_vec = c(16, 18, 14),
    n_vec = c(50, 48, 55),
    noise_sd = 0.03
  )

  fit <- aligned_mcca(sim$blocks, sim$idx, N = 70, ncomp = 2, ridge = 1e-6)

  perm <- sample.int(70)
  idx_perm <- lapply(sim$idx, .permute_reference_labels, perm = perm)
  fit_perm <- aligned_mcca(sim$blocks, idx_perm, N = 70, ncomp = 2, ridge = 1e-6)

  S1 <- multivarious::scores(fit)
  S2 <- multivarious::scores(fit_perm)[perm, , drop = FALSE]
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  rel <- norm(P1 - P2, type = "F") / (norm(P1, type = "F") + 1e-12)
  expect_lt(rel, 1e-10)
})

test_that("aligned_mcca loses recovery quality when row_index is randomized", {
  set.seed(33)
  sim <- .sim_aligned_mcca(
    N = 75,
    k = 2,
    p_vec = c(15, 17, 14),
    n_vec = c(55, 52, 58),
    noise_sd = 0.02
  )

  fit_true <- aligned_mcca(sim$blocks, sim$idx, N = 75, ncomp = 2, ridge = 1e-6)
  idx_wrong <- lapply(sim$idx, sample)
  fit_wrong <- aligned_mcca(sim$blocks, idx_wrong, N = 75, ncomp = 2, ridge = 1e-6)

  cc_true <- cancor(
    scale(sim$Z, center = TRUE, scale = FALSE),
    scale(multivarious::scores(fit_true), center = TRUE, scale = FALSE)
  )$cor
  cc_wrong <- cancor(
    scale(sim$Z, center = TRUE, scale = FALSE),
    scale(multivarious::scores(fit_wrong), center = TRUE, scale = FALSE)
  )$cor

  expect_gt(mean(cc_true[1:2]), mean(cc_wrong[1:2]) + 0.2)
  expect_gt(mean(cc_true[1:2]), 0.65)
})
