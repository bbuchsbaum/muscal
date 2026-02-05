library(testthat)
library(muscal)

.sim_mcca_blocks <- function(n, k, p_vec, noise_sd = 0.05) {
  Z <- matrix(rnorm(n * k), nrow = n, ncol = k)
  blocks <- lapply(p_vec, function(p) {
    A <- matrix(rnorm(p * k), nrow = p, ncol = k)
    signal <- Z %*% t(A)
    signal + matrix(rnorm(n * p, sd = noise_sd), nrow = n, ncol = p)
  })
  list(blocks = blocks, Z = Z)
}

test_that("mcca.list produces expected structure", {
  set.seed(1)
  blocks <- replicate(3, matrix(rnorm(30), nrow = 5), simplify = FALSE)
  res <- mcca(blocks, ncomp = 2)

  expect_s3_class(res, "mcca")
  expect_equal(length(res$block_indices), 3)

  ncols <- sapply(blocks, ncol)
  expect_equal(res$block_indices[[1]], 1:ncols[1])
  expect_equal(res$block_indices[[2]], (ncols[1] + 1):(ncols[1] + ncols[2]))
  expect_equal(res$block_indices[[3]], (ncols[1] + ncols[2] + 1):sum(ncols))

  expect_equal(multivarious::ncomp(res), 2)
  expect_equal(dim(multivarious::scores(res)), c(5, 2))
  expect_equal(length(res$partial_scores), 3)
  expect_equal(dim(res$partial_scores[[1]]), c(5, 2))
})

test_that("mcca handles p > n blocks without error", {
  set.seed(2)
  blocks <- replicate(3, matrix(rnorm(5 * 30), nrow = 5), simplify = FALSE)
  res <- mcca(blocks, ncomp = 2)
  expect_s3_class(res, "mcca")

  X_concat <- do.call(cbind, blocks)
  proj <- project(res, X_concat)
  expect_true(is.matrix(proj))
  expect_equal(dim(proj), c(5, 2))

  # projecting the training data should reproduce stored scores (up to sign)
  s_fit <- multivarious::scores(res)
  expect_equal(abs(cor(proj[, 1], s_fit[, 1])), 1, tolerance = 1e-6)
})

test_that("mcca passes basic mathematical sanity checks on synthetic low-rank data", {
  set.seed(10)
  sim <- .sim_mcca_blocks(n = 80, k = 3, p_vec = c(20, 25, 30), noise_sd = 0.03)
  res <- mcca(sim$blocks, ncomp = 3, ridge = 1e-6)

  S <- multivarious::scores(res)
  expect_equal(dim(S), c(80, 3))

  # Orthonormal eigenvectors imply crossprod(scores) == diag(lambda)
  XtX <- crossprod(S)
  expect_lt(max(abs(XtX - diag(diag(XtX)))), 1e-6)
  expect_equal(
    unname(diag(XtX)),
    unname(res$lambda[seq_len(3)]),
    tolerance = 1e-6
  )
  expect_equal(
    unname(res$sdev[seq_len(3)]^2),
    unname(res$lambda[seq_len(3)]),
    tolerance = 1e-10
  )

  # Sum of block partial scores reproduces the compromise scores
  S_sum <- Reduce(`+`, res$partial_scores)
  expect_equal(unname(S_sum), unname(S), tolerance = 1e-6)

  # Recover latent space (up to rotation/sign) via canonical correlations
  cc <- cancor(scale(sim$Z, center = TRUE, scale = FALSE),
               scale(S, center = TRUE, scale = FALSE))$cor
  expect_gt(min(cc[seq_len(3)]), 0.9)

  # Internal consistency between v and canonical_weights
  for (i in seq_along(res$block_indices)) {
    idx <- res$block_indices[[i]]
    V_block <- res$v[idx, , drop = FALSE]
    W_raw <- res$canonical_weights[[i]]
    expect_equal(
      V_block %*% diag(res$sdev),
      res$block_weights[i] * W_raw,
      tolerance = 1e-6
    )
  }
})

test_that("mcca is invariant to block rescaling (with scaled ridge)", {
  set.seed(11)
  sim <- .sim_mcca_blocks(n = 60, k = 2, p_vec = c(25, 25, 25), noise_sd = 0.05)
  blocks1 <- sim$blocks
  blocks2 <- sim$blocks
  blocks2[[1]] <- 10 * blocks2[[1]]

  res1 <- mcca(blocks1, ncomp = 2, ridge = 1e-6)
  res2 <- mcca(blocks2, ncomp = 2, ridge = 1e-6)

  S1 <- scale(multivarious::scores(res1), center = TRUE, scale = FALSE)
  S2 <- scale(multivarious::scores(res2), center = TRUE, scale = FALSE)
  cc <- cancor(S1, S2)$cor
  expect_gt(min(cc[seq_len(2)]), 0.99)
})

test_that("mcca block_weights can drop a block cleanly", {
  set.seed(12)
  sim <- .sim_mcca_blocks(n = 50, k = 2, p_vec = c(20, 20, 20), noise_sd = 0.04)
  blocks <- sim$blocks

  res_drop <- mcca(blocks, ncomp = 2, ridge = 1e-6, block_weights = c(1, 0, 1))
  expect_lt(max(abs(res_drop$partial_scores[[2]])), 1e-10)
  idx2 <- res_drop$block_indices[[2]]
  expect_lt(max(abs(res_drop$v[idx2, , drop = FALSE])), 1e-10)

  res_reduced <- mcca(blocks[c(1, 3)], ncomp = 2, ridge = 1e-6)
  S_drop <- scale(multivarious::scores(res_drop), center = TRUE, scale = FALSE)
  S_red <- scale(multivarious::scores(res_reduced), center = TRUE, scale = FALSE)
  cc <- cancor(S_drop, S_red)$cor
  expect_gt(min(cc[seq_len(2)]), 0.99)
})

test_that("mcca ridge=0 warns and still returns a fit for centered blocks", {
  set.seed(13)
  blocks <- replicate(3, matrix(rnorm(20 * 30), nrow = 20), simplify = FALSE)
  # Force rank deficiency deterministically so K is not SPD when ridge=0
  blocks <- lapply(blocks, function(X) {
    X[nrow(X), ] <- X[1, ]
    X
  })
  expect_warning(
    res <- mcca(blocks, ncomp = 2, ridge = 0),
    regexp = "ridge=0 yields a singular system"
  )
  expect_s3_class(res, "mcca")
})

test_that("mcca use_future path matches serial", {
  skip_if_not_installed("future")
  skip_if_not_installed("furrr")

  op <- future::plan()
  on.exit(future::plan(op), add = TRUE)
  future::plan(future::sequential)

  set.seed(21)
  blocks <- replicate(3, matrix(rnorm(60 * 15), nrow = 60), simplify = FALSE)

  res_serial <- mcca(blocks, ncomp = 2, ridge = 1e-6, use_future = FALSE)
  res_future <- mcca(blocks, ncomp = 2, ridge = 1e-6, use_future = TRUE)

  S1 <- multivarious::scores(res_serial)
  S2 <- multivarious::scores(res_future)
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
  expect_lt(rel, 1e-12)
})

test_that("mcca covers edge cases (zero/NA blocks, ridge=0 without centering, eigen fallback)", {
  set.seed(30)
  n <- 10

  X_zero <- matrix(0, nrow = n, ncol = 5)
  X_dense <- matrix(rnorm(n * 6), nrow = n, ncol = 6)

  fit <- mcca(list(X_zero = X_zero, X_dense = X_dense, X_more = matrix(rnorm(n * 7), n, 7)),
              preproc = multivarious::pass(), ncomp = 2, ridge = 1e-6)
  expect_s3_class(fit, "mcca")
  expect_equal(dim(multivarious::scores(fit)), c(n, 2))

  X_bad <- X_dense
  X_bad[1, 1] <- NA_real_
  expect_error(
    mcca(list(B1 = X_bad, B2 = X_dense), preproc = multivarious::pass(), ncomp = 2),
    regexp = "Non-finite values"
  )

  # With p >= n and no centering, ridge=0 should be feasible without warnings.
  X_full1 <- matrix(rnorm(n * 15), n, 15)
  X_full2 <- matrix(rnorm(n * 12), n, 12)
  expect_warning(
    fit0 <- mcca(list(B1 = X_full1, B2 = X_full2),
                 preproc = multivarious::pass(), ncomp = 2, ridge = 0),
    regexp = NA
  )
  expect_s3_class(fit0, "mcca")

  expect_error(mcca(list(B1 = X_dense, B2 = matrix(rnorm(n * 8), n, 8)), ncomp = 1.5))

  # ncomp == n triggers eigen() fallback (skips RSpectra::eigs_sym)
  n2 <- 5
  blocks2 <- list(
    B1 = matrix(rnorm(n2 * 6), n2, 6),
    B2 = matrix(rnorm(n2 * 7), n2, 7)
  )
  fit_full <- mcca(blocks2, ncomp = n2, ridge = 1e-6, preproc = multivarious::pass())
  expect_equal(ncol(multivarious::scores(fit_full)), n2)
})
