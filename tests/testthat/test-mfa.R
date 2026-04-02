library(testthat)
library(muscal)

# Basic MFA run should produce a valid object

.sim_mfa_blocks <- function(n, k, p_vec, noise_sd = 0) {
  Z <- scale(matrix(rnorm(n * k), nrow = n, ncol = k), center = TRUE, scale = FALSE)
  blocks <- lapply(p_vec, function(p) {
    A <- matrix(rnorm(p * k), nrow = p, ncol = k)
    signal <- Z %*% t(A)
    signal + matrix(rnorm(n * p, sd = noise_sd), nrow = n, ncol = p)
  })
  list(blocks = blocks, Z = Z)
}

.subspace_rel_diff <- function(S1, S2) {
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
}

test_that("mfa.list produces expected structure", {
  set.seed(1)
  blocks <- replicate(3, matrix(rnorm(20), nrow = 5), simplify = FALSE)
  res <- mfa(blocks, ncomp = 2)
  expect_s3_class(res, "mfa")
  expect_equal(length(res$block_indices), 3)
  ncols <- sapply(blocks, ncol)
  expect_equal(res$block_indices[[1]], 1:ncols[1])
  expect_equal(res$block_indices[[2]], (ncols[1]+1):(ncols[1]+ncols[2]))
  expect_equal(res$block_indices[[3]], (ncols[1]+ncols[2]+1):sum(ncols))
  expect_length(res$alpha, 3)
  expect_true(all(res$alpha > 0))
  expect_equal(multivarious::ncomp(res), 2)
})

# Error conditions in mfa.multiblock

test_that("mfa.multiblock validates input", {
  x1 <- matrix(rnorm(10), nrow = 5)
  expect_error(mfa(list(x1)))

  x2 <- matrix(rnorm(12), nrow = 6)
  expect_error(mfa(list(x1, x2)))

  expect_error(mfa(list(x1, x1), normalization = "custom"))
})

# Normalization factor computation for MFA mode

test_that("normalization_factors MFA matches svd", {
  X1 <- matrix(1:4, nrow = 2)
  X2 <- matrix(2:5, nrow = 2)
  nf <- muscal:::normalization_factors(list(X1, X2), type = "MFA")
  expect_equal(nf[1], 1/(svd(X1)$d[1]^2))
  expect_equal(nf[2], 1/(svd(X2)$d[1]^2))
})

test_that("mfa exposes a standard out-of-sample contract", {
  set.seed(4)
  blocks <- list(
    X1 = matrix(rnorm(60), nrow = 10),
    X2 = matrix(rnorm(50), nrow = 10)
  )
  fit <- mfa(blocks, ncomp = 2)

  expect_equal(fit$task, "reconstruction")
  expect_equal(fit$fit_spec$method, "mfa")
  expect_true(fit$fit_spec$refit_supported)
  expect_true(is.list(fit$fit_spec$refit))
  expect_setequal(fit$oos_types, c("scores", "reconstruction"))

  Xnew <- do.call(cbind, blocks) + matrix(rnorm(110, sd = 0.05), nrow = 10)
  scores_new <- predict(fit, Xnew, type = "scores")
  recon_new <- predict(fit, Xnew, type = "reconstruction")

  expect_equal(scores_new, multivarious::project(fit, Xnew))
  expect_equal(dim(scores_new), c(nrow(Xnew), multivarious::ncomp(fit)))
  expect_equal(dim(recon_new), dim(Xnew))
  expect_equal(predict(fit, type = "scores"), multivarious::scores(fit))
  expect_error(predict(fit, type = "reconstruction"), "`new_data` must be supplied")
})

test_that("mfa refit contract can refit the stored training data", {
  set.seed(41)
  blocks <- list(
    X1 = matrix(rnorm(48), nrow = 12),
    X2 = matrix(rnorm(36), nrow = 12)
  )
  fit <- mfa(blocks, ncomp = 2)

  refit <- fit$fit_spec$refit$fit_fn(fit$fit_spec$refit$data)

  expect_s3_class(refit, "mfa")
  expect_equal(dim(refit$v), dim(fit$v))
  expect_equal(length(refit$block_indices), length(fit$block_indices))
  expect_true(all(is.finite(refit$sdev)))
})

test_that("mfa exactly recovers noiseless low-rank structure and reconstruction improves with ncomp", {
  set.seed(43)
  sim <- .sim_mfa_blocks(n = 60, k = 3, p_vec = c(10, 12, 8), noise_sd = 0)
  X_concat <- do.call(cbind, sim$blocks)

  fits <- lapply(1:3, function(k) {
    suppressMessages(
      mfa(
        sim$blocks,
        ncomp = k,
        preproc = multivarious::pass(),
        normalization = "None"
      )
    )
  })

  mse <- vapply(fits, function(fit) {
    Xhat <- predict(fit, X_concat, type = "reconstruction")
    mean((Xhat - X_concat)^2)
  }, numeric(1))

  expect_true(all(diff(mse) <= 1e-10))
  expect_lt(mse[[3]], 1e-10)

  S_fit <- multivarious::scores(fits[[3]])
  P_fit <- S_fit %*% solve(crossprod(S_fit), t(S_fit))
  P_true <- sim$Z %*% solve(crossprod(sim$Z), t(sim$Z))
  rel <- norm(P_fit - P_true, type = "F") / (norm(P_true, type = "F") + 1e-12)
  expect_lt(rel, 1e-8)
})

test_that("mfa latent recovery degrades as noise increases", {
  set.seed(44)
  low_noise <- .sim_mfa_blocks(n = 70, k = 2, p_vec = c(9, 11, 10), noise_sd = 0.02)
  high_noise <- .sim_mfa_blocks(n = 70, k = 2, p_vec = c(9, 11, 10), noise_sd = 0.35)

  fit_low <- suppressMessages(
    mfa(low_noise$blocks, ncomp = 2, preproc = multivarious::pass(), normalization = "None")
  )
  fit_high <- suppressMessages(
    mfa(high_noise$blocks, ncomp = 2, preproc = multivarious::pass(), normalization = "None")
  )

  rel_low <- .subspace_rel_diff(multivarious::scores(fit_low), low_noise$Z)
  rel_high <- .subspace_rel_diff(multivarious::scores(fit_high), high_noise$Z)

  expect_lt(rel_low, rel_high)
  expect_lt(rel_low, 0.15)
})
