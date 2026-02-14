library(testthat)
library(muscal)

.sim_ipca_blocks <- function(n, r, p_vec, noise_sd = 0.05) {
  Z <- matrix(rnorm(n * r), nrow = n, ncol = r)
  blocks <- lapply(p_vec, function(p) {
    A <- matrix(rnorm(p * r), nrow = p, ncol = r)
    Z %*% t(A) + matrix(rnorm(n * p, sd = noise_sd), nrow = n, ncol = p)
  })
  names(blocks) <- paste0("B", seq_along(blocks))
  blocks
}

test_that("ipca returns expected structure", {
  set.seed(1)
  blocks <- .sim_ipca_blocks(n = 60, r = 3, p_vec = c(25, 30, 20), noise_sd = 0.03)
  fit <- ipca(blocks, ncomp = 3, lambda = 1, method = "gram", max_iter = 60, tol = 1e-6)

  expect_s3_class(fit, "ipca")
  expect_s3_class(fit, "multiblock_biprojector")
  expect_equal(dim(multivarious::scores(fit)), c(60, 3))
  expect_equal(dim(fit$v), c(sum(vapply(blocks, ncol, integer(1))), 3))
  expect_equal(length(fit$block_indices), 3)
  expect_equal(length(fit$Delta_leading_eigenvalues), 3)
  expect_true(all(is.finite(fit$Sigma_eigenvalues)))
  expect_equal(mean(fit$Sigma_eigenvalues), 1, tolerance = 1e-6)
  expect_true(fit$iter <= 60)
})

test_that("ipca validates dimensions and lambda", {
  set.seed(2)
  X1 <- matrix(rnorm(20 * 8), 20, 8)
  X2 <- matrix(rnorm(21 * 7), 21, 7)
  expect_error(ipca(list(X1 = X1, X2 = X2), ncomp = 2, lambda = 1))

  X2_ok <- matrix(rnorm(20 * 7), 20, 7)
  expect_error(ipca(list(X1 = X1, X2 = X2_ok), ncomp = 2, lambda = 0))
  expect_error(ipca(list(X1 = X1, X2 = X2_ok), ncomp = 2, lambda = c(1, -1)))
  expect_error(ipca(list(X1 = X1, X2 = X2_ok), ncomp = 1.5, lambda = 1))
})

test_that("ipca auto method uses gram when p > n", {
  set.seed(3)
  blocks <- .sim_ipca_blocks(n = 20, r = 2, p_vec = c(40, 35, 30), noise_sd = 0.05)
  fit <- ipca(blocks, ncomp = 2, lambda = 1, method = "auto", max_iter = 40, tol = 1e-5)
  expect_identical(fit$method_used, "gram")
})

test_that("project.ipca supports matrix and list inputs", {
  set.seed(31)
  blocks <- .sim_ipca_blocks(n = 50, r = 2, p_vec = c(20, 25, 15), noise_sd = 0.04)
  fit <- ipca(blocks, ncomp = 2, lambda = 1, method = "gram", max_iter = 50, tol = 1e-6)

  pred_mat <- project(fit, do.call(cbind, blocks))
  pred_list <- project(fit, blocks)

  expect_true(is.matrix(pred_mat))
  expect_equal(dim(pred_mat), c(50, 2))
  expect_equal(unname(pred_mat), unname(pred_list), tolerance = 1e-8)

  S_fit <- multivarious::scores(fit)
  P1 <- pred_mat %*% solve(crossprod(pred_mat), t(pred_mat))
  P2 <- S_fit %*% solve(crossprod(S_fit), t(S_fit))
  rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
  expect_lt(rel, 5e-2)
})

test_that("ipca dense and gram methods agree on score subspace", {
  set.seed(4)
  blocks <- .sim_ipca_blocks(n = 50, r = 2, p_vec = c(12, 15, 10), noise_sd = 0.04)

  fit_dense <- ipca(blocks, ncomp = 2, lambda = 0.8, method = "dense", max_iter = 80, tol = 1e-7)
  fit_gram <- ipca(blocks, ncomp = 2, lambda = 0.8, method = "gram", max_iter = 80, tol = 1e-7)

  S1 <- multivarious::scores(fit_dense)
  S2 <- multivarious::scores(fit_gram)
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
  expect_lt(rel, 1e-3)
})

test_that("ipca handles high-dimensional blocks and returns finite outputs", {
  set.seed(5)
  blocks <- .sim_ipca_blocks(n = 30, r = 3, p_vec = c(120, 90, 110), noise_sd = 0.08)
  fit <- ipca(blocks, ncomp = 3, lambda = c(1, 0.8, 1.2), method = "gram", max_iter = 50, tol = 1e-5)

  expect_true(all(is.finite(multivarious::scores(fit))))
  expect_true(all(is.finite(fit$v)))
  expect_equal(ncol(multivarious::scores(fit)), 3)
})

test_that("ipca_tune_alpha returns valid selection and fit", {
  set.seed(6)
  blocks <- .sim_ipca_blocks(n = 40, r = 2, p_vec = c(25, 35, 30), noise_sd = 0.05)

  tuned <- ipca_tune_alpha(
    data = blocks,
    ncomp = 2,
    alpha_grid = c(0.1, 1, 10),
    tie = "inv_p",
    holdout_frac = 0.1,
    n_masks = 2,
    method = "gram",
    max_iter = 40,
    tol = 1e-5,
    seed = 123
  )

  expect_true(is.list(tuned))
  expect_true(tuned$best_alpha %in% c(0.1, 1, 10))
  expect_equal(length(tuned$best_lambda), 3)
  expect_equal(nrow(tuned$results), 3)
  expect_true(all(is.finite(tuned$results$holdout_mse)))
  expect_true(all(tuned$results$n_success >= 1))
  expect_true(all(tuned$results$converged_rate >= 0 & tuned$results$converged_rate <= 1))
  expect_s3_class(tuned$fit, "ipca")
})

test_that("ipca_tune_alpha tie schemes produce expected lambdas", {
  set.seed(7)
  blocks <- .sim_ipca_blocks(n = 30, r = 2, p_vec = c(10, 20, 40), noise_sd = 0.08)
  p_vec <- vapply(blocks, ncol, integer(1))

  tuned_invp <- ipca_tune_alpha(
    data = blocks,
    ncomp = 2,
    alpha_grid = c(2),
    tie = "inv_p",
    holdout_frac = 0.05,
    n_masks = 2,
    method = "gram",
    max_iter = 20,
    tol = 1e-4,
    seed = 42
  )
  expect_equal(unname(as.numeric(tuned_invp$best_lambda)), unname(2 / p_vec), tolerance = 1e-12)

  tuned_equal <- ipca_tune_alpha(
    data = blocks,
    ncomp = 2,
    alpha_grid = c(2),
    tie = "equal",
    holdout_frac = 0.05,
    n_masks = 2,
    method = "gram",
    max_iter = 20,
    tol = 1e-4,
    seed = 42
  )
  expect_equal(unname(as.numeric(tuned_equal$best_lambda)), unname(rep(2, length(p_vec))), tolerance = 1e-12)
})

test_that("ipca_tune_alpha warm start preserves alpha selection", {
  set.seed(8)
  blocks <- .sim_ipca_blocks(n = 40, r = 2, p_vec = c(25, 35, 30), noise_sd = 0.07)

  tuned_cold <- ipca_tune_alpha(
    data = blocks,
    ncomp = 2,
    alpha_grid = c(0.1, 1, 10),
    holdout_frac = 0.1,
    n_masks = 2,
    warm_start = FALSE,
    method = "gram",
    max_iter = 40,
    tol = 1e-5,
    seed = 11
  )

  tuned_warm <- ipca_tune_alpha(
    data = blocks,
    ncomp = 2,
    alpha_grid = c(0.1, 1, 10),
    holdout_frac = 0.1,
    n_masks = 2,
    warm_start = TRUE,
    method = "gram",
    max_iter = 40,
    tol = 1e-5,
    seed = 11
  )

  expect_true(all(is.finite(tuned_warm$results$holdout_mse)))
  expect_true(all(is.finite(tuned_cold$results$holdout_mse)))
  expect_lte(
    min(tuned_warm$results$holdout_mse),
    min(tuned_cold$results$holdout_mse) + 0.1
  )
  expect_true(isTRUE(tuned_warm$warm_start))
})

test_that("ipca truncated eigensolver matches dense full subspace", {
  set.seed(9)
  blocks <- .sim_ipca_blocks(n = 70, r = 3, p_vec = c(20, 24, 18), noise_sd = 0.05)

  fit_full <- ipca(
    blocks,
    ncomp = 3,
    lambda = 1,
    method = "dense",
    max_iter = 60,
    tol = 1e-6,
    eig_solver = "full"
  )

  fit_trunc <- ipca(
    blocks,
    ncomp = 3,
    lambda = 1,
    method = "dense",
    max_iter = 60,
    tol = 1e-6,
    eig_solver = "truncated",
    eig_rank = 10,
    eig_trunc_min_n = 2
  )

  expect_identical(fit_trunc$eig_solver_used, "truncated")
  S1 <- multivarious::scores(fit_full)
  S2 <- multivarious::scores(fit_trunc)
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
  expect_lt(rel, 1e-3)
})
