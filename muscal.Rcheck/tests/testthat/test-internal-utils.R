library(testthat)
library(muscal)

set.seed(1)
make_cov <- function(seed) {
  set.seed(seed)
  A <- matrix(rnorm(9), 3, 3)
  tcrossprod(A)
}
subject_data <- lapply(1:3, make_cov)
res_none <- covstatis(subject_data, ncomp = 2, norm_method = "none", dcenter = FALSE)
res_mfa  <- covstatis(subject_data, ncomp = 2, norm_method = "mfa", dcenter = FALSE)

# --- .pre_process_new_cov ----------------------------------------------------

test_that(".pre_process_new_cov validates symmetry", {
  bad <- matrix(1:9, 3, 3)
  expect_error(muscal:::.pre_process_new_cov(res_none, bad), "symmetric")
})

test_that(".pre_process_new_cov validates dimensions", {
  bad <- matrix(0, 4, 4)
  expect_error(muscal:::.pre_process_new_cov(res_none, bad), "same dimensions")
})

test_that(".pre_process_new_cov checks eigenvalue for MFA normalization", {
  zero_mat <- matrix(0, 3, 3)
  expect_error(muscal:::.pre_process_new_cov(res_mfa, zero_mat), "too small")
})

# --- ls_ridge ---------------------------------------------------------------

test_that("ls_ridge handles zero predictors and ridge solution", {
  B0 <- muscal:::ls_ridge(matrix(nrow = 5, ncol = 0), matrix(1:10, 5, 2))
  expect_equal(dim(B0), c(0, 2))

  set.seed(2)
  Z <- matrix(rnorm(20), 5, 4)
  X <- matrix(rnorm(10), 5, 2)
  lambda <- 0.5
  coef_est <- muscal:::ls_ridge(Z, X, lambda = lambda)
  coef_manual <- solve(t(Z) %*% Z + lambda * diag(4), t(Z) %*% X)
  expect_equal(coef_est, coef_manual, tolerance = 1e-8)
})

# --- project_covariate ------------------------------------------------------

test_that("project_covariate beta scaling matches manual formula", {
  y <- rnorm(length(res_none$partial_scores))
  beta_res <- project_covariate(res_none, y, what = "dimension", scale = "beta")
  G <- muscal:::.get_G_scores(res_none)
  manual <- as.numeric(t(G) %*% y / colSums(G^2))
  expect_equal(as.numeric(beta_res), manual)
})

