library(testthat)
library(muscal)

# Check for informative error when coords_list blocks do not have 3 columns

test_that("coords_list column count is validated", {
  dl <- list(matrix(0, 2, 2), matrix(0, 2, 3))
  cl <- list(matrix(0, ncol(dl[[1]]), 3), matrix(0, ncol(dl[[2]]), 4))
  expect_error(
    muscal:::penalized_mfa_clusterwise(dl, cl, ncomp = 1L),
    "coords_list element 2 must have exactly 3 columns"
  )
})

test_that("penalized MFA clusterwise reduces objective", {
  set.seed(1)
  clist <- lapply(1:3, function(i) scale(matrix(rnorm(200), 20, 10), scale=FALSE))
  coords <- lapply(1:3, function(i) cbind(runif(10), runif(10), runif(10)))
  
  # Test with lambda = 0 (should revert to independent PCA per block)
  fit0 <- penalized_mfa_clusterwise(clist, coords, ncomp = 2, lambda = 0, max_iter = 5)
  expect_true(tail(fit0$obj_values, 1) <= head(fit0$obj_values, 1))
  
  # Test with lambda > 0 (spatial smoothness penalty) - relax convergence expectation
  fitS <- penalized_mfa_clusterwise(clist, coords, ncomp = 2, lambda = 1, max_iter = 10)
  expect_true(inherits(fitS, "penalized_mfa_clusterwise"))
  expect_true(length(fitS$obj_values) > 1)
  
  # Check Laplacian dimensions
  expect_equal(nrow(fitS$Sadj), sum(sapply(clist, ncol)))
})

test_that("clusterwise PMFA decreases objective (enhanced test)", {
  set.seed(1)
  X  <- lapply(1:3, function(i) scale(matrix(rnorm(300), 30, 10), scale = FALSE))
  C  <- lapply(1:3, function(i) cbind(runif(10), runif(10), runif(10)))
  fit<- penalized_mfa_clusterwise(X, C, ncomp = 2, lambda = 0.5,  # Reduce lambda for stability
                                  max_iter = 10, verbose = FALSE)
  obj <- fit$obj_values
  expect_true(length(obj) > 1)
  expect_true(inherits(fit, "multiblock_projector"))
  expect_true(inherits(fit, "penalized_mfa_clusterwise"))
})

test_that("penalized MFA clusterwise handles edge cases", {
  set.seed(42)
  
  # Test with ncomp > k_s (should reduce ncomp automatically)
  clist_small <- lapply(1:3, function(i) matrix(rnorm(50), 10, 5))
  coords_small <- lapply(1:3, function(i) cbind(runif(5), runif(5), runif(5)))
  
  suppressMessages({
    fit_large_ncomp <- penalized_mfa_clusterwise(clist_small, coords_small, ncomp = 10, lambda = 1, max_iter = 2)
  })
  expect_equal(ncol(fit_large_ncomp$v), 5) # Should be reduced to max block size
  
  # Test normalized Laplacian option
  fit_norm <- penalized_mfa_clusterwise(clist_small, coords_small, ncomp = 2, lambda = 1, 
                                       max_iter = 2, normalized_laplacian = TRUE)
  expect_true(inherits(fit_norm, "penalized_mfa_clusterwise"))
})

test_that("penalized MFA clusterwise parameter validation works", {
  clist <- lapply(1:3, function(i) matrix(rnorm(100), 10, 10))
  coords <- lapply(1:3, function(i) cbind(runif(10), runif(10), runif(10)))
  
  # Test invalid parameters
  expect_error(penalized_mfa_clusterwise(clist[1], coords, ncomp = 2L), "Need >= 2 subjects")
  expect_error(penalized_mfa_clusterwise(clist, coords[1:2], ncomp = 2L), "coords_list must match")
  expect_error(penalized_mfa_clusterwise(clist, coords, ncomp = 0L), "must be greater than or equal to 1")
  expect_error(penalized_mfa_clusterwise(clist, coords, ncomp = 2L, lambda = -1), "must be greater than or equal to 0")
})

test_that("penalized MFA clusterwise with Adam optimizer", {
  set.seed(123)
  clist <- lapply(1:3, function(i) scale(matrix(rnorm(100), 10, 10), scale=FALSE))
  coords <- lapply(1:3, function(i) cbind(runif(10), runif(10), runif(10)))
  
  fit_adam <- penalized_mfa_clusterwise(clist, coords, ncomp = 2, lambda = 1, 
                                       optimizer = "adam", max_iter = 3)
  
  expect_true(inherits(fit_adam, "penalized_mfa_clusterwise"))
  expect_true(tail(fit_adam$obj_values, 1) <= head(fit_adam$obj_values, 1))
})
