library(testthat)
library(musca)

# Basic MFA run should produce a valid object

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
  nf <- musca:::normalization_factors(list(X1, X2), type = "MFA")
  expect_equal(nf[1], 1/(svd(X1)$d[1]^2))
  expect_equal(nf[2], 1/(svd(X2)$d[1]^2))
})
