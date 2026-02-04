library(testthat)
library(muscal)

# Test data setup
set.seed(123)
X1 <- scale(matrix(rnorm(100 * 4), 100, 4), scale = FALSE)
X2 <- scale(matrix(rnorm(100 * 16), 100, 16), scale = FALSE)
data_list_uneven <- list(B1 = X1, B2 = X2)
data_list_even <- list(B1 = X1, B2 = X1)

test_that("penalized_mfa.list handles uneven and even blocks correctly", {
  # Test with uneven blocks: should warn and produce no consensus matrix
  expect_warning(
    res_uneven <- penalized_mfa.list(data_list_uneven, ncomp = 2, lambda = 1, compute_consensus = TRUE, verbose = FALSE, penalty_method = "projection"),
    "Cannot compute consensus"
  )
  expect_s3_class(res_uneven, "penalized_mfa")
  expect_null(attr(res_uneven, "consensus"))

  # Test with even blocks: should produce a consensus matrix
  res_even <- penalized_mfa.list(data_list_even, ncomp = 2, lambda = 1, compute_consensus = TRUE, verbose = FALSE, penalty_method = "projection")
  expect_s3_class(res_even, "penalized_mfa")
  expect_true(!is.null(attr(res_even, "consensus")))
  expect_equal(dim(attr(res_even, "consensus")), c(4, 2))
})

test_that("penalized_mfa works with multiblock input", {
  mb <- multiblock(data_list_uneven)
  # Should warn about disabling penalty when blocks have different dimensions
  expect_warning(
    res <- penalized_mfa(mb, ncomp = 2, lambda = 1, penalty_method = "projection"),
    "Cannot use 'projection' penalty with blocks of different dimensions"
  )
  expect_s3_class(res, "penalized_mfa")
})

test_that("penalized_mfa stops with fewer than 2 blocks", {
  expect_error(penalized_mfa.list(list(X1), ncomp = 2), "requires at least 2 data blocks")
})
