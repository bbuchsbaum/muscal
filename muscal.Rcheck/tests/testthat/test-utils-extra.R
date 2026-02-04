library(testthat)
library(muscal)

# Test compute_sim_mat with a simple inner product similarity

test_that("compute_sim_mat returns correct symmetric matrix", {
  blocks <- list(
    matrix(1:4, nrow = 2),
    matrix(2:5, nrow = 2),
    matrix(3:6, nrow = 2)
  )
  simfun <- function(a, b) sum(a * b)
  M <- muscal:::compute_sim_mat(blocks, simfun)
  expected <- matrix(0, 3, 3)
  expected[1,2] <- expected[2,1] <- sum(blocks[[1]] * blocks[[2]])
  expected[1,3] <- expected[3,1] <- sum(blocks[[1]] * blocks[[3]])
  expected[2,3] <- expected[3,2] <- sum(blocks[[2]] * blocks[[3]])
  expect_equal(M, expected)
})

# Test normalization_factors for Frob and None cases

test_that("normalization_factors handles Frob and None", {
  blocks <- list(matrix(1:4, nrow = 2), matrix(2:5, nrow = 2))
  nf_frob <- muscal:::normalization_factors(blocks, type = "Frob")
  nf_none <- muscal:::normalization_factors(blocks, type = "None")
  expect_equal(nf_frob, c(sum(blocks[[1]]^2), sum(blocks[[2]]^2)))
  expect_equal(nf_none, c(1, 1))
})

# Test structure of significant_components when RMT check is disabled

test_that("significant_components basic structure", {
  V_list <- list(diag(2), diag(2))
  fit <- list(V_list = V_list, sdev = c(1, 0.5))
  res <- muscal:::significant_components(fit, n = 10, k_vec = c(2, 2), check_rmt = FALSE)
  expect_type(res, "list")
  expect_true(all(c("keep", "rmt_pass", "icc_pass") %in% names(res)))
  expect_length(res$rmt_pass, 2)
  expect_length(res$icc_pass, 2)
})

# Diagnostics for penalized_mfa_clusterwise argument checks

test_that("coords_list length mismatch triggers error", {
  dl <- list(matrix(0, 2, 2), matrix(0, 2, 2))
  cl <- list(matrix(0, 2, 3))
  expect_error(penalized_mfa_clusterwise(dl, cl, ncomp = 1L), "coords_list")
})

test_that("coords_list row count mismatch triggers error", {
  dl <- list(matrix(0, 2, 2), matrix(0, 3, 2))
  cl <- list(matrix(0, 2, 3), matrix(0, 4, 3))
  expect_error(
    penalized_mfa_clusterwise(dl, cl, ncomp = 1L),
    "rows in each coords_list element"
  )
})
