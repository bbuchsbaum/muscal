library(testthat)
library(musca)

skip_if_not_installed("abind")

# Tests for summarize_matrix

test_that("summarize_matrix computes element-wise function", {
  mats <- list(matrix(1:4, nrow = 2), matrix(2:5, nrow = 2))
  res_mean <- musca:::summarize_matrix(mats, mean)
  expected_mean <- (mats[[1]] + mats[[2]]) / 2
  expect_equal(res_mean, expected_mean)
})

# Tests for summarize_boot

test_that("summarize_boot aggregates bootstrap results", {
  boot_i <- list(matrix(1:4, nrow = 2), matrix(2:5, nrow = 2))
  boot_j <- list(matrix(10:13, nrow = 2), matrix(11:14, nrow = 2))
  boot_ret <- list(boot_i = boot_i, boot_j = boot_j)
  res <- musca:::summarize_boot(boot_ret, alpha = 0.25)

  mean_i <- (boot_i[[1]] + boot_i[[2]]) / 2
  sd_val <- matrix(sqrt(0.5), 2, 2)
  upper_i <- mean_i + 0.5
  lower_i <- mean_i - 0.5

  mean_j <- (boot_j[[1]] + boot_j[[2]]) / 2
  upper_j <- mean_j + 0.5
  lower_j <- mean_j - 0.5

  expect_equal(res$boot_scores_mean, mean_i)
  expect_equal(res$boot_scores_sd, sd_val)
  expect_equal(res$boot_scores_upper, upper_i)
  expect_equal(res$boot_scores_lower, lower_i)
  expect_equal(res$boot_lds_mean, mean_j)
  expect_equal(res$boot_lds_upper, upper_j)
  expect_equal(res$boot_lds_lower, lower_j)
})

# Tests for significant_components with RMT

test_that("significant_components performs RMT and ICC", {
  V1 <- matrix(c(1,0,0,0,1,0), 3, 2)
  V2 <- matrix(c(0,0,1,0,0,1), 3, 2)
  fit <- list(V_list = list(V1, V2), sdev = c(1, 0.5))
  res <- significant_components(fit, n = 10, k_vec = c(3, 3), check_rmt = TRUE)
  expect_equal(res$rmt_pass, c(TRUE, FALSE))
  expect_equal(res$keep, integer(0))
  expect_equal(res$icc, c(0, 0))
  expect_equal(res$mp_edge, 0.625, tolerance = 1e-6)
})

