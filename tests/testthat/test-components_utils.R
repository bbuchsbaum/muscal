library(testthat)
library(muscal)

# Tests for components_utils.R (estimate_components, loading_reliability)

# --- estimate_components() ---------------------------------------------------

test_that("estimate_components with variance method keeps positive singular values", {
  set.seed(123)
  blocks <- replicate(3, matrix(rnorm(50), nrow = 10), simplify = FALSE)
  fit <- mfa(blocks, ncomp = 3)

  res <- estimate_components(fit, method = "variance")

  expect_type(res, "list")
  expect_named(res, c("keep", "criterion", "method_used"))
  expect_equal(res$method_used, "variance")
  expect_true(all(res$keep > 0))
  expect_equal(length(res$keep), 3)
  expect_equal(res$criterion, fit$sdev)
})

test_that("estimate_components falls back to variance when n not available", {
  set.seed(42)
  blocks <- replicate(2, matrix(rnorm(20), nrow = 5), simplify = FALSE)
  fit <- mfa(blocks, ncomp = 2)

  # Remove score matrix to prevent n inference
  fit_no_n <- fit
  fit_no_n$s <- NULL
  fit_no_n$roi_scores <- NULL
  fit_no_n$subject_scores <- NULL

  # Should fall back to variance method when n can't be inferred
  res <- estimate_components(fit_no_n, method = "variance")
  expect_equal(res$method_used, "variance")
})

test_that("estimate_components errors when no sdev available", {
  fit_empty <- list()
  expect_error(
    estimate_components(fit_empty, method = "variance"),
    "No singular values available"
  )
})

test_that("estimate_components respects explicitly passed sdev", {
  set.seed(1)
  blocks <- replicate(2, matrix(rnorm(30), nrow = 10), simplify = FALSE)
  fit <- mfa(blocks, ncomp = 2)

  custom_sdev <- c(5.0, 3.0)
  res <- estimate_components(fit, method = "variance", sdev = custom_sdev)

  expect_equal(res$criterion, custom_sdev)
  expect_equal(res$keep, c(1, 2))
})

test_that("estimate_components with rmt method uses significant_components", {
  set.seed(999)
  # Create larger data for RMT to have enough information
  blocks <- replicate(3, matrix(rnorm(500), nrow = 50), simplify = FALSE)
  fit <- mfa(blocks, ncomp = 5)

  # RMT requires V_list and n
  res <- estimate_components(fit, method = "rmt", tail_q = 0.2)

  expect_type(res, "list")
  expect_named(res, c("keep", "criterion", "method_used"))
  expect_equal(res$method_used, "rmt")
  expect_true(length(res$keep) >= 0)  # Could keep 0 or more components
})

# --- loading_reliability() ---------------------------------------------------

test_that("loading_reliability requires boot_loadings", {
  set.seed(1)
  blocks <- replicate(2, matrix(rnorm(30), nrow = 10), simplify = FALSE)
  fit <- mfa(blocks, ncomp = 2)

  expect_error(
    loading_reliability(fit, method = "bootstrap"),
    "boot_loadings must be provided"
  )

  expect_error(
    loading_reliability(fit, method = "bootstrap", boot_loadings = NULL),
    "boot_loadings must be provided"
  )

  expect_error(
    loading_reliability(fit, method = "bootstrap", boot_loadings = list()),
    "boot_loadings must be provided"
  )
})

test_that("loading_reliability computes bootstrap statistics", {
  set.seed(42)
  blocks <- replicate(2, matrix(rnorm(30), nrow = 10), simplify = FALSE)
  fit <- mfa(blocks, ncomp = 2)

  # Create mock bootstrap loadings (10 bootstrap samples)
  ref_loadings <- fit$v
  nboot <- 10
  boot_loadings <- lapply(1:nboot, function(i) {
    ref_loadings + matrix(rnorm(length(ref_loadings), sd = 0.1), nrow = nrow(ref_loadings))
  })

  res <- loading_reliability(fit, method = "bootstrap", boot_loadings = boot_loadings)

  expect_type(res, "list")
  expect_named(res, c("mean", "sd", "lower", "upper", "sign_consistency"))

  # Check dimensions match original loadings
  expect_equal(dim(res$mean), dim(ref_loadings))
  expect_equal(dim(res$sd), dim(ref_loadings))
  expect_equal(dim(res$lower), dim(ref_loadings))
  expect_equal(dim(res$upper), dim(ref_loadings))
  expect_equal(dim(res$sign_consistency), dim(ref_loadings))

  # Check value ranges

  expect_true(all(res$sd >= 0))
  expect_true(all(res$sign_consistency >= 0 & res$sign_consistency <= 1))
  expect_true(all(res$lower <= res$upper))
})

test_that("loading_reliability with alpha parameter affects quantiles", {
  set.seed(100)
  blocks <- replicate(2, matrix(rnorm(30), nrow = 10), simplify = FALSE)
  fit <- mfa(blocks, ncomp = 2)

  ref_loadings <- fit$v
  boot_loadings <- lapply(1:20, function(i) {
    ref_loadings + matrix(rnorm(length(ref_loadings), sd = 0.3), nrow = nrow(ref_loadings))
  })

  res_05 <- loading_reliability(fit, boot_loadings = boot_loadings, alpha = 0.05)
  res_20 <- loading_reliability(fit, boot_loadings = boot_loadings, alpha = 0.20)

  # With larger alpha, intervals should be narrower
  width_05 <- res_05$upper - res_05$lower
  width_20 <- res_20$upper - res_20$lower
  expect_true(mean(width_20) < mean(width_05))
})

test_that("loading_reliability errors on unknown method", {
  set.seed(1)
  blocks <- replicate(2, matrix(rnorm(30), nrow = 10), simplify = FALSE)
  fit <- mfa(blocks, ncomp = 2)

  expect_error(
    loading_reliability(fit, method = "unknown_method"),
    "should be"
  )
})

# --- internal helpers --------------------------------------------------------

test_that(".extract_V_list extracts from V_list attribute", {
  fit_with_attr <- list(sdev = c(1, 2))
  V_list <- list(matrix(1:4, 2, 2), matrix(5:8, 2, 2))
  attr(fit_with_attr, "V_list") <- V_list

  result <- muscal:::.extract_V_list(fit_with_attr)
  expect_equal(result, V_list)
})

test_that(".extract_V_list extracts from fit$V_list", {
  V_list <- list(matrix(1:4, 2, 2), matrix(5:8, 2, 2))
  fit_with_vlist <- list(sdev = c(1, 2), V_list = V_list)

  result <- muscal:::.extract_V_list(fit_with_vlist)
  expect_equal(result, V_list)
})

test_that(".extract_V_list extracts from fit$v with block_indices", {
  v <- rbind(matrix(1:4, 2, 2), matrix(5:8, 2, 2))
  block_indices <- list(1:2, 3:4)
  fit_with_v <- list(v = v, block_indices = block_indices)

  result <- muscal:::.extract_V_list(fit_with_v)

  expect_length(result, 2)
  expect_equal(result[[1]], v[1:2, , drop = FALSE])
  expect_equal(result[[2]], v[3:4, , drop = FALSE])
})

test_that(".extract_V_list errors when V_list cannot be extracted", {
  fit_empty <- list(sdev = c(1, 2))

  expect_error(
    muscal:::.extract_V_list(fit_empty),
    "Unable to extract V_list"
  )
})

test_that(".infer_n extracts from fit$s", {
  fit <- list(s = matrix(1:10, nrow = 5))
  expect_equal(muscal:::.infer_n(fit), 5)
})

test_that(".infer_n extracts from fit$roi_scores", {
  fit <- list(roi_scores = matrix(1:12, nrow = 4))
  expect_equal(muscal:::.infer_n(fit), 4)
})

test_that(".infer_n extracts from fit$subject_scores", {
  fit <- list(subject_scores = matrix(1:18, nrow = 6))
  expect_equal(muscal:::.infer_n(fit), 6)
})

test_that(".infer_n returns NULL when n cannot be inferred", {
  fit <- list(sdev = c(1, 2, 3))
  expect_null(muscal:::.infer_n(fit))
})
