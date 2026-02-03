library(testthat)
library(muscal)

# Tests for misc utility functions

# --- %||% operator -----------------------------------------------------------

test_that("%||% returns second arg when first is NULL", {
  expect_equal(NULL %||% 5, 5)
  expect_equal(NULL %||% "default", "default")
  expect_equal(NULL %||% list(a = 1), list(a = 1))
})

test_that("%||% returns first arg when not NULL", {
  expect_equal(3 %||% 4, 3)
  expect_equal("value" %||% "default", "value")
  expect_equal(list(a = 1) %||% list(b = 2), list(a = 1))
  expect_equal(0 %||% 5, 0)  # Zero is not NULL
  expect_equal(FALSE %||% TRUE, FALSE)  # FALSE is not NULL
})

# --- adam_update_block -------------------------------------------------------

test_that("adam_update_block computes correct momentum update", {
  set.seed(1)
  V <- matrix(rnorm(4), 2, 2)
  G <- matrix(rnorm(4), 2, 2)
  M <- matrix(0, 2, 2)
  V2 <- matrix(0, 2, 2)

  res <- muscal:::adam_update_block(V, G, M, V2,
    step_count = 1, beta1 = 0.9, beta2 = 0.999,
    adam_epsilon = 1e-8, learning_rate = 0.01
  )

  # Manual computation
  exp_M <- 0.1 * G
  exp_V2 <- 0.001 * (G * G)
  M_hat <- exp_M / (1 - 0.9)
  V_hat <- exp_V2 / (1 - 0.999)
  step <- 0.01 * M_hat / (sqrt(V_hat) + 1e-8)
  exp_V_new <- V - step

  expect_equal(res$M, exp_M)
  expect_equal(res$V2, exp_V2)
  expect_equal(res$V, exp_V_new)
})

test_that("adam_update_block handles multiple steps", {
  set.seed(42)
  V <- matrix(rnorm(4), 2, 2)
  G <- matrix(rnorm(4), 2, 2)
  M <- matrix(0.1, 2, 2)  # Non-zero starting momentum
  V2 <- matrix(0.01, 2, 2)

  res <- muscal:::adam_update_block(V, G, M, V2,
    step_count = 5, beta1 = 0.9, beta2 = 0.999,
    adam_epsilon = 1e-8, learning_rate = 0.001
  )

  expect_true(is.matrix(res$V))
  expect_true(is.matrix(res$M))
  expect_true(is.matrix(res$V2))
  expect_equal(dim(res$V), c(2, 2))
})

# --- within_class_scatter ----------------------------------------------------

test_that("within_class_scatter computes correct scatter matrix", {
  set.seed(123)
  # Simple example: 2 classes, 4 observations each
  X <- rbind(
    matrix(rnorm(8, mean = 0), 4, 2),
    matrix(rnorm(8, mean = 5), 4, 2)
  )
  classes <- factor(c(rep("A", 4), rep("B", 4)))

  Sw <- muscal:::within_class_scatter(X, classes)

  expect_true(is.matrix(Sw))
  expect_equal(dim(Sw), c(2, 2))
  expect_true(isSymmetric(Sw))
  # Sw should be positive semi-definite
  expect_true(all(eigen(Sw)$values >= -1e-10))
})

test_that("within_class_scatter handles single class", {
  set.seed(1)
  X <- matrix(rnorm(10), 5, 2)
  classes <- factor(rep("A", 5))

  Sw <- muscal:::within_class_scatter(X, classes)

  expect_equal(dim(Sw), c(2, 2))
  # For a single class, it's just the scatter around that class mean
  centered <- scale(X, center = TRUE, scale = FALSE)
  expected <- crossprod(centered)
  # Compare values, ignoring attributes (dimnames, etc.)
  expect_equal(as.numeric(Sw), as.numeric(expected))
})

# --- between_class_scatter ---------------------------------------------------

test_that("between_class_scatter computes correct scatter matrix", {
  set.seed(123)
  X <- rbind(
    matrix(rnorm(8, mean = 0), 4, 2),
    matrix(rnorm(8, mean = 5), 4, 2)
  )
  classes <- factor(c(rep("A", 4), rep("B", 4)))
  grand_mean <- colMeans(X)

  Sb <- muscal:::between_class_scatter(X, classes, grand_mean)

  expect_true(is.matrix(Sb))
  expect_equal(dim(Sb), c(2, 2))
  expect_true(isSymmetric(Sb))
  # Sb should be positive semi-definite
  expect_true(all(eigen(Sb)$values >= -1e-10))
})

test_that("between_class_scatter is zero for single class", {
  set.seed(1)
  X <- matrix(rnorm(10), 5, 2)
  classes <- factor(rep("A", 5))
  grand_mean <- colMeans(X)

  Sb <- muscal:::between_class_scatter(X, classes, grand_mean)

  # With one class, class mean = grand mean, so Sb should be 0
  expect_equal(Sb, matrix(0, 2, 2))
})

# --- significant_components --------------------------------------------------

test_that("significant_components validates inputs", {
  set.seed(1)
  # Create a fit with explicit V_list
  V_list <- list(matrix(rnorm(20), 5, 4), matrix(rnorm(20), 5, 4), matrix(rnorm(20), 5, 4))
  fit <- list(sdev = c(2, 1.5, 1, 0.5), V_list = V_list)

  # Invalid n
  expect_error(
    significant_components(fit, n = 0),
    "greater than"
  )

  # Mismatched k_vec length
  expect_error(
    significant_components(fit, n = 20, k_vec = c(5, 5)),
    "does not match"
  )
})

test_that("significant_components returns expected structure", {
  set.seed(42)
  # Create a fit object with explicit V_list (since mfa doesn't expose it directly)
  V_list <- list(
    matrix(rnorm(20), 5, 4),
    matrix(rnorm(20), 5, 4),
    matrix(rnorm(20), 5, 4)
  )
  fit <- list(
    sdev = c(3.5, 2.0, 1.0, 0.5),
    V_list = V_list
  )

  result <- significant_components(fit, n = 40, check_rmt = TRUE)

  expect_type(result, "list")
  expect_named(result, c(
    "keep", "rmt_pass", "icc_pass", "icc", "icc_pvalue",
    "mp_edge", "sigma2_est", "lambda", "n", "k_vec", "alpha"
  ))

  expect_type(result$keep, "integer")
  expect_type(result$rmt_pass, "logical")
  expect_type(result$icc_pass, "logical")
  expect_length(result$rmt_pass, 4)
  expect_length(result$icc_pass, 4)
  expect_equal(result$n, 40)
  expect_equal(result$alpha, 0.05)
})

test_that("significant_components skips RMT when sdev not available", {
  # Create fit without sdev
  fit_no_sdev <- list(
    V_list = list(matrix(rnorm(20), 5, 4), matrix(rnorm(20), 5, 4))
  )

  expect_warning(
    result <- significant_components(fit_no_sdev, n = 20, check_rmt = TRUE),
    "sdev"
  )

  # RMT should be skipped, all pass by default
  expect_true(all(result$rmt_pass))
})

test_that("significant_components with check_rmt = FALSE", {
  set.seed(1)
  # Create a fit object with explicit V_list
  V_list <- list(
    matrix(rnorm(15), 5, 3),
    matrix(rnorm(15), 5, 3)
  )
  fit <- list(
    sdev = c(2.0, 1.0, 0.5),
    V_list = V_list
  )

  expect_message(
    result <- significant_components(fit, n = 10, check_rmt = FALSE),
    "skipped"
  )

  # All should pass RMT when check disabled
  expect_true(all(result$rmt_pass))
})

test_that("significant_components errors when V_list missing", {
  fit_empty <- list(sdev = c(1, 2, 3))

  expect_error(
    significant_components(fit_empty, n = 10),
    "V_list"
  )
})

test_that("significant_components handles V_list as attribute", {
  set.seed(1)
  V_list <- list(matrix(rnorm(12), 4, 3), matrix(rnorm(12), 4, 3))
  fit <- list(sdev = c(2, 1, 0.5))
  attr(fit, "V_list") <- V_list

  result <- significant_components(fit, n = 20)

  expect_type(result, "list")
  expect_length(result$rmt_pass, 3)
})
