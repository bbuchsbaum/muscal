library(testthat)
library(muscal)
library(multidesign)

skip_if_not_installed("multidesign")

# --- bada.multidesign basic tests --------------------------------------------

test_that("bada.multidesign works with basic splitting", {
  set.seed(123)
  x <- matrix(rnorm(20), nrow = 10)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = 5),
    y = rep(c("A", "B"), 5)
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 1)

  expect_s3_class(res, "bada")
  expect_equal(levels(res$subjects), c("S1", "S2"))
  expect_equal(levels(res$labels), c("A", "B"))
  expect_length(res$block_indices, 2)
  expect_true(multivarious::ncomp(res) >= 1)
})

test_that("bada.multidesign produces valid projector", {
  set.seed(42)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 5), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 2)

  expect_s3_class(res, "bada")
  expect_true(inherits(res, "projector"))
  expect_true(!is.null(res$v))
  expect_true(!is.null(res$s))
  expect_true(!is.null(res$sdev))
  expect_equal(nrow(res$s), n_obs)
})

test_that("bada.multidesign with multiple classes", {
  set.seed(1)
  n_obs <- 30
  x <- matrix(rnorm(n_obs * 6), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2", "S3"), each = n_obs / 3),
    y = factor(rep(c("A", "B", "C"), n_obs / 3))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 2)

  expect_s3_class(res, "bada")
  expect_equal(length(levels(res$labels)), 3)
  expect_equal(length(levels(res$subjects)), 3)
  # Should have barycenters for each class
  expect_equal(nrow(res$barycenters), 3)
})

# --- bada.multidesign with residual analysis ---------------------------------

test_that("bada.multidesign with residual components", {
  set.seed(99)
  n_obs <- 40
  x <- matrix(rnorm(n_obs * 10), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  # Using rescomp > 0 enables residual analysis
  res <- bada(md, y = y, subject = subj_id, ncomp = 1, resdim = 5, rescomp = 2)

  expect_s3_class(res, "bada")
  expect_equal(res$resdim, 5)
  expect_equal(res$rescomp, 2)
  # Total components = ncomp (class) + rescomp (residual)
  expect_equal(ncol(res$v), 3)  # 1 + 2
})

# --- reprocess.bada tests ----------------------------------------------------

test_that("reprocess.bada with no block applies averaged preprocessing", {
  set.seed(123)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 8), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 1)

  # Create new data with same dimensions
  new_data <- matrix(rnorm(5 * 8), nrow = 5)
  reprocessed <- reprocess(res, new_data)

  expect_true(is.matrix(reprocessed))
  expect_equal(nrow(reprocessed), 5)
  expect_equal(ncol(reprocessed), ncol(x))
})

test_that("reprocess.bada with block uses block-specific preprocessing", {
  set.seed(456)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 8), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 1)

  new_data <- matrix(rnorm(5 * 8), nrow = 5)
  reprocessed <- reprocess(res, new_data, block = "S1")

  expect_true(is.matrix(reprocessed))
  expect_equal(nrow(reprocessed), 5)
})

test_that("reprocess.bada validates block name", {
  set.seed(789)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 8), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 1)
  new_data <- matrix(rnorm(5 * 8), nrow = 5)

  expect_error(
    reprocess(res, new_data, block = "INVALID"),
    class = "error"
  )
})

test_that("reprocess.bada validates data dimensions", {
  set.seed(111)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 8), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 1)

  # Wrong number of columns should error
  new_data_wrong <- matrix(rnorm(5 * 3), nrow = 5)  # 3 cols instead of 8
  expect_error(reprocess(res, new_data_wrong))
})

# --- project.bada tests ------------------------------------------------------

test_that("project.bada projects new data", {
  set.seed(222)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 8), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 2)

  new_data <- matrix(rnorm(5 * 8), nrow = 5)
  projected <- project(res, new_data)

  expect_true(is.matrix(projected))
  expect_equal(nrow(projected), 5)
  # ncomp may be capped by number of classes (2) minus 1
  expect_true(ncol(projected) >= 1)
})

test_that("project.bada with block parameter", {
  set.seed(333)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 8), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 2)

  new_data <- matrix(rnorm(5 * 8), nrow = 5)
  projected <- project(res, new_data, block = "S2")

  expect_true(is.matrix(projected))
  expect_equal(nrow(projected), 5)
  # ncomp may be capped by number of classes
  expect_true(ncol(projected) >= 1)
})

# --- summarize_boot internal function ----------------------------------------

test_that("summarize_boot computes correct statistics", {
  set.seed(100)

  # Create mock bootstrap result structure
  boot_ret <- data.frame(i = 1:5)
  boot_ret$boot_i <- lapply(1:5, function(x) matrix(rnorm(6), 2, 3))
  boot_ret$boot_j <- lapply(1:5, function(x) matrix(rnorm(6), 2, 3))

  result <- muscal:::summarize_boot(boot_ret, alpha = 0.05)

  expect_type(result, "list")
  expect_named(result, c(
    "boot_scores_mean", "boot_scores_sd", "boot_scores_upper", "boot_scores_lower",
    "boot_lds_mean", "boot_lds_sd", "boot_lds_upper", "boot_lds_lower"
  ))

  # Check dimensions
  expect_equal(dim(result$boot_scores_mean), c(2, 3))
  expect_equal(dim(result$boot_lds_mean), c(2, 3))

  # Check that sd is non-negative
  expect_true(all(result$boot_scores_sd >= 0))
  expect_true(all(result$boot_lds_sd >= 0))
})

# --- summarize_matrix internal function --------------------------------------

test_that("summarize_matrix applies function across matrices", {
  matrices <- list(
    matrix(1:4, 2, 2),
    matrix(5:8, 2, 2),
    matrix(9:12, 2, 2)
  )

  result <- muscal:::summarize_matrix(matrices, mean)

  expected <- matrix(c(5, 6, 7, 8), 2, 2)
  expect_equal(result, expected)
})

test_that("summarize_matrix with sd function", {
  set.seed(42)
  matrices <- lapply(1:10, function(x) matrix(rnorm(4), 2, 2))

  result <- muscal:::summarize_matrix(matrices, sd)

  expect_equal(dim(result), c(2, 2))
  expect_true(all(result >= 0))
})

# --- bootstrap.bada test (simplified) ----------------------------------------

test_that("bootstrap.bada requires data argument", {
  set.seed(444)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 6), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 1)

  expect_error(
    multivarious::bootstrap(res, nboot = 10),
    "data.*required"
  )
})
