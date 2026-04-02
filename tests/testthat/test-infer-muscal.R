library(testthat)
library(muscal)

test_that("infer_muscal bootstraps component summaries across method families", {
  set.seed(31)
  X1 <- matrix(rnorm(48), nrow = 12)
  X2 <- matrix(rnorm(36), nrow = 12)
  Y <- matrix(rnorm(36), nrow = 12)

  mfa_fit <- mfa(list(X1 = X1, X2 = X2), ncomp = 2)
  anchored_fit <- anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = 1:12, X2 = 1:12),
    ncomp = 2
  )

  mfa_boot <- infer_muscal(mfa_fit, method = "bootstrap", statistic = "sdev", nrep = 6, seed = 101)
  anchored_boot <- infer_muscal(anchored_fit, method = "bootstrap", statistic = "sdev", nrep = 6, seed = 101)

  expect_s3_class(mfa_boot, "muscal_bootstrap_result")
  expect_s3_class(anchored_boot, "muscal_bootstrap_result")
  expect_true(all(c("component", "observed", "mean", "sd", "lower", "upper") %in% names(mfa_boot$summary)))
  expect_true(all(c("component", "observed", "mean", "sd", "lower", "upper") %in% names(anchored_boot$summary)))
  expect_equal(nrow(mfa_boot$summary), length(mfa_fit$sdev))
  expect_equal(nrow(anchored_boot$summary), length(anchored_fit$sdev))
  expect_equal(mfa_boot$n_success, 6)
  expect_equal(anchored_boot$n_success, 6)
})

test_that("infer_muscal can summarize bootstrap loading reliability", {
  set.seed(32)
  X1 <- matrix(rnorm(40), nrow = 10)
  X2 <- matrix(rnorm(30), nrow = 10)
  fit <- mfa(list(X1 = X1, X2 = X2), ncomp = 2)

  res <- infer_muscal(fit, method = "bootstrap", statistic = "loadings", nrep = 5, seed = 11)

  expect_s3_class(res, "muscal_bootstrap_result")
  expect_true(all(c("mean", "sd", "lower", "upper", "sign_consistency") %in% names(res$summary)))
  expect_equal(dim(res$summary$mean), dim(fit$v))
  expect_equal(dim(res$summary$sd), dim(fit$v))
})

test_that("infer_muscal permutation results are reproducible", {
  set.seed(33)
  X1 <- matrix(rnorm(48), nrow = 12)
  X2 <- matrix(rnorm(36), nrow = 12)
  fit <- mfa(list(X1 = X1, X2 = X2), ncomp = 2)

  res1 <- infer_muscal(fit, method = "permutation", statistic = "sdev", nrep = 7, seed = 19)
  res2 <- infer_muscal(fit, method = "permutation", statistic = "sdev", nrep = 7, seed = 19)

  expect_s3_class(res1, "muscal_permutation_result")
  expect_equal(res1$perm_values, res2$perm_values)
  expect_equal(res1$component_results$p_value, res2$component_results$p_value)
  expect_equal(res1$component_results$observed, as.numeric(fit$sdev))
})

test_that("infer_muscal permutation distinguishes planted shared signal from null", {
  set.seed(37)
  n <- 40
  k <- 2

  Z <- matrix(rnorm(n * k), nrow = n, ncol = k)
  A1 <- matrix(rnorm(12 * k), nrow = 12, ncol = k)
  A2 <- matrix(rnorm(10 * k), nrow = 10, ncol = k)
  X1_signal <- Z %*% t(A1) + matrix(rnorm(n * 12, sd = 0.05), nrow = n, ncol = 12)
  X2_signal <- Z %*% t(A2) + matrix(rnorm(n * 10, sd = 0.05), nrow = n, ncol = 10)

  signal_fit <- mfa(list(X1 = X1_signal, X2 = X2_signal), ncomp = 2)
  signal_perm <- infer_muscal(
    signal_fit,
    method = "permutation",
    statistic = "sdev",
    nrep = 29,
    seed = 11
  )

  X1_null <- matrix(rnorm(n * 12), nrow = n, ncol = 12)
  X2_null <- matrix(rnorm(n * 10), nrow = n, ncol = 10)
  null_fit <- mfa(list(X1 = X1_null, X2 = X2_null), ncomp = 2)
  null_perm <- infer_muscal(
    null_fit,
    method = "permutation",
    statistic = "sdev",
    nrep = 29,
    seed = 11
  )

  expect_lt(signal_perm$component_results$p_value[[1]], 0.1)
  expect_gt(signal_perm$component_results$observed[[1]], signal_perm$component_results$upper_ci[[1]])
  expect_gt(null_perm$component_results$p_value[[1]], 0.15)
  expect_gt(
    null_perm$component_results$p_value[[1]],
    signal_perm$component_results$p_value[[1]] + 0.1
  )
})

test_that("infer_muscal errors for fits without refit metadata", {
  set.seed(34)
  X1 <- matrix(rnorm(48), nrow = 12)
  X2 <- matrix(rnorm(36), nrow = 12)
  fit <- aligned_mfa(
    list(X1 = X1, X2 = X2),
    row_index = list(1:12, 1:12),
    ncomp = 2
  )

  expect_error(
    infer_muscal(fit, method = "bootstrap", nrep = 3),
    "does not expose refit metadata"
  )
})

test_that("infer_muscal honors explicit refit overrides and records partial failures", {
  set.seed(35)
  X1 <- matrix(rnorm(48), nrow = 12)
  X2 <- matrix(rnorm(36), nrow = 12)
  blocks <- list(X1 = X1, X2 = X2)
  fit <- mfa(blocks, ncomp = 2)

  bootstrap_counter <- local({
    i <- 0L
    function(data) {
      i <<- i + 1L
      if (i %% 2L == 0L) {
        list(blocks = data$blocks, fail = TRUE)
      } else {
        idx <- sample.int(nrow(data$blocks[[1L]]), replace = TRUE)
        list(
          blocks = lapply(data$blocks, function(block) block[idx, , drop = FALSE]),
          fail = FALSE
        )
      }
    }
  })

  refit <- list(
    data = list(blocks = blocks),
    fit_fn = function(data) {
      if (isTRUE(data$fail)) {
        stop("intentional refit failure", call. = FALSE)
      }
      mfa(data$blocks, ncomp = 2)
    },
    bootstrap_fn = bootstrap_counter,
    permutation_fn = function(data) data,
    resample_unit = "rows"
  )

  res <- infer_muscal(
    fit,
    method = "bootstrap",
    statistic = "sdev",
    refit = refit,
    nrep = 4,
    seed = 101
  )

  expect_s3_class(res, "muscal_bootstrap_result")
  expect_equal(res$n_success, 2)
  expect_length(res$failures, 2)
  expect_true(all(vapply(res$failures, `[[`, character(1), "stage") == "fit"))
})

test_that("infer_muscal rejects matrix-valued permutation statistics", {
  set.seed(36)
  X1 <- matrix(rnorm(48), nrow = 12)
  X2 <- matrix(rnorm(36), nrow = 12)
  fit <- mfa(list(X1 = X1, X2 = X2), ncomp = 2)

  expect_error(
    infer_muscal(fit, method = "permutation", statistic = "loadings", nrep = 3),
    "numeric vector statistic"
  )
})
