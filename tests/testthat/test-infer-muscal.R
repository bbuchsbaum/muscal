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
  aligned_fit <- aligned_mfa(
    list(X1 = X1, X2 = X2),
    row_index = list(X1 = 1:12, X2 = 1:12),
    N = 12,
    ncomp = 2
  )
  mcca_fit <- mcca(
    list(X1 = X1, X2 = X2),
    ncomp = 2,
    ridge = 1e-6
  )
  ipca_fit <- ipca(
    list(X1 = X1, X2 = X2),
    ncomp = 2,
    lambda = 1,
    method = "gram",
    max_iter = 40,
    tol = 1e-5
  )
  graph_fit <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = 1:12, X2 = 1:12),
    ncomp = 2,
    graph_lambda = 0
  )
  Y_rrr <- list(
    X1 = X1[, 1:3, drop = FALSE] + matrix(rnorm(12 * 3, sd = 0.05), 12, 3),
    X2 = X2[, 1:3, drop = FALSE] + matrix(rnorm(12 * 3, sd = 0.05), 12, 3)
  )
  rrr_fit <- aligned_rrr(
    Y = Y_rrr,
    X = list(X1 = X1, X2 = X2),
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    ridge = 1e-8,
    max_iter = 60,
    tol = 1e-8
  )

  mfa_boot <- infer_muscal(mfa_fit, method = "bootstrap", statistic = "sdev", nrep = 6, seed = 101)
  anchored_boot <- infer_muscal(anchored_fit, method = "bootstrap", statistic = "sdev", nrep = 6, seed = 101)
  aligned_boot <- infer_muscal(aligned_fit, method = "bootstrap", statistic = "sdev", nrep = 6, seed = 101)
  mcca_boot <- infer_muscal(mcca_fit, method = "bootstrap", statistic = "sdev", nrep = 6, seed = 101)
  ipca_boot <- infer_muscal(ipca_fit, method = "bootstrap", statistic = "sdev", nrep = 6, seed = 101)
  graph_boot <- infer_muscal(graph_fit, method = "bootstrap", statistic = "sdev", nrep = 6, seed = 101)
  rrr_boot <- infer_muscal(rrr_fit, method = "bootstrap", statistic = "sdev", nrep = 6, seed = 101)

  expect_s3_class(mfa_boot, "muscal_bootstrap_result")
  expect_s3_class(anchored_boot, "muscal_bootstrap_result")
  expect_s3_class(aligned_boot, "muscal_bootstrap_result")
  expect_s3_class(mcca_boot, "muscal_bootstrap_result")
  expect_s3_class(ipca_boot, "muscal_bootstrap_result")
  expect_s3_class(graph_boot, "muscal_bootstrap_result")
  expect_s3_class(rrr_boot, "muscal_bootstrap_result")
  expect_true(all(c("component", "observed", "mean", "sd", "lower", "upper") %in% names(mfa_boot$summary)))
  expect_true(all(c("component", "observed", "mean", "sd", "lower", "upper") %in% names(anchored_boot$summary)))
  expect_true(all(c("component", "observed", "mean", "sd", "lower", "upper") %in% names(aligned_boot$summary)))
  expect_true(all(c("component", "observed", "mean", "sd", "lower", "upper") %in% names(mcca_boot$summary)))
  expect_true(all(c("component", "observed", "mean", "sd", "lower", "upper") %in% names(ipca_boot$summary)))
  expect_true(all(c("component", "observed", "mean", "sd", "lower", "upper") %in% names(graph_boot$summary)))
  expect_true(all(c("component", "observed", "mean", "sd", "lower", "upper") %in% names(rrr_boot$summary)))
  expect_equal(nrow(mfa_boot$summary), length(mfa_fit$sdev))
  expect_equal(nrow(anchored_boot$summary), length(anchored_fit$sdev))
  expect_equal(nrow(aligned_boot$summary), length(aligned_fit$sdev))
  expect_equal(nrow(mcca_boot$summary), length(mcca_fit$sdev))
  expect_equal(nrow(ipca_boot$summary), length(ipca_fit$sdev))
  expect_equal(nrow(graph_boot$summary), length(graph_fit$sdev))
  expect_equal(nrow(rrr_boot$summary), length(rrr_fit$sdev))
  expect_equal(mfa_boot$n_success, 6)
  expect_equal(anchored_boot$n_success, 6)
  expect_equal(aligned_boot$n_success, 6)
  expect_equal(mcca_boot$n_success, 6)
  expect_equal(ipca_boot$n_success, 6)
  expect_equal(graph_boot$n_success, 6)
  expect_equal(rrr_boot$n_success, 6)
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

test_that("infer_muscal supports aligned_mfa bootstrap and permutation summaries", {
  set.seed(34)
  X1 <- matrix(rnorm(48), nrow = 12)
  X2 <- matrix(rnorm(36), nrow = 12)
  fit <- aligned_mfa(
    list(X1 = X1, X2 = X2),
    row_index = list(1:12, 1:12),
    N = 12,
    ncomp = 2
  )

  boot <- infer_muscal(fit, method = "bootstrap", statistic = "sdev", nrep = 4, seed = 13)
  perm <- infer_muscal(fit, method = "permutation", statistic = "sdev", nrep = 5, seed = 13)

  expect_s3_class(boot, "muscal_bootstrap_result")
  expect_s3_class(perm, "muscal_permutation_result")
  expect_equal(boot$n_success, 4)
  expect_equal(perm$n_success, 5)
  expect_true(all(is.finite(boot$summary$mean)))
  expect_true(all(perm$component_results$p_value >= 0 & perm$component_results$p_value <= 1))
})

test_that("infer_muscal supports graph_anchored_mfa bootstrap and permutation summaries", {
  set.seed(39)
  Y <- matrix(rnorm(48), nrow = 12)
  X1 <- matrix(rnorm(60), nrow = 12)
  X2 <- matrix(rnorm(60), nrow = 12)
  colnames(X1) <- paste0("f", 1:ncol(X1))
  colnames(X2) <- c(paste0("f", 1:4), "u1")

  fit <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = 1:12, X2 = 1:12),
    ncomp = 2,
    feature_graph = "colnames",
    graph_lambda = 2
  )

  boot <- infer_muscal(fit, method = "bootstrap", statistic = "sdev", nrep = 4, seed = 21)
  perm <- infer_muscal(fit, method = "permutation", statistic = "sdev", nrep = 5, seed = 21)

  expect_s3_class(boot, "muscal_bootstrap_result")
  expect_s3_class(perm, "muscal_permutation_result")
  expect_equal(boot$n_success, 4)
  expect_equal(perm$n_success, 5)
  expect_true(all(is.finite(boot$summary$mean)))
  expect_true(all(perm$component_results$p_value >= 0 & perm$component_results$p_value <= 1))
})

test_that("infer_muscal supports response_aligned_mfa with scalar row-weight inputs", {
  set.seed(42)
  sim <- list(
    X1 = matrix(rnorm(26 * 6), nrow = 26),
    X2 = matrix(rnorm(24 * 5), nrow = 24)
  )
  Z1 <- scale(matrix(rnorm(26 * 2), 26, 2), center = TRUE, scale = FALSE)
  Z2 <- scale(matrix(rnorm(24 * 2), 24, 2), center = TRUE, scale = FALSE)
  B <- matrix(rnorm(3 * 2), 3, 2)
  Y <- list(
    X1 = Z1 %*% t(B) + matrix(rnorm(26 * 3, sd = 0.05), 26, 3),
    X2 = Z2 %*% t(B) + matrix(rnorm(24 * 3, sd = 0.05), 24, 3)
  )
  idx <- list(
    X1 = sample.int(9, nrow(sim$X1), replace = TRUE),
    X2 = sample.int(9, nrow(sim$X2), replace = TRUE)
  )

  fit <- response_aligned_mfa(
    Y = Y,
    X = sim,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    response_weights = c(X1 = 0.7, X2 = 1.3),
    anchor_map = idx,
    anchor_weight = 0.8,
    coupling_lambda = 2,
    max_iter = 60,
    tol = 1e-8,
    ridge = 1e-8
  )

  boot <- infer_muscal(fit, method = "bootstrap", statistic = "sdev", nrep = 4, seed = 17)
  perm <- infer_muscal(fit, method = "permutation", statistic = "sdev", nrep = 4, seed = 17)

  expect_s3_class(boot, "muscal_bootstrap_result")
  expect_s3_class(perm, "muscal_permutation_result")
  expect_equal(boot$n_success, 4)
  expect_equal(perm$n_success, 4)
  expect_true(all(is.finite(boot$summary$mean)))
  expect_true(all(perm$component_results$p_value >= 0 & perm$component_results$p_value <= 1))
})

test_that("infer_muscal supports aligned_rrr bootstrap and permutation summaries", {
  set.seed(52)
  X <- list(
    X1 = matrix(rnorm(28 * 6), nrow = 28),
    X2 = matrix(rnorm(24 * 5), nrow = 24)
  )
  B <- qr.Q(qr(matrix(rnorm(3 * 2), 3, 2)))
  W1 <- matrix(rnorm(6 * 2), 6, 2)
  W2 <- matrix(rnorm(5 * 2), 5, 2)
  Y <- list(
    X1 = X$X1 %*% W1 %*% t(B) + matrix(rnorm(28 * 3, sd = 0.05), 28, 3),
    X2 = X$X2 %*% W2 %*% t(B) + matrix(rnorm(24 * 3, sd = 0.05), 24, 3)
  )

  fit <- aligned_rrr(
    Y = Y,
    X = X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    response_weights = c(X1 = 0.8, X2 = 1.2),
    ridge = 1e-8,
    max_iter = 60,
    tol = 1e-8
  )

  boot <- infer_muscal(fit, method = "bootstrap", statistic = "sdev", nrep = 4, seed = 31)
  perm <- infer_muscal(fit, method = "permutation", statistic = "sdev", nrep = 4, seed = 31)

  expect_s3_class(boot, "muscal_bootstrap_result")
  expect_s3_class(perm, "muscal_permutation_result")
  expect_equal(boot$n_success, 4)
  expect_equal(perm$n_success, 4)
  expect_true(all(is.finite(boot$summary$mean)))
  expect_true(all(perm$component_results$p_value >= 0 & perm$component_results$p_value <= 1))
})

test_that("infer_muscal supports ipca bootstrap and permutation summaries", {
  set.seed(40)
  X1 <- matrix(rnorm(60), nrow = 12)
  X2 <- matrix(rnorm(48), nrow = 12)

  fit <- ipca(
    list(X1 = X1, X2 = X2),
    ncomp = 2,
    lambda = 1,
    method = "gram",
    max_iter = 40,
    tol = 1e-5
  )

  boot <- infer_muscal(fit, method = "bootstrap", statistic = "sdev", nrep = 4, seed = 23)
  perm <- infer_muscal(fit, method = "permutation", statistic = "sdev", nrep = 5, seed = 23)

  expect_s3_class(boot, "muscal_bootstrap_result")
  expect_s3_class(perm, "muscal_permutation_result")
  expect_equal(boot$n_success, 4)
  expect_equal(perm$n_success, 5)
  expect_true(all(is.finite(boot$summary$mean)))
  expect_true(all(perm$component_results$p_value >= 0 & perm$component_results$p_value <= 1))
})

test_that("infer_muscal supports mcca bootstrap and permutation summaries", {
  set.seed(41)
  X1 <- matrix(rnorm(60), nrow = 12)
  X2 <- matrix(rnorm(48), nrow = 12)

  fit <- mcca(
    list(X1 = X1, X2 = X2),
    ncomp = 2,
    ridge = 1e-6
  )

  boot <- infer_muscal(fit, method = "bootstrap", statistic = "sdev", nrep = 4, seed = 29)
  perm <- infer_muscal(fit, method = "permutation", statistic = "sdev", nrep = 5, seed = 29)

  expect_s3_class(boot, "muscal_bootstrap_result")
  expect_s3_class(perm, "muscal_permutation_result")
  expect_equal(boot$n_success, 4)
  expect_equal(perm$n_success, 5)
  expect_true(all(is.finite(boot$summary$mean)))
  expect_true(all(perm$component_results$p_value >= 0 & perm$component_results$p_value <= 1))
})

test_that("infer_muscal aligned_mfa permutation distinguishes shared signal from null", {
  set.seed(38)
  N <- 35
  k <- 2
  S <- scale(matrix(rnorm(N * k), nrow = N, ncol = k), center = TRUE, scale = FALSE)
  idx <- list(
    X1 = sample.int(N, 28, replace = TRUE),
    X2 = sample.int(N, 27, replace = TRUE)
  )

  A1 <- matrix(rnorm(9 * k), nrow = 9, ncol = k)
  A2 <- matrix(rnorm(8 * k), nrow = 8, ncol = k)
  X1_signal <- S[idx$X1, , drop = FALSE] %*% t(A1) + matrix(rnorm(length(idx$X1) * 9, sd = 0.02), nrow = length(idx$X1), ncol = 9)
  X2_signal <- S[idx$X2, , drop = FALSE] %*% t(A2) + matrix(rnorm(length(idx$X2) * 8, sd = 0.02), nrow = length(idx$X2), ncol = 8)

  fit_signal <- aligned_mfa(
    list(X1 = X1_signal, X2 = X2_signal),
    row_index = idx,
    N = N,
    ncomp = 2,
    preproc = multivarious::pass(),
    normalization = "None"
  )
  stat_fn <- function(fit) as.numeric(mean(fit$block_fit$r2))
  perm_signal <- infer_muscal(
    fit_signal,
    method = "permutation",
    statistic_fn = stat_fn,
    nrep = 29,
    seed = 17
  )

  X1_null <- matrix(rnorm(length(idx$X1) * 9), nrow = length(idx$X1), ncol = 9)
  X2_null <- matrix(rnorm(length(idx$X2) * 8), nrow = length(idx$X2), ncol = 8)
  fit_null <- aligned_mfa(
    list(X1 = X1_null, X2 = X2_null),
    row_index = idx,
    N = N,
    ncomp = 2,
    preproc = multivarious::pass(),
    normalization = "None"
  )
  perm_null <- infer_muscal(
    fit_null,
    method = "permutation",
    statistic_fn = stat_fn,
    nrep = 29,
    seed = 17
  )

  expect_lt(perm_signal$component_results$p_value[[1]], 0.1)
  expect_gt(
    perm_null$component_results$p_value[[1]],
    perm_signal$component_results$p_value[[1]]
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
