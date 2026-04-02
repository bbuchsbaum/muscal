library(testthat)
library(muscal)

test_that("metric registry and defaults are task aware", {
  reg <- metric_registry()
  expect_s3_class(reg, "tbl_df")
  expect_true(all(c("task", "metric", "maximize", "description") %in% names(reg)))

  row_reg <- metric_registry("row_alignment")
  expect_true(all(row_reg$task == "retrieval_alignment"))
  expect_setequal(default_metrics("reconstruction"), c("mse", "rmse", "r2"))
  expect_setequal(default_metrics("row_alignment"), c("mean_top1_similarity", "mean_topk_similarity"))
})

test_that("performance_metrics computes reconstruction metrics", {
  truth <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  estimate <- matrix(c(1, 1, 3, 5), nrow = 2, byrow = TRUE)

  res <- performance_metrics(
    task = "reconstruction",
    truth = truth,
    estimate = estimate,
    metrics = c("mse", "rmse", "mae")
  )

  diff <- truth - estimate
  expect_equal(res$mse, mean(diff^2))
  expect_equal(res$rmse, sqrt(mean(diff^2)))
  expect_equal(res$mae, mean(abs(diff)))
})

test_that("performance_metrics computes response prediction metrics", {
  truth <- matrix(c(
    1, 2, 3,
    2, 4, 6
  ), nrow = 2, byrow = TRUE)
  estimate <- truth

  res <- performance_metrics(
    task = "response_prediction",
    truth = truth,
    estimate = estimate,
    metrics = c("mse", "mean_correlation", "mean_cosine_similarity")
  )

  expect_equal(res$mse, 0)
  expect_equal(res$mean_correlation, 1)
  expect_equal(res$mean_cosine_similarity, 1)
})

test_that("performance_metrics computes retrieval and alignment metrics", {
  truth <- matrix(c(
    1, 0,
    0, 1
  ), nrow = 2, byrow = TRUE)

  estimate <- list(
    retrieved_features = list(
      rbind(c(1, 0), c(0, 1)),
      rbind(c(1, 0), c(0, 1))
    ),
    retrieved_ids = list(
      c("a", "b"),
      c("a", "b")
    ),
    oracle_similarity = c(1, 1)
  )

  res <- performance_metrics(
    task = "row_alignment",
    truth = truth,
    estimate = estimate,
    metrics = c("mean_top1_similarity", "mean_topk_similarity", "oracle_gap", "recall_at_k", "mrr"),
    k = 2,
    truth_ids = c("a", "b")
  )

  expect_equal(res$mean_top1_similarity, 0.5)
  expect_equal(res$mean_topk_similarity, 0.5)
  expect_equal(res$oracle_gap, 0.5)
  expect_equal(res$recall_at_k, 1)
  expect_equal(res$mrr, 0.75)
})

test_that("performance metric outputs combine cleanly across folds", {
  truth <- diag(2)
  estimate1 <- truth
  estimate2 <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE)

  fold_scores <- dplyr::bind_rows(
    performance_metrics("response_prediction", truth, estimate1),
    performance_metrics("response_prediction", truth, estimate2)
  )

  expect_equal(nrow(fold_scores), 2)
  expect_true(all(default_metrics("response_prediction") %in% names(fold_scores)))
})

test_that("performance_metrics rejects invalid task and metric combinations", {
  truth <- diag(2)
  estimate <- truth

  expect_error(
    performance_metrics("response_prediction", truth, estimate, metrics = "mean_top1_similarity"),
    "Unsupported metrics"
  )

  expect_error(
    performance_metrics("row_alignment", truth, list(retrieved_features = list(truth[1, , drop = FALSE], truth[2, , drop = FALSE])), metrics = "mrr"),
    "truth_ids"
  )

  expect_error(
    performance_metrics("row_alignment", truth, list(retrieved_features = list(truth[1, , drop = FALSE], truth[2, , drop = FALSE])), metrics = "oracle_gap"),
    "oracle"
  )
})

test_that("performance_metrics aliases and retrieval encodings are equivalent", {
  truth <- matrix(c(
    1, 0,
    0, 1
  ), nrow = 2, byrow = TRUE)

  estimate_matrix <- truth
  estimate_list <- list(
    retrieved_features = list(
      truth[1, , drop = FALSE],
      truth[2, , drop = FALSE]
    )
  )

  res_alias <- performance_metrics("alignment", truth, estimate_matrix, metrics = "mean_top1_similarity")
  res_list <- performance_metrics("row_alignment", truth, estimate_list, metrics = "mean_top1_similarity")

  expect_equal(res_alias$mean_top1_similarity, 1)
  expect_equal(res_alias$mean_top1_similarity, res_list$mean_top1_similarity)
})

test_that("response similarity metrics are invariant to common positive scaling", {
  truth <- matrix(c(
    1, 2, 3,
    2, 3, 5
  ), nrow = 2, byrow = TRUE)
  estimate <- matrix(c(
    1.2, 1.9, 2.8,
    2.1, 2.8, 5.2
  ), nrow = 2, byrow = TRUE)

  base <- performance_metrics(
    "response_prediction",
    truth,
    estimate,
    metrics = c("mean_correlation", "mean_cosine_similarity")
  )
  scaled <- performance_metrics(
    "response_prediction",
    truth * 10,
    estimate * 10,
    metrics = c("mean_correlation", "mean_cosine_similarity")
  )

  expect_equal(base$mean_correlation, scaled$mean_correlation)
  expect_equal(base$mean_cosine_similarity, scaled$mean_cosine_similarity)
})

test_that("performance_metrics validates dimensions and handles zero-norm rows robustly", {
  truth <- rbind(c(0, 0), c(1, 0))
  estimate <- truth

  res <- performance_metrics(
    "response_prediction",
    truth,
    estimate,
    metrics = "mean_cosine_similarity"
  )

  expect_equal(res$mean_cosine_similarity, 1)

  expect_error(
    performance_metrics("response_prediction", truth, estimate[1, , drop = FALSE], metrics = "mse"),
    "shape"
  )

  expect_error(
    performance_metrics("row_alignment", truth, estimate[1, , drop = FALSE], metrics = "mean_top1_similarity"),
    "Top-1 retrieval estimate has 1 rows"
  )
})
