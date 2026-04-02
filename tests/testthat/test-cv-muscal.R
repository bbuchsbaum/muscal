library(testthat)
library(muscal)

test_that("cv_muscal supports explicit row holdout for reconstruction tasks", {
  set.seed(21)
  X1 <- matrix(rnorm(30), nrow = 10)
  X2 <- matrix(rnorm(20), nrow = 10)
  X <- cbind(X1, X2)
  design <- data.frame(group = rep(c("A", "B"), each = 5))
  md <- multidesign::multidesign(X, design)
  folds <- multidesign::cv_rows(md, rows = list(1:2, 6:7), preserve_row_ids = TRUE)

  fit_fn <- function(analysis) {
    Xa <- multidesign::xdata(analysis)
    mfa(
      list(X1 = Xa[, 1:ncol(X1), drop = FALSE], X2 = Xa[, (ncol(X1) + 1):ncol(Xa), drop = FALSE]),
      ncomp = 2
    )
  }

  res <- cv_muscal(
    folds = folds,
    fit_fn = fit_fn,
    estimate_fn = function(model, assessment) {
      predict(model, multidesign::xdata(assessment), type = "reconstruction")
    },
    truth_fn = function(assessment) multidesign::xdata(assessment),
    metrics = c("mse", "rmse")
  )

  expect_s3_class(res, "muscal_cv_result")
  expect_s3_class(res, "cv_result")
  expect_true(all(c(".fold", "mse", "rmse") %in% names(res$scores)))
  expect_equal(nrow(res$scores), 2)
  expect_false(is.null(res$foldframe))
  expect_equal(res$folds[[1]]$held_out$rows, 1:2)
})

test_that("cv_muscal supports grouped holdout for reconstruction tasks", {
  set.seed(22)
  X1 <- matrix(rnorm(36), nrow = 12)
  X2 <- matrix(rnorm(24), nrow = 12)
  X <- cbind(X1, X2)
  design <- data.frame(group = rep(c("A", "B", "C"), each = 4))
  md <- multidesign::multidesign(X, design)
  folds <- multidesign::fold_over(md, group, preserve_row_ids = TRUE)

  fit_fn <- function(analysis) {
    Xa <- multidesign::xdata(analysis)
    mfa(
      list(X1 = Xa[, 1:ncol(X1), drop = FALSE], X2 = Xa[, (ncol(X1) + 1):ncol(Xa), drop = FALSE]),
      ncomp = 2
    )
  }

  res <- cv_muscal(
    folds = folds,
    fit_fn = fit_fn,
    estimate_fn = function(model, assessment) {
      predict(model, multidesign::xdata(assessment), type = "reconstruction")
    },
    truth_fn = function(assessment) multidesign::xdata(assessment),
    metrics = c("mse")
  )

  expect_equal(nrow(res$scores), 3)
  expect_true(all(!is.na(res$scores$mse)))
  expect_true(".splitvar" %in% names(res$foldframe))
})

test_that("cv_muscal supports synchronized hyperdesign row holdout for anchored prediction", {
  set.seed(23)
  N <- 30
  Y <- matrix(rnorm(N * 3), N, 3)
  colnames(Y) <- c("y1", "y2", "y3")

  X1 <- matrix(rnorm(14 * 5), 14, 5)
  X2 <- matrix(rnorm(16 * 4), 16, 4)
  idx1 <- sample.int(N, nrow(X1), replace = FALSE)
  idx2 <- sample.int(N, nrow(X2), replace = FALSE)

  d1 <- multidesign::multidesign(
    X1,
    data.frame(anchor = idx1, Y[idx1, , drop = FALSE])
  )
  d2 <- multidesign::multidesign(
    X2,
    data.frame(anchor = idx2, Y[idx2, , drop = FALSE])
  )
  hd <- multidesign::hyperdesign(list(X1 = d1, X2 = d2), block_names = c("X1", "X2"))

  folds <- multidesign::cv_rows(
    hd,
    rows = list(
      list(X1 = 1:2, X2 = 1:2),
      list(X1 = 3:4, X2 = 3:4)
    ),
    preserve_row_ids = TRUE
  )

  fit_fn <- function(analysis) {
    X_blocks <- multidesign::xdata(analysis)
    idx_blocks <- lapply(multidesign::design(analysis), function(des) des$anchor)
    anchored_mfa(Y = Y, X = X_blocks, row_index = idx_blocks, ncomp = 2)
  }

  estimate_fn <- function(model, assessment) {
    X_blocks <- multidesign::xdata(assessment)
    do.call(rbind, lapply(names(X_blocks), function(block_name) {
      predict(model, X_blocks[[block_name]], block = block_name, type = "response")
    }))
  }

  truth_fn <- function(assessment) {
    des_blocks <- multidesign::design(assessment)
    do.call(rbind, lapply(des_blocks, function(des) {
      as.matrix(des[, c("y1", "y2", "y3"), drop = FALSE])
    }))
  }

  res <- cv_muscal(
    folds = folds,
    fit_fn = fit_fn,
    estimate_fn = estimate_fn,
    truth_fn = truth_fn,
    task = "response_prediction",
    metrics = c("mse", "mean_cosine_similarity")
  )

  expect_equal(nrow(res$scores), 2)
  expect_true(all(c("mse", "mean_cosine_similarity") %in% names(res$scores)))
  expect_true(".block" %in% names(res$foldframe))
  expect_true(inherits(res$folds[[1]]$assessment, "hyperdesign"))
})

test_that("cv_muscal infers the task from the fitted model and matches manual fold scoring", {
  set.seed(24)
  X1 <- matrix(rnorm(24), nrow = 8)
  X2 <- matrix(rnorm(16), nrow = 8)
  X <- cbind(X1, X2)
  md <- multidesign::multidesign(X, data.frame(group = rep(c("A", "B"), each = 4)))
  folds <- multidesign::cv_rows(md, rows = list(1:2, 5:6), preserve_row_ids = TRUE)

  fit_fn <- function(analysis) {
    Xa <- multidesign::xdata(analysis)
    mfa(
      list(X1 = Xa[, 1:ncol(X1), drop = FALSE], X2 = Xa[, (ncol(X1) + 1):ncol(Xa), drop = FALSE]),
      ncomp = 2
    )
  }

  estimate_fn <- function(model, assessment) {
    predict(model, multidesign::xdata(assessment), type = "reconstruction")
  }

  truth_fn <- function(assessment) multidesign::xdata(assessment)

  res <- cv_muscal(
    folds = folds,
    fit_fn = fit_fn,
    estimate_fn = estimate_fn,
    truth_fn = truth_fn
  )

  manual_model <- fit_fn(folds[[1]]$analysis)
  manual_truth <- truth_fn(folds[[1]]$assessment)
  manual_estimate <- estimate_fn(manual_model, folds[[1]]$assessment)
  manual_score <- performance_metrics(manual_model$task, manual_truth, manual_estimate)

  expect_true(all(default_metrics("reconstruction") %in% names(res$scores)))
  expect_equal(res$scores$mse[[1]], manual_score$mse[[1]])
  expect_equal(res$scores$rmse[[1]], manual_score$rmse[[1]])
  expect_equal(res$scores$r2[[1]], manual_score$r2[[1]])
})

test_that("cv_muscal errors when task cannot be resolved", {
  set.seed(25)
  X <- matrix(rnorm(20), nrow = 5)
  md <- multidesign::multidesign(X, data.frame(group = rep(1, 5)))
  folds <- multidesign::cv_rows(md, rows = list(1:2), preserve_row_ids = TRUE)

  expect_error(
    cv_muscal(
      folds = folds,
      fit_fn = function(analysis) list(model = "no-task"),
      estimate_fn = function(model, assessment) multidesign::xdata(assessment),
      truth_fn = function(assessment) multidesign::xdata(assessment),
      metrics = "mse"
    ),
    "does not expose `model\\$task`"
  )
})
