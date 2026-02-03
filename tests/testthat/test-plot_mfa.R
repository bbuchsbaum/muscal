library(testthat)
library(muscal)

# Tests for plot_mfa.R - smoke tests for plotting functions
# These tests verify that plotting functions execute without error

# --- Setup fixture data ------------------------------------------------------

# Create once and reuse - more efficient than recreating in each test
local_mfa_fixture <- function(env = parent.frame()) {
  set.seed(42)
  blocks <- replicate(3, matrix(rnorm(100), nrow = 20), simplify = FALSE)
  names(blocks) <- c("Block1", "Block2", "Block3")
  fit <- mfa(blocks, ncomp = 3)
  assign("mfa_fit", fit, envir = env)
  assign("blocks", blocks, envir = env)
}

local_linked_mfa_fixture <- function(env = parent.frame()) {
  set.seed(123)
  n <- 30
  Y <- matrix(rnorm(n * 10), n, 10)
  X1 <- matrix(rnorm(20 * 8), 20, 8)
  X2 <- matrix(rnorm(25 * 6), 25, 6)
  idx1 <- sample(1:n, 20)
  idx2 <- sample(1:n, 25)

  fit <- linked_mfa(Y, list(X1, X2), row_index = list(idx1, idx2), ncomp = 3)
  assign("linked_fit", fit, envir = env)
}

# --- plot.mfa tests ----------------------------------------------------------

test_that("plot.mfa runs without error", {
  local_mfa_fixture()

  # Suppress graphical output
  pdf(NULL)
  on.exit(dev.off())

  expect_silent(plot(mfa_fit))
  expect_silent(plot(mfa_fit, dims = c(1, 3)))
  expect_silent(plot(mfa_fit, labels = TRUE))
  expect_silent(plot(mfa_fit, col = "red", pch = 1, cex = 0.5))
})

test_that("plot.mfa returns score matrix invisibly", {
  local_mfa_fixture()

  pdf(NULL)
  on.exit(dev.off())

  result <- plot(mfa_fit, dims = c(1, 2))
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 20)
})

# --- plot.linked_mfa tests ---------------------------------------------------

test_that("plot.linked_mfa runs without error", {
  local_linked_mfa_fixture()

  pdf(NULL)
  on.exit(dev.off())

  expect_silent(plot(linked_fit))
  expect_silent(plot(linked_fit, show_coverage = TRUE))
  expect_silent(plot(linked_fit, show_coverage = FALSE))
  expect_silent(plot(linked_fit, labels = TRUE))
  expect_silent(plot(linked_fit, col = "blue"))
})

# --- autoplot.mfa tests ------------------------------------------------------

test_that("autoplot.mfa produces ggplot object", {
  skip_if_not_installed("ggplot2")
  local_mfa_fixture()

  p <- ggplot2::autoplot(mfa_fit)
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.mfa with various options", {
  skip_if_not_installed("ggplot2")
  local_mfa_fixture()

  p1 <- ggplot2::autoplot(mfa_fit, dims = c(1, 3))
  expect_s3_class(p1, "ggplot")

  p2 <- ggplot2::autoplot(mfa_fit, labels = TRUE)
  expect_s3_class(p2, "ggplot")

  color_vec <- rep(c("A", "B"), 10)
  p3 <- ggplot2::autoplot(mfa_fit, color = color_vec)
  expect_s3_class(p3, "ggplot")
})

# --- autoplot.linked_mfa tests -----------------------------------------------

test_that("autoplot.linked_mfa produces ggplot object", {
  skip_if_not_installed("ggplot2")
  local_linked_mfa_fixture()

  p <- ggplot2::autoplot(linked_fit)
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.linked_mfa with coverage coloring", {
  skip_if_not_installed("ggplot2")
  local_linked_mfa_fixture()

  p <- ggplot2::autoplot(linked_fit, color = "coverage")
  expect_s3_class(p, "ggplot")
})

# --- plot_partial_scores tests -----------------------------------------------

test_that("plot_partial_scores.mfa runs without error", {
  skip_if_not_installed("ggplot2")
  local_mfa_fixture()

  # MFA should have partial_scores if it's a recent version
  skip_if(is.null(mfa_fit$partial_scores), "partial_scores not available")

  p <- plot_partial_scores(mfa_fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_partial_scores.mfa with options", {
  skip_if_not_installed("ggplot2")
  local_mfa_fixture()

  skip_if(is.null(mfa_fit$partial_scores), "partial_scores not available")

  p1 <- plot_partial_scores(mfa_fit, dims = c(1, 3))
  expect_s3_class(p1, "ggplot")

  p2 <- plot_partial_scores(mfa_fit, connect = FALSE)
  expect_s3_class(p2, "ggplot")

  p3 <- plot_partial_scores(mfa_fit, show_consensus = FALSE)
  expect_s3_class(p3, "ggplot")

  p4 <- plot_partial_scores(mfa_fit, blocks = 1:2)
  expect_s3_class(p4, "ggplot")
})

test_that("plot_partial_scores.mfa errors when partial_scores not available", {
  local_mfa_fixture()

  mfa_no_partial <- mfa_fit
  mfa_no_partial$partial_scores <- NULL

  expect_error(
    plot_partial_scores(mfa_no_partial),
    "Partial scores are not available"
  )
})

test_that("plot_partial_scores.linked_mfa runs without error", {
  skip_if_not_installed("ggplot2")
  local_linked_mfa_fixture()

  skip_if(is.null(linked_fit$partial_scores), "partial_scores not available")

  p <- plot_partial_scores(linked_fit)
  expect_s3_class(p, "ggplot")
})

# --- plot_block_weights tests ------------------------------------------------

test_that("plot_block_weights.mfa runs without error", {
  skip_if_not_installed("ggplot2")
  local_mfa_fixture()

  p <- plot_block_weights(mfa_fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_block_weights.linked_mfa runs without error", {
  skip_if_not_installed("ggplot2")
  local_linked_mfa_fixture()

  p <- plot_block_weights(linked_fit)
  expect_s3_class(p, "ggplot")
})

# --- plot_block_similarity tests ---------------------------------------------

test_that("plot_block_similarity.mfa runs without error", {
  skip_if_not_installed("ggplot2")
  local_mfa_fixture()

  p1 <- plot_block_similarity(mfa_fit, method = "RV")
  expect_s3_class(p1, "ggplot")

  p2 <- plot_block_similarity(mfa_fit, method = "RV2")
  expect_s3_class(p2, "ggplot")
})

# --- plot_variance tests -----------------------------------------------------

test_that("plot_variance.mfa with bar plot", {
  skip_if_not_installed("ggplot2")
  local_mfa_fixture()

  p <- plot_variance(mfa_fit, type = "bar")
  expect_s3_class(p, "ggplot")
})

test_that("plot_variance.mfa with line plot", {
  skip_if_not_installed("ggplot2")
  local_mfa_fixture()

  p <- plot_variance(mfa_fit, type = "line")
  expect_s3_class(p, "ggplot")
})

# --- plot_loadings tests -----------------------------------------------------

test_that("plot_loadings.mfa with circle type", {
  skip_if_not_installed("ggplot2")
  local_mfa_fixture()

  p <- plot_loadings(mfa_fit, type = "circle")
  expect_s3_class(p, "ggplot")
})

test_that("plot_loadings.mfa with bar type", {
  skip_if_not_installed("ggplot2")
  local_mfa_fixture()

  p <- plot_loadings(mfa_fit, type = "bar")
  expect_s3_class(p, "ggplot")
})

test_that("plot_loadings.mfa with different dims", {
  skip_if_not_installed("ggplot2")
  local_mfa_fixture()

  p <- plot_loadings(mfa_fit, dims = c(1, 3))
  expect_s3_class(p, "ggplot")
})

test_that("plot_loadings.linked_mfa runs without error", {
  skip_if_not_installed("ggplot2")
  local_linked_mfa_fixture()

  p <- plot_loadings(linked_fit)
  expect_s3_class(p, "ggplot")
})

# --- plot_convergence tests --------------------------------------------------

test_that("plot_convergence.linked_mfa runs without error", {
  skip_if_not_installed("ggplot2")
  local_linked_mfa_fixture()

  # Only run if convergence history is available
  skip_if(is.null(linked_fit$convergence) && is.null(linked_fit$objective_history),
          "convergence info not available")

  p <- plot_convergence(linked_fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_convergence with log_scale option", {
  skip_if_not_installed("ggplot2")
  local_linked_mfa_fixture()

  skip_if(is.null(linked_fit$convergence) && is.null(linked_fit$objective_history),
          "convergence info not available")

  p <- plot_convergence(linked_fit, log_scale = TRUE)
  expect_s3_class(p, "ggplot")
})

# --- plot_coverage tests -----------------------------------------------------

test_that("plot_coverage.linked_mfa runs without error", {
  skip_if_not_installed("ggplot2")
  local_linked_mfa_fixture()

  p <- plot_coverage(linked_fit)
  expect_s3_class(p, "ggplot")
})

# --- .plot_scores_base internal helper ---------------------------------------

test_that(".plot_scores_base handles missing sdev", {
  S <- matrix(rnorm(40), nrow = 20, ncol = 2)

  pdf(NULL)
  on.exit(dev.off())

  result <- muscal:::.plot_scores_base(
    S, dims = c(1, 2), labels = FALSE, col = "blue",
    pch = 19, cex = 1, main = "Test", sdev = NULL
  )

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(20, 2))
})

test_that(".plot_scores_base with labels uses rownames", {
  S <- matrix(rnorm(40), nrow = 20, ncol = 2)
  rownames(S) <- paste0("obs", 1:20)

  pdf(NULL)
  on.exit(dev.off())

  # Should not error with labels = TRUE
  expect_silent(
    muscal:::.plot_scores_base(
      S, dims = c(1, 2), labels = TRUE, col = "blue",
      pch = 19, cex = 1, main = "Test", sdev = c(2, 1)
    )
  )
})

# --- .autoplot_scores internal helper ----------------------------------------

test_that(".autoplot_scores creates valid ggplot", {
  skip_if_not_installed("ggplot2")

  S <- matrix(rnorm(40), nrow = 20, ncol = 2)
  rownames(S) <- paste0("obs", 1:20)

  p <- muscal:::.autoplot_scores(
    S, dims = c(1, 2), labels = FALSE, color = NULL,
    alpha = 0.8, size = 2, sdev = c(2, 1), title = "Test"
  )

  expect_s3_class(p, "ggplot")
})

test_that(".autoplot_scores with color variable", {
  skip_if_not_installed("ggplot2")

  S <- matrix(rnorm(40), nrow = 20, ncol = 2)
  color_var <- factor(rep(c("A", "B"), 10))

  p <- muscal:::.autoplot_scores(
    S, dims = c(1, 2), labels = FALSE, color = color_var,
    alpha = 0.8, size = 2, sdev = c(2, 1), title = "Test"
  )

  expect_s3_class(p, "ggplot")
})

test_that(".autoplot_scores with labels", {
  skip_if_not_installed("ggplot2")

  S <- matrix(rnorm(40), nrow = 20, ncol = 2)

  p <- muscal:::.autoplot_scores(
    S, dims = c(1, 2), labels = TRUE, color = NULL,
    alpha = 0.8, size = 2, sdev = c(2, 1), title = "Test"
  )

  expect_s3_class(p, "ggplot")
})

test_that(".autoplot_scores with NULL sdev uses NA variance", {
  skip_if_not_installed("ggplot2")

  S <- matrix(rnorm(40), nrow = 20, ncol = 2)

  p <- muscal:::.autoplot_scores(
    S, dims = c(1, 2), labels = FALSE, color = NULL,
    alpha = 0.8, size = 2, sdev = NULL, title = "Test"
  )

  expect_s3_class(p, "ggplot")
})
