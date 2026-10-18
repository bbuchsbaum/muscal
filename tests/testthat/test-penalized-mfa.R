library(testthat)
library(muscal)

# Test data setup
set.seed(123)
X1 <- scale(matrix(rnorm(100 * 4), 100, 4), scale = FALSE)
X2 <- scale(matrix(rnorm(100 * 16), 100, 16), scale = FALSE)
data_list_uneven <- list(B1 = X1, B2 = X2)
data_list_even <- list(B1 = X1, B2 = X1)

test_that("penalized_mfa.list handles uneven and even blocks correctly", {
  # Test with uneven blocks: should warn and produce no consensus matrix
  expect_warning(
    res_uneven <- penalized_mfa(data_list_uneven, ncomp = 2, lambda = 1, compute_consensus = TRUE, verbose = FALSE, penalty_method = "projection"),
    "Cannot compute consensus"
  )
  expect_s3_class(res_uneven, "penalized_mfa")
  expect_null(attr(res_uneven, "consensus"))

  # Test with even blocks: should produce a consensus matrix
  res_even <- penalized_mfa(data_list_even, ncomp = 2, lambda = 1, compute_consensus = TRUE, verbose = FALSE, penalty_method = "projection")
  expect_s3_class(res_even, "penalized_mfa")
  expect_true(!is.null(attr(res_even, "consensus")))
  expect_equal(dim(attr(res_even, "consensus")), c(4, 2))
})

test_that("penalized_mfa works with multiblock input", {
  mb <- multidesign::multiblock(data_list_uneven)
  # Should warn about disabling penalty when blocks have different dimensions
  expect_warning(
    res <- penalized_mfa(mb, ncomp = 2, lambda = 1, penalty_method = "projection"),
    "Cannot use 'projection' penalty with blocks of different dimensions"
  )
  expect_s3_class(res, "penalized_mfa")
})

test_that("penalized_mfa stops with fewer than 2 blocks", {
  expect_error(penalized_mfa(list(X1), ncomp = 2), "requires at least 2 data blocks")
})

test_that("penalized_mfa supports preproc = NULL after shared preprocessing consolidation", {
  fit <- penalized_mfa(
    data_list_even,
    preproc = NULL,
    ncomp = 2,
    lambda = 1,
    max_iter = 2,
    penalty_method = "projection"
  )

  expect_s3_class(fit, "penalized_mfa")
  expect_s3_class(fit$preproc, "pre_processor")
})

test_that("penalized_mfa rejects unequal observation counts before preprocessing", {
  blocks <- list(short = matrix(1:12, 3, 4), long = matrix(1:20, 5, 4))
  for (consensus in c(FALSE, TRUE)) {
    touched <- FALSE
    expect_error(penalized_mfa(
      blocks, ncomp = 1, compute_consensus = consensus,
      preproc = { touched <- TRUE; multivarious::center() }
    ), "equal positive row counts.*short=3.*long=5")
    expect_false(touched)
  }
  expect_error(penalized_mfa(list(matrix(numeric(), 0, 2), matrix(1, 3, 2))),
               "equal positive row counts")
})

test_that("penalized_mfa delegates shared-row validation from multidesign", {
  md <- multidesign::multidesign(matrix(seq_len(32), 8, 4),
                                 data.frame(subject = c(rep("a", 3), rep("b", 5))))
  expect_error(penalized_mfa(md, subject = subject, ncomp = 1),
               "equal positive row counts")
})

test_that("penalized_mfa rejects preprocessing that changes observation counts", {
  testthat::local_mocked_bindings(
    prepare_block_preprocessors = function(data_list, ...) {
      list(Xp = lapply(data_list, function(x) x[-1, , drop = FALSE]), proclist = NULL)
    }, .package = "muscal"
  )
  blocks <- list(matrix(1:12, 3, 4), matrix(13:24, 3, 4))
  expect_error(penalized_mfa(blocks, ncomp = 1), "preserve.*observation")
})
