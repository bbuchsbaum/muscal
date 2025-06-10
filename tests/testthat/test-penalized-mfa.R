library(testthat)
library(musca)

# --- should_precompute ------------------------------------------------------

test_that("should_precompute respects memory budget", {
  expect_true(musca:::should_precompute(p = 10, mem_mb = 1))
  expect_false(musca:::should_precompute(p = 1000, mem_mb = 1))
})

# --- make_grad_fun ---------------------------------------------------------

test_that("make_grad_fun with XtX and on-the-fly agree", {
  set.seed(1)
  X <- matrix(rnorm(6), nrow = 2)
  XtX <- crossprod(X)
  V <- matrix(rnorm(6), nrow = 3)
  grad_pre <- musca:::make_grad_fun(X, XtX)
  grad_fly <- musca:::make_grad_fun(X)
  expect_equal(grad_pre(V), grad_fly(V))
})

# --- penalized_mfa.list output --------------------------------------------


test_that("penalized_mfa.list returns a projector with expected attributes", {
  skip_if_not_installed("multivarious")
  dl <- list(matrix(rnorm(20), 5, 4), matrix(rnorm(20), 5, 4))
  res <- musca:::penalized_mfa.list(dl, ncomp = 2, lambda = 0.1,
                                     max_iter = 1, nsteps_inner = 1,
                                     learning_rate = 0.01,
                                     compute_consensus = TRUE)
  expect_s3_class(res, "multiblock_projector")
  expect_true(!is.null(attr(res, "obj_values")))
  expect_equal(ncol(res$v), 2)
  expect_length(attr(res, "precompute_info"), 2)
  expect_true(!is.null(attr(res, "consensus")))
})
