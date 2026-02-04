test_that("Adam vs gradient give similar objective drop in penalized_mfa", {
  Xs <- lapply(1:3, function(i) scale(matrix(rnorm(300), 30, 10), scale = FALSE))
  
  # Test with gradient optimizer
  fit_g <- penalized_mfa(Xs, ncomp = 2, lambda = 0.5, optimizer = "gradient", max_iter = 10)
  obj_g <- attr(fit_g, "obj_values")
  
  # Test with adam optimizer
  fit_a <- penalized_mfa(Xs, ncomp = 2, lambda = 0.5, optimizer = "adam", max_iter = 10)
  obj_a <- attr(fit_a, "obj_values")
  
  # Check that objective function decreases for both
  expect_lt(tail(obj_g, 1), head(obj_g, 1))
  expect_lt(tail(obj_a, 1), head(obj_a, 1))
  
  # Check that the final objective values are in a similar ballpark (not a strict test)
  expect_true(abs(tail(obj_g, 1) - tail(obj_a, 1)) / tail(obj_g, 1) < 0.5)
})

test_that("penalized_mfa projection penalty works", {
  Xs <- lapply(1:3, function(i) scale(matrix(rnorm(300), 30, 10), scale = FALSE))
  
  # Test with projection penalty
  fit_p <- penalized_mfa(Xs, ncomp = 2, lambda = 0.5, penalty_method = "projection", max_iter = 10)
  obj_p <- attr(fit_p, "obj_values")
  
  # Check that objective function decreases
  expect_lt(tail(obj_p, 1), head(obj_p, 1))
  
  # Check that it returns the correct class
  expect_s3_class(fit_p, "penalized_mfa")
  expect_s3_class(fit_p, "multiblock_projector")
})

test_that("covstatis duality check from test sketch works", {
  # This test was included in the user's feedback sketch.
  # It checks a property of covstatis, not penalized_mfa, but is included for completeness.
  Xs <- lapply(1:5, function(i) matrix(rnorm(10*10), 10, 10))
  Xs_cov <- lapply(Xs, cov)
  
  # covstatis requires symmetric matrices, cor() or cov() ensures this.
  res <- covstatis(Xs_cov, ncomp=3)
  expect_true(check_duality(res))
}) 