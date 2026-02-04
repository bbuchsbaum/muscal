library(testthat)
library(muscal)

context("Duality in COVSTATIS")

mats <- replicate(6, { A <- matrix(rnorm(40), 5, 8); tcrossprod(A) },
                  simplify = FALSE)

test_that("duality identity holds on random data", {
  set.seed(11)
  
  mod  <- covstatis(mats, ncomp = 4)
  
  # Test helper functions
  expect_true(check_duality(mod))
  
  # Verify energy equality dimension-wise
  F_scores  <- mod$roi_scores
  T_scores  <- muscal:::rv_subject_scores(mod)
  expect_equal(colSums(F_scores^2), colSums(T_scores^2), tolerance = 1e-10)
})

test_that("MFA normalisation preserves duality", {
  set.seed(11)
  mod <- covstatis(mats, ncomp = 4, norm_method = "mfa")
  expect_true(check_duality(mod))
}) 