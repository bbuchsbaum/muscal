library(testthat)
library(muscal)

set.seed(42)
# create a simple set of covariance matrices
gen_cov <- function(seed) {
  set.seed(seed)
  M <- matrix(rnorm(16), 4, 4)
  tcrossprod(M) / 4
}
subject_data <- lapply(1:5, gen_cov)
res <- covstatis(subject_data, ncomp = 2, norm_method = "none", dcenter = FALSE)

covar <- rnorm(5)

test_that("dimension mode returns named numeric vector", {
  dim_cos <- project_covariate(res, covar, what = "dimension", scale = "cosine")
  expect_type(dim_cos, "double")
  expect_length(dim_cos, 2)
  expect_equal(names(dim_cos), c("Dim1", "Dim2"))
  expect_true(all(abs(dim_cos) <= 1))
})

test_that("observation mode returns matrix with labels", {
  roi_beta <- project_covariate(res, covar, what = "observation", scale = "beta")
  expect_equal(dim(roi_beta), c(4, 2))
  expect_equal(rownames(roi_beta), res$labels)
  expect_equal(colnames(roi_beta), c("Dim1", "Dim2"))
})
