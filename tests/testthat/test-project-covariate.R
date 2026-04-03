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
feature_covar <- rnorm(4)

test_that("table-side dimension mode returns named numeric vector", {
  dim_cos <- project_table_covariate(res, covar, what = "dimension", scale = "cosine")
  expect_type(dim_cos, "double")
  expect_length(dim_cos, 2)
  expect_equal(names(dim_cos), c("Dim1", "Dim2"))
  expect_true(all(abs(dim_cos) <= 1))
})

test_that("table-side feature mode returns matrix with labels", {
  feature_beta <- project_table_covariate(res, covar, what = "feature", scale = "beta")
  expect_equal(dim(feature_beta), c(4, 2))
  expect_equal(rownames(feature_beta), res$labels)
  expect_equal(colnames(feature_beta), c("Dim1", "Dim2"))
})

test_that("feature-side projection returns named numeric vector", {
  dim_cos <- project_feature_covariate(res, feature_covar, scale = "cosine")
  expect_type(dim_cos, "double")
  expect_length(dim_cos, 2)
  expect_equal(names(dim_cos), c("Dim1", "Dim2"))
  expect_true(all(abs(dim_cos) <= 1))
})

test_that("generic project_covariate remains backward compatible for table-side input", {
  expect_equal(
    project_covariate(res, covar, what = "dimension", scale = "cosine"),
    project_table_covariate(res, covar, what = "dimension", scale = "cosine")
  )
  expect_equal(
    project_covariate(res, covar, what = "observation", scale = "beta"),
    project_table_covariate(res, covar, what = "feature", scale = "beta")
  )
  expect_equal(
    project_covariate(res, feature_covar, side = "feature", scale = "cosine"),
    project_feature_covariate(res, feature_covar, scale = "cosine")
  )
})
