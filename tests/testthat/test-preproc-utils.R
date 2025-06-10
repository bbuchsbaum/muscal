library(testthat)
library(musca)

# Basic matrices for preprocessing tests
mat1 <- matrix(rnorm(20), nrow = 5)
mat2 <- matrix(rnorm(30), nrow = 5)

test_that("single preprocessor is replicated across blocks", {
  prep_def <- multivarious::center()
  res <- musca:::prepare_block_preprocessors(list(mat1, mat2), prep_def)
  expect_length(res$proclist, 2)
  expect_true(all(vapply(res$proclist, inherits, logical(1), "pre_processor")))
  expect_equal(lapply(res$Xp, ncol), list(ncol(mat1), ncol(mat2)))
})



test_that("list preprocessor length mismatch errors", {
  prep_list <- list(multivarious::center())
  expect_error(
    musca:::prepare_block_preprocessors(list(mat1, mat2), prep_list),
    "length"  # message about length mismatch
  )
})

test_that("inconsistent column counts warn when no preprocessing", {
  expect_warning(
    musca:::prepare_block_preprocessors(list(mat1, mat2[,1:3]), NULL),
    "Input blocks"
  )
})

