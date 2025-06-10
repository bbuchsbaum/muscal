library(testthat)
library(musca)

# Tests for bamfa functions


test_that("bamfa.list assigns default block names", {
  blocks <- list(matrix(0, 2, 2), matrix(1, 2, 2))
  names(blocks) <- NULL
  res <- bamfa(blocks, k_g = 1, k_l = 1, niter = 1, preproc = multivarious::pass())
  expect_equal(names(res$block_indices), c("Block_1", "Block_2"))
})


test_that("bamfa.multiblock checks row consistency", {
  skip_if_not_installed("multidesign")
  mb <- multidesign::multiblock(list(A = matrix(0, 2, 2), B = matrix(0, 3, 2)))
  expect_error(bamfa(mb, k_g = 1, k_l = 1, niter = 1, preproc = multivarious::pass()))
})


test_that("bamfa.multidesign requires subject argument", {
  skip_if_not_installed("multidesign")
  x <- matrix(rnorm(10), nrow = 5)
  design <- data.frame(subj = rep(c("S1", "S2"), each = 5))
  md <- multidesign::multidesign(x, design)
  expect_error(bamfa(md, k_g = 1, k_l = 1, niter = 1, preproc = multivarious::pass()))
})


test_that("bamfa.list returns projector with expected structure", {
  set.seed(123)
  blocks <- list(A = matrix(rnorm(20), 5, 4), B = matrix(rnorm(20), 5, 4))
  res <- bamfa(blocks, k_g = 2, k_l = 1, niter = 2, preproc = multivarious::pass())
  expect_s3_class(res, "bamfa")
  expect_s3_class(res, "multiblock_projector")
  expect_equal(res$k_g, 2)
  expect_equal(res$k_l, 1)
  expect_length(res$B_list, 2)
  expect_length(res$block_indices, 2)
})

