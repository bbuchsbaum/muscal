library(testthat)
library(musca)

# Check for informative error when coords_list blocks do not have 3 columns

test_that("coords_list column count is validated", {
  dl <- list(matrix(0, 2, 2), matrix(0, 2, 3))
  cl <- list(matrix(0, ncol(dl[[1]]), 3), matrix(0, ncol(dl[[2]]), 4))
  expect_error(
    penalized_mfa_clusterwise(dl, cl, ncomp = 1),
    "coords_list element 2 must have exactly 3 columns"
  )
})
