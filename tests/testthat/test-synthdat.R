library(testthat)
library(musca)

# Tests for synthetic_multiblock

set.seed(1)
res <- synthetic_multiblock(S = 3, n = 10, p = 5, r = 2, sigma = 0, sphere = FALSE)

# check basic structure
expect_equal(length(res$data_list), 3)
expect_null(res$coords_list)
expect_equal(dim(res$F_true), c(10,2))

# data matrices should equal F_true %*% t(V_true[[i]]) when sigma=0
for(i in 1:3) {
  expect_equal(dim(res$V_true[[i]]), c(5,2))
  expect_equal(res$data_list[[i]], res$F_true %*% t(res$V_true[[i]]))
  expect_equal(colMeans(res$data_list[[i]]), rep(0,5))
  expect_equal(crossprod(res$V_true[[i]]), diag(2), tolerance = 1e-12)
}

# orthonormality of shared factors
expect_equal(crossprod(res$F_true), diag(2), tolerance = 1e-12)
