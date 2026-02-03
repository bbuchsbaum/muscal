library(testthat)
library(muscal)

# Tests for synthetic_multiblock()

test_that("synthetic_multiblock creates correct structure", {
  set.seed(1)
  res <- synthetic_multiblock(S = 3, n = 10, p = 5, r = 2, sigma = 0, sphere = FALSE)

  expect_type(res, "list")
  expect_named(res, c("data_list", "coords_list", "V_true", "F_true", "Sadj"))

  expect_length(res$data_list, 3)
  expect_null(res$coords_list)
  expect_null(res$Sadj)
  expect_equal(dim(res$F_true), c(10, 2))
})

test_that("synthetic_multiblock data equals F_true %*% t(V_true) when sigma=0", {
  set.seed(1)
  res <- synthetic_multiblock(S = 3, n = 10, p = 5, r = 2, sigma = 0, sphere = FALSE)

  for (i in 1:3) {
    expect_equal(dim(res$V_true[[i]]), c(5, 2))
    expect_equal(res$data_list[[i]], res$F_true %*% t(res$V_true[[i]]))
  }
})

test_that("synthetic_multiblock data is column-centered", {
  set.seed(1)
  res <- synthetic_multiblock(S = 3, n = 10, p = 5, r = 2, sigma = 0, sphere = FALSE)

  for (i in 1:3) {
    col_means <- colMeans(res$data_list[[i]])
    expect_equal(col_means, rep(0, 5), tolerance = 1e-10)
  }
})

test_that("synthetic_multiblock V_true are orthonormal", {
  set.seed(1)
  res <- synthetic_multiblock(S = 3, n = 10, p = 5, r = 2, sigma = 0, sphere = FALSE)

  for (i in 1:3) {
    cross <- crossprod(res$V_true[[i]])
    expect_equal(cross, diag(2), tolerance = 1e-12)
  }
})

test_that("synthetic_multiblock F_true is orthonormal", {
  set.seed(1)
  res <- synthetic_multiblock(S = 3, n = 10, p = 5, r = 2, sigma = 0, sphere = FALSE)

  cross <- crossprod(res$F_true)
  expect_equal(cross, diag(2), tolerance = 1e-12)
})

test_that("synthetic_multiblock respects seed for reproducibility", {
  res1 <- synthetic_multiblock(S = 2, n = 10, p = 5, r = 2, seed = 42)
  res2 <- synthetic_multiblock(S = 2, n = 10, p = 5, r = 2, seed = 42)

  expect_equal(res1$F_true, res2$F_true)
  expect_equal(res1$V_true, res2$V_true)
  expect_equal(res1$data_list, res2$data_list)
})

test_that("synthetic_multiblock with different seeds gives different results", {
  res1 <- synthetic_multiblock(S = 2, n = 10, p = 5, r = 2, seed = 1)
  res2 <- synthetic_multiblock(S = 2, n = 10, p = 5, r = 2, seed = 2)

  expect_false(identical(res1$F_true, res2$F_true))
  expect_false(identical(res1$data_list, res2$data_list))
})

test_that("synthetic_multiblock handles vector p argument", {
  set.seed(1)
  res <- synthetic_multiblock(S = 3, n = 10, p = c(5, 8, 10), r = 2, sigma = 0)

  expect_equal(ncol(res$data_list[[1]]), 5)
  expect_equal(ncol(res$data_list[[2]]), 8)
  expect_equal(ncol(res$data_list[[3]]), 10)
  expect_equal(nrow(res$V_true[[1]]), 5)
  expect_equal(nrow(res$V_true[[2]]), 8)
  expect_equal(nrow(res$V_true[[3]]), 10)
})

test_that("synthetic_multiblock with sigma adds noise", {
  set.seed(1)
  res_clean <- synthetic_multiblock(S = 2, n = 10, p = 5, r = 2, sigma = 0)
  res_noisy <- synthetic_multiblock(S = 2, n = 10, p = 5, r = 2, sigma = 0.5, seed = 1)

  # Data should differ when sigma > 0
  expected_clean <- res_clean$F_true %*% t(res_clean$V_true[[1]])
  actual_noisy <- res_noisy$data_list[[1]]

  # The difference should be the noise
  noise <- actual_noisy - res_noisy$F_true %*% t(res_noisy$V_true[[1]])
  expect_true(sd(noise) > 0)
})

test_that("synthetic_multiblock with sphere=TRUE generates spatial coordinates", {
  skip_if_not_installed("RANN")
  skip_if_not_installed("Matrix")

  set.seed(42)
  res <- synthetic_multiblock(S = 2, n = 20, p = 10, r = 2, sigma = 0,
                               sphere = TRUE, k_nn = 5)

  # coords_list should be populated
  expect_length(res$coords_list, 2)
  expect_equal(dim(res$coords_list[[1]]), c(10, 3))  # p x 3 for 3D sphere
  expect_equal(dim(res$coords_list[[2]]), c(10, 3))

  # Sadj (spatial adjacency/Laplacian) should be a sparse matrix for all blocks combined
  expect_s4_class(res$Sadj, "Matrix")
  # Total dimension is sum of all p values
  expect_equal(nrow(res$Sadj), 20)  # 10 + 10 = 20
  expect_equal(ncol(res$Sadj), 20)
})

test_that("synthetic_multiblock sphere coordinates are on unit sphere", {
  skip_if_not_installed("RANN")
  skip_if_not_installed("Matrix")

  set.seed(42)
  res <- synthetic_multiblock(S = 2, n = 20, p = 10, r = 2, sigma = 0,
                               sphere = TRUE, k_nn = 5)

  for (i in 1:2) {
    coords <- res$coords_list[[i]]
    # Each row should have norm = 1 (on unit sphere)
    norms <- apply(coords, 1, function(x) sqrt(sum(x^2)))
    expect_equal(norms, rep(1, nrow(coords)), tolerance = 1e-10)
  }
})

test_that("synthetic_multiblock sphere Laplacian has correct structure", {
  skip_if_not_installed("RANN")
  skip_if_not_installed("Matrix")

  set.seed(42)
  res <- synthetic_multiblock(S = 2, n = 20, p = 10, r = 2, sigma = 0,
                               sphere = TRUE, k_nn = 3)

  # Laplacian should be symmetric
  expect_true(Matrix::isSymmetric(res$Sadj))
  # Diagonal should be positive (degree of each vertex)
  expect_true(all(Matrix::diag(res$Sadj) >= 0))
})

test_that("synthetic_multiblock with r=1 creates rank-1 data", {
  set.seed(1)
  res <- synthetic_multiblock(S = 2, n = 10, p = 5, r = 1, sigma = 0)

  expect_equal(dim(res$F_true), c(10, 1))
  expect_equal(dim(res$V_true[[1]]), c(5, 1))
  expect_equal(dim(res$V_true[[2]]), c(5, 1))

  # Data matrix should have rank 1 (all singular values except first are 0)
  svd_res <- svd(res$data_list[[1]])
  expect_equal(svd_res$d[2:5], rep(0, 4), tolerance = 1e-10)
})

test_that("synthetic_multiblock with large r", {
  set.seed(1)
  res <- synthetic_multiblock(S = 2, n = 20, p = 10, r = 5, sigma = 0)

  expect_equal(dim(res$F_true), c(20, 5))
  expect_equal(dim(res$V_true[[1]]), c(10, 5))

  # Check orthonormality still holds
  expect_equal(crossprod(res$F_true), diag(5), tolerance = 1e-10)
  expect_equal(crossprod(res$V_true[[1]]), diag(5), tolerance = 1e-10)
})

test_that("synthetic_multiblock output can be used with mfa", {
  set.seed(1)
  res <- synthetic_multiblock(S = 3, n = 20, p = 10, r = 3, sigma = 0.1)

  # Should be able to fit mfa on the generated data
  fit <- mfa(res$data_list, ncomp = 3)

  expect_s3_class(fit, "mfa")
  expect_length(fit$block_indices, 3)
})
