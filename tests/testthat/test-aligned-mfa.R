library(testthat)
library(muscal)

test_that("aligned_mfa validates row_index and returns expected structure", {
  set.seed(1)
  N <- 30
  X1 <- matrix(rnorm(10 * 6), 10, 6)
  X2 <- matrix(rnorm(12 * 5), 12, 5)
  idx1 <- sample.int(N, nrow(X1), replace = FALSE)
  idx2 <- sample.int(N, nrow(X2), replace = FALSE)

  expect_error(aligned_mfa(list(X1, X2), list(idx1)))
  expect_error(aligned_mfa(list(X1, X2), list(idx1, c(idx2, 1L))))
  expect_error(aligned_mfa(list(X1, X2), list(idx1, rep(999L, nrow(X2))), N = N))
  expect_error(aligned_mfa(list(X1, X2), list(idx1, c(idx2[1], NA_integer_, idx2[-c(1, 2)]))))
  expect_error(aligned_mfa(list(X1, X2), list(idx1, idx2), N = 5))

  fit <- aligned_mfa(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), N = N, ncomp = 2)
  expect_s3_class(fit, "aligned_mfa")
  expect_true(inherits(fit, "multiblock_biprojector"))
  expect_equal(length(fit$block_indices), 2)
  expect_equal(names(fit$block_indices), c("X1", "X2"))
  expect_equal(nrow(multivarious::scores(fit)), N)
  expect_equal(multivarious::ncomp(fit), 2)
  expect_equal(fit$N, N)
})

test_that("aligned_mfa agrees with MFA when all blocks share rows", {
  set.seed(2)
  N <- 40
  X1 <- matrix(rnorm(N * 8), N, 8)
  X2 <- matrix(rnorm(N * 5), N, 5)
  blocks <- list(X1 = X1, X2 = X2)
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  fit_aligned <- aligned_mfa(
    blocks,
    idx,
    ncomp = 2,
    normalization = "MFA",
    ridge = 1e-10,
    max_iter = 200,
    tol = 1e-10
  )
  fit_mfa <- mfa(blocks, ncomp = 2, normalization = "MFA")

  S1 <- multivarious::scores(fit_aligned)
  S2 <- multivarious::scores(fit_mfa)
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))

  rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
  expect_lt(rel, 1e-2)
})

test_that("aligned_mfa handles p > n blocks without error", {
  set.seed(3)
  N <- 25
  X1 <- matrix(rnorm(10 * 60), 10, 60)
  X2 <- matrix(rnorm(12 * 55), 12, 55)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit <- aligned_mfa(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), N = N, ncomp = 3)
  expect_s3_class(fit, "aligned_mfa")
  expect_equal(dim(multivarious::scores(fit)), c(N, 3))
})

test_that("feature_groups='colnames' shrinks within-group loading differences in aligned_mfa", {
  set.seed(4)
  N <- 50
  p_shared <- 6
  p1_unique <- 4
  p2_unique <- 3
  cn_shared <- paste0("f", seq_len(p_shared))
  cn1 <- c(cn_shared, paste0("u1_", seq_len(p1_unique)))
  cn2 <- c(cn_shared, paste0("u2_", seq_len(p2_unique)))

  X1 <- matrix(rnorm(N * length(cn1)), N, length(cn1))
  X2 <- matrix(rnorm(N * length(cn2)), N, length(cn2))
  colnames(X1) <- cn1
  colnames(X2) <- cn2
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  fit0 <- aligned_mfa(list(X1 = X1, X2 = X2), idx, ncomp = 2, feature_groups = "colnames", feature_lambda = 0)
  fit1 <- aligned_mfa(list(X1 = X1, X2 = X2), idx, ncomp = 2, feature_groups = "colnames", feature_lambda = 5)

  V0_1 <- fit0$V_list$X1[cn_shared, , drop = FALSE]
  V0_2 <- fit0$V_list$X2[cn_shared, , drop = FALSE]
  V1_1 <- fit1$V_list$X1[cn_shared, , drop = FALSE]
  V1_2 <- fit1$V_list$X2[cn_shared, , drop = FALSE]

  d0 <- rowSums((V0_1 - V0_2)^2)
  d1 <- rowSums((V1_1 - V1_2)^2)
  expect_lt(mean(d1), mean(d0))
})
