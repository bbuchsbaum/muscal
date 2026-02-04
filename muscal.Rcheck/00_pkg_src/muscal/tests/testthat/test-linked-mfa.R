library(testthat)
library(muscal)

test_that("linked_mfa validates row_index and returns expected structure", {
  set.seed(1)
  Y <- matrix(rnorm(30 * 4), 30, 4)
  X1 <- matrix(rnorm(10 * 6), 10, 6)
  X2 <- matrix(rnorm(12 * 5), 12, 5)
  idx1 <- sample.int(nrow(Y), nrow(X1), replace = FALSE)
  idx2 <- sample.int(nrow(Y), nrow(X2), replace = FALSE)

  expect_error(linked_mfa(Y, list(X1, X2), list(idx1)))
  expect_error(linked_mfa(Y, list(X1, X2), list(idx1, c(idx2, 1L))))
  expect_error(linked_mfa(Y, list(X1, X2), list(idx1, rep(999L, nrow(X2)))))
  expect_error(linked_mfa(Y, list(X1, X2), list(idx1, c(idx2[1], NA_integer_, idx2[-c(1, 2)]))))

  fit <- linked_mfa(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
  expect_s3_class(fit, "linked_mfa")
  expect_true(inherits(fit, "multiblock_biprojector"))
  expect_equal(length(fit$block_indices), 3)
  expect_equal(names(fit$block_indices), c("Y", "X1", "X2"))
  expect_equal(nrow(multivarious::scores(fit)), nrow(Y))
  expect_equal(multivarious::ncomp(fit), 2)
})

test_that("linked_mfa agrees with MFA when all blocks share rows", {
  set.seed(2)
  N <- 40
  Y <- matrix(rnorm(N * 6), N, 6)
  X1 <- matrix(rnorm(N * 8), N, 8)
  X2 <- matrix(rnorm(N * 5), N, 5)
  blocks <- list(Y = Y, X1 = X1, X2 = X2)

  # Full overlap mapping
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  fit_lmfa <- linked_mfa(
    Y,
    list(X1 = X1, X2 = X2),
    idx,
    ncomp = 2,
    normalization = "MFA",
    ridge = 1e-10,
    max_iter = 200,
    tol = 1e-10
  )
  fit_mfa <- mfa(blocks, ncomp = 2, normalization = "MFA")

  S1 <- multivarious::scores(fit_lmfa)
  S2 <- multivarious::scores(fit_mfa)

  # Compare score subspaces (rotation/sign invariant) via projection matrices.
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))

  rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
  expect_lt(rel, 1e-2)
})

test_that("feature_groups='colnames' shrinks within-group loading differences", {
  set.seed(3)
  N <- 50
  Y <- matrix(rnorm(N * 3), N, 3)

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

  fit0 <- linked_mfa(Y, list(X1 = X1, X2 = X2), idx, ncomp = 2, feature_groups = "colnames", feature_lambda = 0)
  fit1 <- linked_mfa(Y, list(X1 = X1, X2 = X2), idx, ncomp = 2, feature_groups = "colnames", feature_lambda = 5)

  V0_1 <- fit0$V_list$X1[cn_shared, , drop = FALSE]
  V0_2 <- fit0$V_list$X2[cn_shared, , drop = FALSE]
  V1_1 <- fit1$V_list$X1[cn_shared, , drop = FALSE]
  V1_2 <- fit1$V_list$X2[cn_shared, , drop = FALSE]

  d0 <- rowSums((V0_1 - V0_2)^2)
  d1 <- rowSums((V1_1 - V1_2)^2)

  expect_lt(mean(d1), mean(d0))
})
