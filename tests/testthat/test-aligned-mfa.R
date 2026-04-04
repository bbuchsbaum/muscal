library(testthat)
library(muscal)

.permute_reference_labels <- function(idx, perm) {
  perm[idx]
}

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

test_that("aligned_mfa exposes a standard out-of-sample contract", {
  set.seed(11)
  N <- 45
  X1 <- matrix(rnorm(18 * 6), 18, 6)
  X2 <- matrix(rnorm(20 * 5), 20, 5)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit <- aligned_mfa(
    list(X1 = X1, X2 = X2),
    list(X1 = idx1, X2 = idx2),
    N = N,
    ncomp = 2
  )

  expect_equal(fit$task, "row_alignment")
  expect_equal(fit$fit_spec$method, "aligned_mfa")
  expect_true(fit$fit_spec$refit_supported)
  expect_setequal(fit$oos_types, c("scores", "reconstruction"))
  expect_true(all(c("X1", "X2") %in% names(fit$block_preproc)))

  new_rows <- X2[1:4, , drop = FALSE] + matrix(rnorm(20, sd = 0.05), nrow = 4)
  scores_new <- project(fit, new_rows, block = "X2")
  xhat <- predict(fit, new_rows, block = "X2", type = "reconstruction")

  expect_equal(scores_new, predict(fit, new_rows, block = "X2", type = "scores"))
  expect_equal(dim(scores_new), c(nrow(new_rows), multivarious::ncomp(fit)))
  expect_equal(dim(xhat), dim(new_rows))
})

test_that("aligned_mfa objective trace is finite, non-increasing, and invariant to reference relabeling", {
  set.seed(45)
  N <- 55
  k <- 2
  S <- scale(matrix(rnorm(N * k), N, k), center = TRUE, scale = FALSE)
  idx <- list(
    X1 = sample.int(N, 40, replace = TRUE),
    X2 = sample.int(N, 38, replace = TRUE)
  )
  V1 <- matrix(rnorm(7 * k), 7, k)
  V2 <- matrix(rnorm(6 * k), 6, k)
  X1 <- S[idx$X1, , drop = FALSE] %*% t(V1) + matrix(rnorm(length(idx$X1) * 7, sd = 0.02), length(idx$X1), 7)
  X2 <- S[idx$X2, , drop = FALSE] %*% t(V2) + matrix(rnorm(length(idx$X2) * 6, sd = 0.02), length(idx$X2), 6)

  fit <- aligned_mfa(
    list(X1 = X1, X2 = X2),
    idx,
    N = N,
    ncomp = k,
    preproc = multivarious::pass(),
    normalization = "None",
    max_iter = 60,
    tol = 0,
    ridge = 1e-10
  )

  obj <- fit$objective_trace
  expect_gt(length(obj), 0)
  expect_true(all(is.finite(obj)))
  expect_lte(max(diff(obj)), 1e-8)

  perm <- sample.int(N)
  idx_perm <- lapply(idx, .permute_reference_labels, perm = perm)
  fit_perm <- aligned_mfa(
    list(X1 = X1, X2 = X2),
    idx_perm,
    N = N,
    ncomp = k,
    preproc = multivarious::pass(),
    normalization = "None",
    max_iter = 60,
    tol = 0,
    ridge = 1e-10
  )

  S1 <- multivarious::scores(fit)
  S2 <- multivarious::scores(fit_perm)[perm, , drop = FALSE]
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  rel <- norm(P1 - P2, type = "F") / (norm(P1, type = "F") + 1e-12)
  expect_lt(rel, 1e-8)
})
