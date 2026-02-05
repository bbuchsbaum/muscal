library(testthat)
library(muscal)

.expand_rowsum_to_N <- function(x, idx, N) {
  rs <- rowsum(x, idx, reorder = FALSE)
  out <- matrix(0, nrow = N, ncol = ncol(x))
  out[as.integer(rownames(rs)), ] <- rs
  out
}

.sim_aligned_mcca <- function(N, k, p_vec, n_vec, noise_sd = 0.05) {
  Z <- matrix(rnorm(N * k), nrow = N, ncol = k)
  idx_list <- lapply(n_vec, function(nk) sample.int(N, nk, replace = TRUE))
  blocks <- Map(function(p, idx) {
    A <- matrix(rnorm(p * k), nrow = p, ncol = k)
    signal <- Z[idx, , drop = FALSE] %*% t(A)
    signal + matrix(rnorm(length(idx) * p, sd = noise_sd), nrow = length(idx), ncol = p)
  }, p_vec, idx_list)
  list(blocks = blocks, idx = idx_list, Z = Z)
}

test_that("aligned_mcca validates row_index and returns expected structure", {
  set.seed(1)
  N <- 30
  X1 <- matrix(rnorm(10 * 6), 10, 6)
  X2 <- matrix(rnorm(12 * 5), 12, 5)
  idx1 <- sample.int(N, nrow(X1), replace = FALSE)
  idx2 <- sample.int(N, nrow(X2), replace = FALSE)

  expect_error(aligned_mcca(list(X1, X2), list(idx1)))
  expect_error(aligned_mcca(list(X1, X2), list(idx1, c(idx2, 1L))))
  expect_error(aligned_mcca(list(X1, X2), list(idx1, rep(999L, nrow(X2))), N = N))
  expect_error(aligned_mcca(list(X1, X2), list(idx1, c(idx2[1], NA_integer_, idx2[-c(1, 2)]))))
  expect_error(aligned_mcca(list(X1, X2), list(idx1, idx2), N = 5))

  fit <- aligned_mcca(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), N = N, ncomp = 2)
  expect_s3_class(fit, "aligned_mcca")
  expect_true(inherits(fit, "multiblock_biprojector"))
  expect_equal(length(fit$block_indices), 2)
  expect_equal(names(fit$block_indices), c("X1", "X2"))
  expect_equal(nrow(multivarious::scores(fit)), N)
  expect_equal(multivarious::ncomp(fit), 2)
  expect_equal(fit$N, N)
  expect_equal(length(fit$partial_scores), 2)
  expect_equal(length(fit$canonical_weights), 2)
})

test_that("aligned_mcca infers N when not supplied", {
  set.seed(1)
  X1 <- matrix(rnorm(10 * 6), 10, 6)
  X2 <- matrix(rnorm(12 * 5), 12, 5)
  idx1 <- c(1:9, 15)
  idx2 <- sample.int(15, nrow(X2), replace = TRUE)

  fit <- aligned_mcca(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
  expect_equal(fit$N, 15)
  expect_equal(nrow(multivarious::scores(fit)), 15)

  idx_bad <- idx2
  idx_bad[1] <- 0L
  expect_error(aligned_mcca(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx_bad)))
})

test_that("aligned_mcca supports per-block preproc lists and validates ncomp", {
  set.seed(2)
  N <- 20
  X1 <- matrix(rnorm(N * 6), N, 6)
  X2 <- matrix(rnorm(N * 5), N, 5)
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  fit <- aligned_mcca(
    list(X1 = X1, X2 = X2),
    idx,
    N = N,
    ncomp = 2,
    preproc = list(multivarious::center(), multivarious::center())
  )
  expect_s3_class(fit, "aligned_mcca")

  expect_error(aligned_mcca(list(X1 = X1, X2 = X2), idx, N = N, ncomp = 1.5))
})

test_that("aligned_mcca agrees with mcca when all blocks share rows", {
  set.seed(2)
  N <- 40
  X1 <- matrix(rnorm(N * 8), N, 8)
  X2 <- matrix(rnorm(N * 5), N, 5)
  blocks <- list(X1 = X1, X2 = X2)
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  fit_aligned <- aligned_mcca(blocks, idx, N = N, ncomp = 2, ridge = 1e-6)
  fit_mcca <- mcca(blocks, ncomp = 2, ridge = 1e-6)

  S1 <- multivarious::scores(fit_aligned)
  S2 <- multivarious::scores(fit_mcca)
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))

  rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
  expect_lt(rel, 1e-10)
})

test_that("aligned_mcca block_weights can drop a block cleanly", {
  set.seed(7)
  N <- 35
  X1 <- matrix(rnorm(N * 15), N, 15)
  X2 <- matrix(rnorm(N * 10), N, 10)
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  res_drop <- aligned_mcca(list(X1 = X1, X2 = X2), idx, N = N, ncomp = 2, block_weights = c(1, 0))
  expect_lt(max(abs(res_drop$partial_scores[[2]])), 1e-10)
  idx2 <- res_drop$block_indices[[2]]
  expect_lt(max(abs(res_drop$v[idx2, , drop = FALSE])), 1e-10)
})

test_that("aligned_mcca uses eigen() fallback when ncomp == N", {
  set.seed(8)
  N <- 6
  X1 <- matrix(rnorm(N * 8), N, 8)
  X2 <- matrix(rnorm(N * 7), N, 7)
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  fit <- aligned_mcca(list(X1 = X1, X2 = X2), idx, N = N, ncomp = N, preproc = multivarious::pass())
  expect_equal(ncol(multivarious::scores(fit)), N)
})

test_that("aligned_mcca passes mathematical sanity checks on synthetic aligned low-rank data", {
  set.seed(10)
  N <- 80
  k <- 3
  Z <- matrix(rnorm(N * k), nrow = N, ncol = k)
  idx <- list(
    B1 = seq_len(N),                       # full coverage of reference rows
    B2 = sample.int(N, N, replace = FALSE), # permutation (row order differs)
    B3 = sample(rep(seq_len(N), 2))         # balanced duplicates (each row twice)
  )

  p_vec <- c(B1 = 25, B2 = 30, B3 = 20)
  blocks <- Map(function(p, idxk) {
    A <- matrix(rnorm(p * k), nrow = p, ncol = k)
    signal <- Z[idxk, , drop = FALSE] %*% t(A)
    signal + matrix(rnorm(length(idxk) * p, sd = 0.02), nrow = length(idxk), ncol = p)
  }, p_vec, idx)

  res <- aligned_mcca(blocks, idx, N = N, ncomp = k, ridge = 1e-6)

  S <- multivarious::scores(res)
  expect_equal(dim(S), c(N, k))

  XtX <- crossprod(S)
  expect_lt(max(abs(XtX - diag(diag(XtX)))), 1e-6)
  expect_equal(unname(diag(XtX)), unname(res$lambda[seq_len(k)]), tolerance = 1e-6)
  expect_equal(unname(res$sdev[seq_len(k)]^2), unname(res$lambda[seq_len(k)]), tolerance = 1e-10)

  # Sum of block partial scores mapped into reference space reproduces compromise scores
  S_sum <- Reduce(`+`, lapply(seq_along(res$partial_scores), function(k) {
    .expand_rowsum_to_N(res$partial_scores[[k]], res$row_index[[k]], res$N)
  }))
  expect_equal(unname(S_sum), unname(S), tolerance = 1e-6)

  # Recover latent space (up to rotation/sign) via canonical correlations
  cc <- cancor(scale(Z, center = TRUE, scale = FALSE),
               scale(S, center = TRUE, scale = FALSE))$cor
  expect_gt(min(cc[seq_len(k)]), 0.98)

  # Internal consistency between v and canonical_weights
  for (k in seq_along(res$block_indices)) {
    idx_feat <- res$block_indices[[k]]
    V_block <- res$v[idx_feat, , drop = FALSE]
    W_raw <- res$canonical_weights[[k]]
    expect_equal(
      V_block %*% diag(res$sdev),
      res$block_weights[k] * W_raw,
      tolerance = 1e-6
    )
  }
})

test_that("aligned_mcca handles p > n blocks without error", {
  set.seed(3)
  N <- 25
  X1 <- matrix(rnorm(10 * 60), 10, 60)
  X2 <- matrix(rnorm(12 * 55), 12, 55)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit <- aligned_mcca(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), N = N, ncomp = 3)
  expect_s3_class(fit, "aligned_mcca")
  expect_equal(dim(multivarious::scores(fit)), c(N, 3))
})

test_that("aligned_mcca ridge=0 warns and still returns a fit for rank-deficient blocks", {
  set.seed(5)
  N <- 20
  X_singular <- matrix(rnorm(N * 5), nrow = N, ncol = 5)   # K is singular (rank <= 5)
  X_fullrank <- matrix(rnorm(N * 30), nrow = N, ncol = 30) # typically full row rank
  idx <- list(B1 = seq_len(N), B2 = seq_len(N))

  expect_warning(
    res <- aligned_mcca(list(B1 = X_singular, B2 = X_fullrank), idx, N = N, ncomp = 2, ridge = 0),
    regexp = "ridge=0 yields a singular system"
  )
  expect_s3_class(res, "aligned_mcca")
})

test_that("anchored_mcca matches aligned_mcca with Y included", {
  set.seed(4)
  N <- 30
  Y <- matrix(rnorm(N * 5), N, 5)
  X1 <- matrix(rnorm(20 * 10), 20, 10)
  X2 <- matrix(rnorm(15 * 8), 15, 8)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit_anchor <- anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
  expect_s3_class(fit_anchor, "anchored_mcca")
  expect_s3_class(fit_anchor, "aligned_mcca")

  fit_aligned <- aligned_mcca(
    X = c(list(Y = Y), list(X1 = X1, X2 = X2)),
    row_index = c(list(Y = seq_len(N)), list(X1 = idx1, X2 = idx2)),
    N = N,
    ncomp = 2
  )

  S1 <- multivarious::scores(fit_anchor)
  S2 <- multivarious::scores(fit_aligned)
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
  expect_lt(rel, 1e-12)
})

test_that("anchored_mcca adapts block_weights for Y and validates names/lengths", {
  set.seed(6)
  N <- 30
  Y <- matrix(rnorm(N * 5), N, 5)
  X1 <- matrix(rnorm(20 * 10), 20, 10)
  X2 <- matrix(rnorm(15 * 8), 15, 8)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  # X-only weights -> prepend Y weight 1
  fit_xonly <- anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                             ncomp = 2, block_weights = c(0.5, 2))
  expect_equal(fit_xonly$block_weights, c(1, 0.5, 2))

  # Named X-only weights -> still prepend Y=1 and order as (Y, X1, X2)
  fit_named <- anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                             ncomp = 2, block_weights = c(X1 = 0.5, X2 = 2))
  expect_equal(fit_named$block_weights, c(1, 0.5, 2))

  # Wrong length
  expect_error(
    anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                  ncomp = 2, block_weights = 1)
  )

  # Missing X block name
  expect_error(
    anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2),
                  ncomp = 2, block_weights = c(X1 = 1))
  )
})

test_that("anchored_mcca assigns default X names and validates row_index values", {
  set.seed(9)
  N <- 25
  Y <- matrix(rnorm(N * 4), N, 4)
  X1 <- matrix(rnorm(10 * 6), 10, 6)
  X2 <- matrix(rnorm(12 * 5), 12, 5)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit <- anchored_mcca(Y, list(X1, X2), list(idx1, idx2), ncomp = 2)
  expect_equal(fit$names, c("Y", "X1", "X2"))

  idx_bad <- idx2
  idx_bad[1] <- NA_integer_
  expect_error(anchored_mcca(Y, list(X1, X2), list(idx1, idx_bad), ncomp = 2))

  idx_oob <- idx2
  idx_oob[1] <- N + 1L
  expect_error(anchored_mcca(Y, list(X1, X2), list(idx1, idx_oob), ncomp = 2))
})

test_that("aligned_mcca use_future path matches serial", {
  skip_if_not_installed("future")
  skip_if_not_installed("furrr")

  op <- future::plan()
  on.exit(future::plan(op), add = TRUE)
  future::plan(future::sequential)

  set.seed(22)
  N <- 40
  X1 <- matrix(rnorm(N * 10), N, 10)
  X2 <- matrix(rnorm(N * 12), N, 12)
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  res_serial <- aligned_mcca(list(X1 = X1, X2 = X2), idx, N = N, ncomp = 2, use_future = FALSE)
  res_future <- aligned_mcca(list(X1 = X1, X2 = X2), idx, N = N, ncomp = 2, use_future = TRUE)

  S1 <- multivarious::scores(res_serial)
  S2 <- multivarious::scores(res_future)
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
  expect_lt(rel, 1e-12)
})
