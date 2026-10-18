# Synthetic tests for aligned_mfa(), anchored_mfa()/linked_mfa(), and aligned_mcca().
#
# These are intended for a package testthat suite. A few are deliberate
# specification tests: they encode behaviors the implementation should satisfy,
# even if the current code does not yet pass all of them.

library(testthat)

orthonormal_scores <- function(N, K) {
  S <- matrix(rnorm(N * K), N, K)
  S <- qr.Q(qr(S), complete = FALSE)
  if (ncol(S) > K) S <- S[, seq_len(K), drop = FALSE]
  S
}

subspace_rmse <- function(A, B) {
  QA <- qr.Q(qr(A), complete = FALSE)
  QB <- qr.Q(qr(B), complete = FALSE)
  if (ncol(QA) > ncol(QB)) QA <- QA[, seq_len(ncol(QB)), drop = FALSE]
  if (ncol(QB) > ncol(QA)) QB <- QB[, seq_len(ncol(QA)), drop = FALSE]
  P1 <- QA %*% t(QA)
  P2 <- QB %*% t(QB)
  sqrt(mean((P1 - P2)^2))
}

mean_matched_loading_distance <- function(V_list) {
  stopifnot(length(V_list) >= 2L)
  cn <- Reduce(intersect, lapply(V_list, rownames))
  stopifnot(length(cn) > 0)
  pairs <- combn(seq_along(V_list), 2, simplify = FALSE)
  d <- c()
  for (pr in pairs) {
    A <- V_list[[pr[1]]][cn, , drop = FALSE]
    B <- V_list[[pr[2]]][cn, , drop = FALSE]
    d <- c(d, sqrt(rowSums((A - B)^2)))
  }
  mean(d)
}

sim_anchored_mfa <- function(N = 80,
                             q = 6,
                             p = c(12, 12),
                             K = 3,
                             noise_y = 0.05,
                             noise_x = 0.10,
                             reps = 1,
                             matched_loading_sd = 0.10) {
  S_true <- orthonormal_scores(N, K)
  B_true <- matrix(rnorm(q * K), q, K)

  template <- matrix(rnorm(max(p) * K), max(p), K)
  V_true <- lapply(p, function(pk) {
    template[seq_len(pk), , drop = FALSE] +
      matrix(rnorm(pk * K, sd = matched_loading_sd), pk, K)
  })

  Y <- S_true %*% t(B_true) + matrix(rnorm(N * q, sd = noise_y), N, q)
  colnames(Y) <- paste0("y", seq_len(q))

  X <- row_index <- vector("list", length(p))
  names(X) <- names(row_index) <- paste0("X", seq_along(p))

  for (k in seq_along(p)) {
    idx <- rep(seq_len(N), each = reps)
    X[[k]] <- S_true[idx, , drop = FALSE] %*% t(V_true[[k]]) +
      matrix(rnorm(length(idx) * p[k], sd = noise_x), length(idx), p[k])
    colnames(X[[k]]) <- paste0("roi", seq_len(p[k]))
    rownames(V_true[[k]]) <- colnames(X[[k]])
    row_index[[k]] <- idx
  }

  list(Y = Y, X = X, row_index = row_index, S_true = S_true, B_true = B_true, V_true = V_true)
}

sim_aligned_mfa <- function(N = 80,
                            p = c(12, 12),
                            K = 3,
                            noise_x = 0.10,
                            reps = 1,
                            matched_loading_sd = 0.10) {
  S_true <- orthonormal_scores(N, K)

  template <- matrix(rnorm(max(p) * K), max(p), K)
  V_true <- lapply(p, function(pk) {
    template[seq_len(pk), , drop = FALSE] +
      matrix(rnorm(pk * K, sd = matched_loading_sd), pk, K)
  })

  X <- row_index <- vector("list", length(p))
  names(X) <- names(row_index) <- paste0("X", seq_along(p))

  for (k in seq_along(p)) {
    idx <- rep(seq_len(N), each = reps)
    X[[k]] <- S_true[idx, , drop = FALSE] %*% t(V_true[[k]]) +
      matrix(rnorm(length(idx) * p[k], sd = noise_x), length(idx), p[k])
    colnames(X[[k]]) <- paste0("roi", seq_len(p[k]))
    rownames(V_true[[k]]) <- colnames(X[[k]])
    row_index[[k]] <- idx
  }

  list(X = X, row_index = row_index, N = N, S_true = S_true, V_true = V_true)
}

sim_aligned_mcca <- function(N = 70,
                             p = c(14, 12, 10),
                             q = 8,
                             K = 3,
                             noise = 0.12,
                             keep_frac = c(1.0, 0.85, 0.75),
                             duplicate_frac = 0.15) {
  S_true <- orthonormal_scores(N, K)
  A_y <- matrix(rnorm(q * K), q, K)
  Y <- S_true %*% t(A_y) + matrix(rnorm(N * q, sd = noise), N, q)
  colnames(Y) <- paste0("y", seq_len(q))

  X <- row_index <- vector("list", length(p))
  names(X) <- names(row_index) <- paste0("B", seq_along(p))

  for (k in seq_along(p)) {
    A_k <- matrix(rnorm(p[k] * K), p[k], K)
    keep <- sort(sample.int(N, size = max(2L, floor(keep_frac[k] * N)), replace = FALSE))
    idx <- keep
    n_dup <- floor(length(idx) * duplicate_frac)
    if (n_dup > 0L) idx <- c(idx, sample(idx, size = n_dup, replace = TRUE))
    X[[k]] <- S_true[idx, , drop = FALSE] %*% t(A_k) +
      matrix(rnorm(length(idx) * p[k], sd = noise), length(idx), p[k])
    colnames(X[[k]]) <- paste0("v", seq_len(p[k]))
    row_index[[k]] <- idx
  }

  list(Y = Y, X = X, row_index = row_index, N = N, S_true = S_true)
}

test_that("anchored_mfa recovers the shared anchored row space on synthetic data", {
  skip_if_not(exists("anchored_mfa", mode = "function"))
  set.seed(1)
  sim <- sim_anchored_mfa()

  fit <- anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-6,
    max_iter = 100
  )

  expect_equal(dim(fit$s), c(nrow(sim$Y), 3))
  expect_lt(subspace_rmse(fit$s, sim$S_true), 0.20)
  expect_true(all(fit$block_fit$r2 > 0.70))
})

test_that("linked_mfa is a strict alias of anchored_mfa", {
  skip_if_not(exists("anchored_mfa", mode = "function"))
  skip_if_not(exists("linked_mfa", mode = "function"))
  set.seed(2)
  sim <- sim_anchored_mfa()

  fit1 <- anchored_mfa(sim$Y, sim$X, sim$row_index, ncomp = 3, normalization = "None", ridge = 1e-6)
  fit2 <- linked_mfa(sim$Y, sim$X, sim$row_index, ncomp = 3, normalization = "None", ridge = 1e-6)

  expect_equal(fit1$s, fit2$s, tolerance = 1e-10)
  expect_equal(fit1$v, fit2$v, tolerance = 1e-10)
})

test_that("feature grouping shrinks matched loadings in anchored_mfa", {
  skip_if_not(exists("anchored_mfa", mode = "function"))
  set.seed(3)
  sim <- sim_anchored_mfa(noise_x = 0.20, matched_loading_sd = 0.25)

  fit_free <- anchored_mfa(
    sim$Y, sim$X, sim$row_index,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-6,
    feature_groups = NULL,
    feature_lambda = 0,
    max_iter = 80
  )
  fit_grouped <- anchored_mfa(
    sim$Y, sim$X, sim$row_index,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-6,
    feature_groups = "colnames",
    feature_lambda = 25,
    max_iter = 80
  )

  expect_lt(
    mean_matched_loading_distance(fit_grouped$V_list),
    mean_matched_loading_distance(fit_free$V_list)
  )
})

test_that("aligned_mfa recovers the shared row space without a privileged anchor", {
  skip_if_not(exists("aligned_mfa", mode = "function"))
  set.seed(4)
  sim <- sim_aligned_mfa()

  fit <- aligned_mfa(
    X = sim$X,
    row_index = sim$row_index,
    N = sim$N,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-6,
    max_iter = 100
  )

  expect_equal(dim(fit$s), c(sim$N, 3))
  expect_lt(subspace_rmse(fit$s, sim$S_true), 0.20)
  expect_true(all(fit$block_fit$r2 > 0.70))
})

test_that("aligned_mfa should match row_index by names, not by list position [spec]", {
  skip_if_not(exists("aligned_mfa", mode = "function"))
  set.seed(5)
  sim <- sim_aligned_mfa()

  fit_ref <- aligned_mfa(
    X = sim$X,
    row_index = sim$row_index,
    N = sim$N,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-6,
    max_iter = 80
  )

  idx_swapped <- sim$row_index[c("X2", "X1")]
  fit_named <- aligned_mfa(
    X = sim$X,
    row_index = idx_swapped,
    N = sim$N,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-6,
    max_iter = 80
  )

  expect_lt(subspace_rmse(fit_ref$s, fit_named$s), 1e-8)
})

test_that("aligned_mfa should allow ncomp larger than min block feature count when ridge is present [spec]", {
  skip_if_not(exists("aligned_mfa", mode = "function"))
  set.seed(6)
  sim <- sim_aligned_mfa(N = 60, p = c(1, 6), K = 3, noise_x = 0.05)

  fit <- aligned_mfa(
    X = sim$X,
    row_index = sim$row_index,
    N = sim$N,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-4,
    max_iter = 50
  )

  expect_equal(ncol(fit$s), 3L)
})

test_that("zero-weight blocks should not influence weighted blocks through feature priors [spec]", {
  skip_if_not(exists("anchored_mfa", mode = "function"))
  set.seed(7)
  sim <- sim_anchored_mfa(N = 70, p = c(10, 10), K = 2, noise_x = 0.15)

  X_alt <- sim$X
  X_alt[[2]] <- matrix(rnorm(length(X_alt[[2]])), nrow(X_alt[[2]]), ncol(X_alt[[2]]))
  colnames(X_alt[[2]]) <- colnames(sim$X[[2]])

  alpha <- c(1, 1, 0) # Y, X1, X2
  fit_a <- anchored_mfa(
    sim$Y, sim$X, sim$row_index,
    ncomp = 2,
    normalization = "custom",
    alpha = alpha,
    feature_groups = "colnames",
    feature_lambda = 25,
    ridge = 1e-6,
    max_iter = 80
  )
  fit_b <- anchored_mfa(
    sim$Y, X_alt, sim$row_index,
    ncomp = 2,
    normalization = "custom",
    alpha = alpha,
    feature_groups = "colnames",
    feature_lambda = 25,
    ridge = 1e-6,
    max_iter = 80
  )

  expect_lt(subspace_rmse(fit_a$s, fit_b$s), 1e-8)
  expect_equal(fit_a$V_list[[1]], fit_b$V_list[[1]], tolerance = 1e-6)
})

test_that("aligned_mcca recovers a meaningful shared row space on synthetic data", {
  skip_if_not(exists("aligned_mcca", mode = "function"))
  set.seed(8)
  sim <- sim_aligned_mcca()

  fit <- aligned_mcca(
    X = sim$X,
    row_index = sim$row_index,
    N = sim$N,
    ncomp = 3,
    ridge = 1e-4
  )

  expect_equal(dim(fit$s), c(sim$N, 3))
  expect_lt(subspace_rmse(fit$s, sim$S_true), 0.30)

  mapped_cor <- vapply(seq_along(sim$X), function(k) {
    C <- suppressWarnings(cor(fit$partial_scores[[k]], sim$S_true[sim$row_index[[k]], , drop = FALSE]))
    mean(abs(diag(C)))
  }, numeric(1))
  expect_gt(mean(mapped_cor), 0.40)
})

test_that("anchored_mcca matches aligned_mcca with Y inserted as an ordinary aligned block", {
  skip_if_not(exists("aligned_mcca", mode = "function"))
  skip_if_not(exists("anchored_mcca", mode = "function"))
  set.seed(9)
  sim <- sim_aligned_mcca()

  fit_aligned <- aligned_mcca(
    X = c(list(Y = sim$Y), sim$X),
    row_index = c(list(Y = seq_len(sim$N)), sim$row_index),
    N = sim$N,
    ncomp = 3,
    ridge = 1e-4
  )
  fit_anchored <- anchored_mcca(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-4
  )

  expect_lt(subspace_rmse(fit_aligned$s, fit_anchored$s), 1e-8)
})

test_that("aligned_mcca should reject all-zero block weights [spec]", {
  skip_if_not(exists("aligned_mcca", mode = "function"))
  set.seed(10)
  sim <- sim_aligned_mcca(p = c(8, 8), K = 2)

  expect_error(
    aligned_mcca(
      X = sim$X,
      row_index = sim$row_index,
      N = sim$N,
      ncomp = 2,
      ridge = 1e-4,
      normalization = "custom",
      alpha = c(0, 0)
    ),
    regexp = "at least one|positive|non-zero",
    ignore.case = TRUE
  )
})

test_that("anchored_mfa objective trace is non-increasing without feature priors", {
  skip_if_not(exists("anchored_mfa", mode = "function"))
  set.seed(11)
  sim <- sim_anchored_mfa(noise_x = 0.15)

  fit <- anchored_mfa(
    sim$Y, sim$X, sim$row_index,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-6,
    feature_groups = NULL,
    feature_lambda = 0,
    max_iter = 80
  )

  diffs <- diff(fit$objective_trace)
  expect_true(all(diffs <= 1e-8))
})

test_that("aligned_mfa objective trace is non-increasing without feature priors", {
  skip_if_not(exists("aligned_mfa", mode = "function"))
  set.seed(12)
  sim <- sim_aligned_mfa(noise_x = 0.15)

  fit <- aligned_mfa(
    X = sim$X,
    row_index = sim$row_index,
    N = sim$N,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-6,
    feature_groups = NULL,
    feature_lambda = 0,
    max_iter = 80
  )

  diffs <- diff(fit$objective_trace)
  expect_true(all(diffs <= 1e-8))
})
