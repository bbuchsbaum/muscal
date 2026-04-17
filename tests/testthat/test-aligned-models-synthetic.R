library(testthat)
library(muscal)

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
    row_index[[k]] <- idx
  }

  list(Y = Y, X = X, row_index = row_index, S_true = S_true)
}

sim_aligned_mfa <- function(N = 80,
                            p = c(12, 12),
                            K = 3,
                            noise_x = 0.10,
                            reps = 1) {
  S_true <- orthonormal_scores(N, K)
  template <- matrix(rnorm(max(p) * K), max(p), K)

  X <- row_index <- vector("list", length(p))
  names(X) <- names(row_index) <- paste0("X", seq_along(p))
  for (k in seq_along(p)) {
    idx <- rep(seq_len(N), each = reps)
    Vk <- template[seq_len(p[k]), , drop = FALSE]
    X[[k]] <- S_true[idx, , drop = FALSE] %*% t(Vk) +
      matrix(rnorm(length(idx) * p[k], sd = noise_x), length(idx), p[k])
    colnames(X[[k]]) <- paste0("roi", seq_len(p[k]))
    row_index[[k]] <- idx
  }

  list(X = X, row_index = row_index, N = N, S_true = S_true)
}

sim_aligned_mcca <- function(N = 70,
                             p = c(14, 12, 10),
                             q = 8,
                             K = 3,
                             noise = 0.12) {
  S_true <- orthonormal_scores(N, K)
  A_y <- matrix(rnorm(q * K), q, K)
  Y <- S_true %*% t(A_y) + matrix(rnorm(N * q, sd = noise), N, q)

  X <- row_index <- vector("list", length(p))
  names(X) <- names(row_index) <- paste0("B", seq_along(p))
  for (k in seq_along(p)) {
    A_k <- matrix(rnorm(p[k] * K), p[k], K)
    idx <- sample.int(N, size = max(2L, floor(0.8 * N)), replace = TRUE)
    X[[k]] <- S_true[idx, , drop = FALSE] %*% t(A_k) +
      matrix(rnorm(length(idx) * p[k], sd = noise), length(idx), p[k])
    row_index[[k]] <- idx
  }

  list(Y = Y, X = X, row_index = row_index, N = N, S_true = S_true)
}

test_that("anchored_mfa matches row_index by names, not list position", {
  set.seed(201)
  sim <- sim_anchored_mfa()

  fit_ref <- anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-6
  )

  fit_named <- anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index[c("X2", "X1")],
    ncomp = 3,
    normalization = "None",
    ridge = 1e-6
  )

  expect_equal(names(fit_named$row_index), names(sim$row_index))
  expect_lt(subspace_rmse(fit_ref$s, fit_named$s), 1e-6)
})

test_that("aligned_mfa matches row_index by names, not list position", {
  set.seed(202)
  sim <- sim_aligned_mfa()

  fit_ref <- aligned_mfa(
    X = sim$X,
    row_index = sim$row_index,
    N = sim$N,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-6
  )

  fit_named <- aligned_mfa(
    X = sim$X,
    row_index = sim$row_index[c("X2", "X1")],
    N = sim$N,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-6
  )

  expect_equal(names(fit_named$row_index), names(sim$row_index))
  expect_lt(subspace_rmse(fit_ref$s, fit_named$s), 1e-6)
})

test_that("aligned_mcca matches row_index by names, not list position", {
  set.seed(203)
  sim <- sim_aligned_mcca()

  fit_ref <- aligned_mcca(
    X = sim$X,
    row_index = sim$row_index,
    N = sim$N,
    ncomp = 3,
    ridge = 1e-4
  )

  fit_named <- aligned_mcca(
    X = sim$X,
    row_index = sim$row_index[rev(names(sim$row_index))],
    N = sim$N,
    ncomp = 3,
    ridge = 1e-4
  )

  expect_equal(names(fit_named$row_index), names(sim$row_index))
  expect_lt(subspace_rmse(fit_ref$s, fit_named$s), 1e-8)
})

test_that("zero-weight blocks do not influence weighted blocks through feature priors", {
  set.seed(204)
  sim <- sim_anchored_mfa(N = 70, p = c(10, 10), K = 2, noise_x = 0.15)

  X_alt <- sim$X
  X_alt[[2]] <- matrix(rnorm(length(X_alt[[2]])), nrow(X_alt[[2]]), ncol(X_alt[[2]]))
  colnames(X_alt[[2]]) <- colnames(sim$X[[2]])

  alpha <- c(1, 1, 0)
  fit_a <- anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = 2,
    normalization = "custom",
    alpha = alpha,
    feature_groups = "colnames",
    feature_lambda = 25,
    ridge = 1e-6,
    max_iter = 80
  )
  fit_b <- anchored_mfa(
    Y = sim$Y,
    X = X_alt,
    row_index = sim$row_index,
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
  expect_lt(max(abs(fit_a$V_list[[2]])), 1e-8)
  expect_lt(max(abs(fit_b$V_list[[2]])), 1e-8)
})

test_that("aligned_mfa does not force ncomp down to the narrowest block width", {
  set.seed(205)
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

test_that("aligned_mcca rejects all-zero block weights", {
  set.seed(206)
  sim <- sim_aligned_mcca(p = c(8, 8), K = 2)

  expect_error(
    aligned_mcca(
      X = sim$X[1:2],
      row_index = sim$row_index[1:2],
      N = sim$N,
      ncomp = 2,
      ridge = 1e-4,
      block_weights = c(0, 0)
    ),
    regexp = "strictly positive|non-zero|positive",
    ignore.case = TRUE
  )
})

test_that("aligned_mcca recovers a meaningful shared row space on synthetic aligned data", {
  set.seed(207)
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
})
