library(testthat)
library(Matrix)
library(muscal)

sim_gamfa <- function(N = 80,
                      q = 6,
                      p = c(12, 12),
                      K = 3,
                      noise_y = 0.05,
                      noise_x = 0.10,
                      smooth_scores = FALSE,
                      reps = 1,
                      matched_loading_sd = 0.15) {
  stopifnot(length(p) >= 1L)

  if (smooth_scores) {
    tt <- seq(0, 1, length.out = N)
    S_true <- cbind(
      cos(2 * pi * tt),
      sin(2 * pi * tt),
      scale(tt, center = TRUE, scale = FALSE)
    )
    if (K > ncol(S_true)) {
      S_true <- cbind(S_true, matrix(rnorm(N * (K - ncol(S_true))), N, K - ncol(S_true)))
    }
    S_true <- S_true[, seq_len(K), drop = FALSE]
  } else {
    S_true <- matrix(rnorm(N * K), N, K)
  }

  S_true <- qr.Q(qr(S_true), complete = FALSE)
  if (ncol(S_true) > K) S_true <- S_true[, seq_len(K), drop = FALSE]

  B_true <- matrix(rnorm(q * K), q, K)
  V_template <- matrix(rnorm(max(p) * K), max(p), K)
  V_true <- lapply(p, function(pk) {
    V_template[seq_len(pk), , drop = FALSE] +
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

  list(
    Y = Y,
    X = X,
    row_index = row_index,
    S_true = S_true
  )
}

procrustes_align <- function(A, B) {
  sv <- svd(crossprod(A, B))
  R <- sv$u %*% t(sv$v)
  A %*% R
}

subspace_rmse <- function(A, B) {
  sqrt(mean((procrustes_align(A, B) - B)^2))
}

dirichlet_energy <- function(S, L) {
  sum((L %*% S) * S)
}

chain_laplacian <- function(n, weight = 1) {
  A <- sparseMatrix(
    i = c(seq_len(n - 1L), seq(2L, n)),
    j = c(seq(2L, n), seq_len(n - 1L)),
    x = weight,
    dims = c(n, n)
  )
  Diagonal(n = n, x = rowSums(A)) - A
}

block_recons <- function(fit, X, row_index) {
  out <- vector("list", length(X) + 1L)
  names(out) <- c("Y", names(X))
  out[["Y"]] <- fit$s %*% t(fit$B)
  for (k in seq_along(X)) {
    out[[names(X)[k]]] <- fit$s[row_index[[k]], , drop = FALSE] %*% t(fit$V_list[[k]])
  }
  out
}

test_that("graph-off synthetic reconstruction agrees with anchored_mfa", {
  set.seed(1)
  dat <- sim_gamfa(N = 60, p = c(10, 14), K = 3)

  fit_base <- anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-8,
    max_iter = 100,
    tol = 1e-8
  )

  fit_graph <- graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 3,
    normalization = "None",
    graph_lambda = 0,
    score_graph_lambda = 0,
    ridge = 1e-8,
    max_iter = 100,
    tol = 1e-8
  )

  rec0 <- block_recons(fit_base, dat$X, dat$row_index)
  rec1 <- block_recons(fit_graph, dat$X, dat$row_index)

  expect_lt(norm(rec0$Y - rec1$Y, type = "F"), 1e-5)
  for (nm in names(dat$X)) {
    expect_lt(norm(rec0[[nm]] - rec1[[nm]], type = "F"), 1e-5)
  }
})

test_that("synthetic score graph lowers Dirichlet energy and improves smooth-score recovery", {
  set.seed(3)
  dat <- sim_gamfa(
    N = 90,
    p = c(12, 12),
    K = 3,
    smooth_scores = TRUE,
    noise_y = 0.30,
    noise_x = 0.25
  )
  L <- chain_laplacian(nrow(dat$Y))

  fit0 <- graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 3,
    normalization = "None",
    graph_lambda = 0,
    score_graph_lambda = 0,
    ridge = 1e-6
  )

  fitS <- graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 3,
    normalization = "None",
    score_graph = L,
    score_graph_form = "laplacian",
    score_graph_lambda = 1,
    ridge = 1e-6
  )

  expect_lt(dirichlet_energy(fitS$s, L), dirichlet_energy(fit0$s, L))
  expect_lt(subspace_rmse(fitS$s, dat$S_true), subspace_rmse(fit0$s, dat$S_true))
})

test_that("rows sharing the same anchor row have partial scores near the same anchor score", {
  set.seed(4)
  dat <- sim_gamfa(N = 50, p = c(10, 10), K = 2, reps = 3, noise_x = 0.05)

  fit <- graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    graph_lambda = 0,
    score_graph_lambda = 0,
    ridge = 1e-6
  )

  for (k in seq_along(dat$X)) {
    idx <- dat$row_index[[k]]
    err <- fit$partial_scores[[k]] - fit$s[idx, , drop = FALSE]
    expect_lt(mean(rowSums(err^2)), 0.10)
  }
})

test_that("name-misaligned row_index is reordered by name rather than position", {
  set.seed(8)
  dat <- sim_gamfa(N = 40, p = c(6, 7), K = 2)
  idx_bad <- dat$row_index[c("X2", "X1")]

  fit_good <- graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    ridge = 1e-8
  )
  fit_reordered <- graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = idx_bad,
    ncomp = 2,
    normalization = "None",
    ridge = 1e-8
  )

  expect_equal(names(fit_reordered$row_index), names(dat$row_index))
  expect_lt(norm(fit_good$s - fit_reordered$s, type = "F"), 1e-8)
})

test_that("invalid adjacency matrices are rejected clearly", {
  set.seed(9)
  dat <- sim_gamfa(N = 30, p = c(5, 5), K = 2)
  bad_adj <- Diagonal(nrow(dat$Y), x = 0)
  bad_adj[1, 2] <- -1
  bad_adj[2, 1] <- -1

  expect_error(
    graph_anchored_mfa(
      Y = dat$Y,
      X = dat$X,
      row_index = dat$row_index,
      ncomp = 2,
      normalization = "None",
      score_graph = bad_adj,
      score_graph_form = "adjacency",
      score_graph_lambda = 1
    ),
    regexp = "non-negative|adjacency weights",
    ignore.case = TRUE
  )
})

test_that("component count is not forced down by the smallest auxiliary block width", {
  set.seed(7)
  dat <- sim_gamfa(N = 60, p = c(1, 12), K = 3)

  fit <- graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 3,
    normalization = "None",
    ridge = 1e-6
  )

  expect_equal(ncol(fit$s), 3)
})

test_that("custom anchor weight matches equivalent rescaling of Y", {
  set.seed(10)
  dat <- sim_gamfa(N = 50, p = c(8, 9), K = 2, noise_y = 0.08, noise_x = 0.08)
  alpha_y <- 4

  fit_alpha <- graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    preproc = multivarious::pass(),
    normalization = "custom",
    alpha = c(alpha_y, 1, 1),
    ridge = 1e-6,
    max_iter = 80,
    tol = 1e-8
  )

  fit_scaled <- graph_anchored_mfa(
    Y = sqrt(alpha_y) * dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    preproc = multivarious::pass(),
    normalization = "custom",
    alpha = c(1, 1, 1),
    ridge = 1e-6,
    max_iter = 80,
    tol = 1e-8
  )

  expect_lt(norm(fit_alpha$s - fit_scaled$s, type = "F"), 1e-6)
  expect_lt(norm(fit_alpha$V_list$X1 - fit_scaled$V_list$X1, type = "F"), 3e-6)
  expect_lt(norm(fit_alpha$V_list$X2 - fit_scaled$V_list$X2, type = "F"), 3e-6)
})
