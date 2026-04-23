library(testthat)
library(muscal)

.projected_grad_norm <- function(S, G) {
  PG <- G - S %*% ((crossprod(S, G) + t(crossprod(S, G))) / 2)
  sqrt(sum(PG^2))
}

sim_cgamfa <- function(N = 60,
                       q = 5,
                       p = c(10, 10),
                       K = 2,
                       smooth_scores = FALSE,
                       block_shift_sd = 0.35,
                       matched_loading_sd = 0.15,
                       noise_y = 0.05,
                       noise_x = 0.08) {
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

  X <- row_index <- Z_true <- vector("list", length(p))
  names(X) <- names(row_index) <- names(Z_true) <- paste0("X", seq_along(p))

  for (k in seq_along(p)) {
    idx <- seq_len(N)
    block_shift <- matrix(rnorm(N * K, sd = block_shift_sd), N, K)
    Zk <- S_true + block_shift
    Xk <- Zk %*% t(V_true[[k]]) + matrix(rnorm(N * p[k], sd = noise_x), N, p[k])

    colnames(Xk) <- paste0("roi", seq_len(p[k]))
    X[[k]] <- Xk
    row_index[[k]] <- idx
    Z_true[[k]] <- Zk
  }

  list(
    Y = Y,
    X = X,
    row_index = row_index,
    S_true = S_true,
    Z_true = Z_true,
    V_true = V_true
  )
}

mean_coupling_gap <- function(fit) {
  mean(vapply(seq_along(fit$Z_list), function(k) {
    idx <- fit$row_index[[k]]
    mean(rowSums((fit$Z_list[[k]] - fit$s[idx, , drop = FALSE])^2))
  }, numeric(1)))
}

block_sse <- function(fit, X) {
  sum(vapply(seq_along(X), function(k) {
    Zk <- fit$Z_list[[k]]
    Vk <- fit$V_list[[k]]
    sum((X[[k]] - Zk %*% t(Vk))^2)
  }, numeric(1)))
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
  A <- Matrix::sparseMatrix(
    i = c(seq_len(n - 1L), seq(2L, n)),
    j = c(seq(2L, n), seq_len(n - 1L)),
    x = weight,
    dims = c(n, n)
  )
  Matrix::Diagonal(n = n, x = Matrix::rowSums(A)) - A
}

matched_loading_gap <- function(fit, feature_ids = c(1, 2)) {
  blocks <- fit$V_list[feature_ids]
  p <- min(vapply(blocks, nrow, integer(1)))
  mean(vapply(seq_len(p), function(j) {
    sum((blocks[[1]][j, ] - blocks[[2]][j, ])^2)
  }, numeric(1)))
}

test_that("coupled_graph_anchored_mfa returns the expected fitted structure", {
  set.seed(100)
  dat <- sim_cgamfa()

  fit <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    coupling_lambda = 1,
    ridge = 1e-6
  )

  expect_s3_class(fit, "coupled_graph_anchored_mfa")
  expect_s3_class(fit, "linked_mfa")
  expect_length(fit$Z_list, length(dat$X))
  expect_equal(names(fit$Z_list), names(dat$X))
  expect_equal(dim(fit$s), c(nrow(dat$Y), 2))
  expect_equal(dim(fit$Z_list[[1]]), c(nrow(dat$X[[1]]), 2))
  expect_equal(fit$score_representation, "anchor_scores")
  expect_equal(fit$S, fit$s)
  expect_equal(names(fit$score_index), names(dat$X))
  expect_identical(fit$partial_scores, fit$Z_list)
  expect_true(all(is.finite(fit$s)))
  expect_true(all(vapply(fit$Z_list, function(z) all(is.finite(z)), logical(1))))
  expect_true(all(is.finite(fit$objective_trace)))
})

test_that("larger coupling_lambda pulls block scores toward shared anchor scores", {
  set.seed(101)
  dat <- sim_cgamfa(block_shift_sd = 0.45)

  fit_lo <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    coupling_lambda = 0.05,
    ridge = 1e-6,
    max_iter = 80,
    tol = 1e-8
  )

  fit_hi <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    coupling_lambda = 20,
    ridge = 1e-6,
    max_iter = 80,
    tol = 1e-8
  )

  expect_lt(mean_coupling_gap(fit_hi), mean_coupling_gap(fit_lo))
})

test_that("on noiseless misspecified-shared data, the coupled model dramatically outperforms the shared model", {
  set.seed(104)
  dat <- sim_cgamfa(
    N = 40,
    q = 4,
    p = c(8, 8),
    K = 2,
    block_shift_sd = 0.25,
    matched_loading_sd = 0.05,
    noise_y = 0,
    noise_x = 0
  )

  fit_shared <- graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    preproc = multivarious::pass(),
    normalization = "None",
    graph_lambda = 0,
    score_graph_lambda = 0,
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  fit_coupled <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    preproc = multivarious::pass(),
    normalization = "None",
    coupling_lambda = 1,
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  sse_shared <- sum(vapply(seq_along(dat$X), function(k) {
    Sk <- fit_shared$s[dat$row_index[[k]], , drop = FALSE]
    Vk <- fit_shared$V_list[[k]]
    sum((dat$X[[k]] - Sk %*% t(Vk))^2)
  }, numeric(1)))
  sse_coupled <- block_sse(fit_coupled, dat$X)

  expect_lt(sse_coupled, 0.1 * sse_shared)
})

test_that("coupled model improves block reconstruction when blocks deviate from the shared stimulus scores", {
  set.seed(102)
  dat <- sim_cgamfa(block_shift_sd = 0.40, noise_x = 0.06)

  fit_shared <- graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    graph_lambda = 0,
    score_graph_lambda = 0,
    ridge = 1e-6,
    max_iter = 100,
    tol = 1e-8
  )

  fit_coupled <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    coupling_lambda = 1,
    ridge = 1e-6,
    max_iter = 100,
    tol = 1e-8
  )

  sse_shared <- sum(vapply(seq_along(dat$X), function(k) {
    Sk <- fit_shared$s[dat$row_index[[k]], , drop = FALSE]
    Vk <- fit_shared$V_list[[k]]
    sum((dat$X[[k]] - Sk %*% t(Vk))^2)
  }, numeric(1)))
  sse_coupled <- block_sse(fit_coupled, dat$X)

  expect_lt(sse_coupled, 0.9 * sse_shared)
})

test_that("when there are no block-specific deviations, coupled and shared models agree closely", {
  set.seed(105)
  dat <- sim_cgamfa(block_shift_sd = 0, noise_x = 0.05, noise_y = 0.05)

  fit_shared <- graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    graph_lambda = 0,
    score_graph_lambda = 0,
    ridge = 1e-6,
    max_iter = 100,
    tol = 1e-8
  )

  fit_coupled <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    coupling_lambda = 10,
    ridge = 1e-6,
    max_iter = 100,
    tol = 1e-8
  )

  sse_shared <- sum(vapply(seq_along(dat$X), function(k) {
    Sk <- fit_shared$s[dat$row_index[[k]], , drop = FALSE]
    Vk <- fit_shared$V_list[[k]]
    sum((dat$X[[k]] - Sk %*% t(Vk))^2)
  }, numeric(1)))
  sse_coupled <- block_sse(fit_coupled, dat$X)

  expect_lt(abs(sse_coupled - sse_shared) / max(sse_shared, 1e-8), 0.10)
  expect_lt(mean_coupling_gap(fit_coupled), 1e-3)
})

test_that("correct score graph improves smooth shared-score recovery in the coupled model", {
  set.seed(106)
  dat <- sim_cgamfa(
    N = 80,
    p = c(10, 10),
    K = 3,
    smooth_scores = TRUE,
    block_shift_sd = 0.15,
    noise_y = 0.20,
    noise_x = 0.18
  )
  L <- chain_laplacian(nrow(dat$Y))

  fit0 <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 3,
    normalization = "None",
    coupling_lambda = 1,
    score_graph_lambda = 0,
    ridge = 1e-6
  )

  fitS <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 3,
    normalization = "None",
    coupling_lambda = 1,
    score_graph = L,
    score_graph_form = "laplacian",
    score_graph_lambda = 1,
    ridge = 1e-6
  )

  expect_lt(dirichlet_energy(fitS$s, L), dirichlet_energy(fit0$s, L))
  expect_lt(subspace_rmse(fitS$s, dat$S_true), subspace_rmse(fit0$s, dat$S_true))
})

test_that("correct feature graph aligns matched loadings in the coupled model", {
  set.seed(107)
  dat <- sim_cgamfa(
    N = 70,
    p = c(8, 8),
    K = 2,
    block_shift_sd = 0.25,
    matched_loading_sd = 0.02,
    noise_y = 0.08,
    noise_x = 0.10
  )

  edges <- do.call(rbind, lapply(seq_len(ncol(dat$X[[1]])), function(j) {
    data.frame(
      block1 = "X1",
      feature1 = j,
      block2 = "X2",
      feature2 = j,
      weight = 1
    )
  }))

  fit0 <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    coupling_lambda = 1,
    graph_lambda = 0,
    ridge = 1e-6
  )

  fitG <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    coupling_lambda = 1,
    feature_graph = edges,
    graph_lambda = 2,
    ridge = 1e-6
  )

  expect_lt(matched_loading_gap(fitG), matched_loading_gap(fit0))
})

test_that("coupled fit is invariant to block-list order after name alignment", {
  set.seed(108)
  dat <- sim_cgamfa()

  fit_a <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    coupling_lambda = 1,
    ridge = 1e-6
  )

  fit_b <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X[c("X2", "X1")],
    row_index = dat$row_index[c("X2", "X1")],
    ncomp = 2,
    normalization = "None",
    coupling_lambda = 1,
    ridge = 1e-6
  )

  expect_equal(names(fit_b$V_list), c("X2", "X1"))
  expect_lt(norm(fit_a$s - fit_b$s, type = "F"), 1e-6)
  expect_lt(norm(fit_a$Z_list$X1 - fit_b$Z_list$X1, type = "F"), 1e-6)
  expect_lt(norm(fit_a$Z_list$X2 - fit_b$Z_list$X2, type = "F"), 1e-6)
})

test_that("coupled objective trace is non-increasing up to numerical tolerance", {
  set.seed(109)
  dat <- sim_cgamfa(block_shift_sd = 0.30, noise_x = 0.10, noise_y = 0.10)

  fit <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    coupling_lambda = 1,
    score_graph = "knn",
    score_graph_k = 6,
    score_graph_lambda = 0.5,
    ridge = 1e-6,
    max_iter = 60,
    tol = 1e-10
  )

  diffs <- diff(fit$objective_trace)
  expect_lt(tail(fit$objective_trace, 1), fit$objective_trace[[1]])
  expect_lte(max(diffs), 1e-2)
  expect_gte(mean(diffs <= 1e-6), 0.75)
})

test_that("predict.coupled_graph_anchored_mfa returns scores, reconstructions, and responses", {
  set.seed(103)
  dat <- sim_cgamfa()

  fit <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    coupling_lambda = 1,
    ridge = 1e-6
  )

  new_rows <- dat$X[[1]][1:5, , drop = FALSE]
  scores <- predict(fit, new_rows, block = "X1", type = "scores")
  recon <- predict(fit, new_rows, block = "X1", type = "reconstruction")
  resp <- predict(fit, new_rows, block = "X1", type = "response")

  expect_equal(dim(scores), c(5, 2))
  expect_equal(dim(recon), c(5, ncol(dat$X[[1]])))
  expect_equal(dim(resp), c(5, ncol(dat$Y)))
})

test_that("coupled_graph_anchored_mfa supports orthonormal score constraint", {
  set.seed(110)
  dat <- sim_cgamfa()

  fit <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    score_constraint = "orthonormal",
    coupling_lambda = 1,
    max_iter = 15
  )

  expect_equal(crossprod(fit$s), diag(2), tolerance = 1e-5)
  expect_true(all(is.finite(fit$objective_trace)))
})

test_that("coupled_graph_anchored_mfa orthonormal path has a small projected gradient residual", {
  set.seed(111)
  dat <- sim_cgamfa(smooth_scores = TRUE)

  fit <- coupled_graph_anchored_mfa(
    Y = dat$Y,
    X = dat$X,
    row_index = dat$row_index,
    ncomp = 2,
    normalization = "None",
    preproc = multivarious::pass(),
    score_constraint = "orthonormal",
    score_graph = chain_laplacian(nrow(dat$Y)),
    score_graph_form = "laplacian",
    score_graph_lambda = 0.5,
    coupling_lambda = 1,
    ridge = 1e-6,
    max_iter = 40,
    tol = 1e-8
  )

  grad <- muscal:::.cgamfa_score_gradient(
    S = fit$s,
    Y = dat$Y,
    B = fit$B,
    Z_list = fit$Z_list,
    row_index = fit$row_index,
    alpha_y = fit$alpha_blocks[[1]],
    coupling_lambda = fit$coupling_lambda,
    score_graph_laplacian = fit$score_graph_laplacian,
    score_graph_lambda = fit$score_graph_lambda,
    ridge = fit$ridge
  )
  pg_rel <- .projected_grad_norm(fit$s, grad) / (sqrt(sum(grad^2)) + 1e-12)

  expect_lt(pg_rel, 2e-1)
  expect_lte(max(diff(fit$objective_trace)), 1e-6)
})
