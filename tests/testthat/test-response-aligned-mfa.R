library(testthat)
library(muscal)

.sim_response_aligned_mfa <- function(n = c(50, 44),
                                      p = c(8, 7),
                                      q = 4,
                                      K = 2,
                                      noise_x = 0.04,
                                      noise_y = 0.04,
                                      seed = 1) {
  set.seed(seed)
  B <- matrix(rnorm(q * K), q, K)
  V_list <- lapply(p, function(pk) matrix(rnorm(pk * K), pk, K))
  Z_list <- lapply(n, function(nk) scale(matrix(rnorm(nk * K), nk, K), center = TRUE, scale = FALSE))

  X <- lapply(seq_along(n), function(k) {
    Z_list[[k]] %*% t(V_list[[k]]) + matrix(rnorm(n[[k]] * p[[k]], sd = noise_x), n[[k]], p[[k]])
  })
  Y <- lapply(seq_along(n), function(k) {
    Z_list[[k]] %*% t(B) + matrix(rnorm(n[[k]] * q, sd = noise_y), n[[k]], q)
  })

  names(X) <- paste0("X", seq_along(X))
  names(Y) <- names(X)
  names(V_list) <- names(X)
  names(Z_list) <- names(X)

  list(X = X, Y = Y, V_list = V_list, B = B, Z_list = Z_list)
}

test_that("response_aligned_mfa validates inputs and exposes the standard contract", {
  set.seed(1)
  X1 <- matrix(rnorm(30), nrow = 10)
  X2 <- matrix(rnorm(40), nrow = 10)
  Y1 <- matrix(rnorm(20), nrow = 10)
  Y2 <- matrix(rnorm(20), nrow = 10)

  expect_error(response_aligned_mfa(Y = list(Y1), X = list(X1, X2)))
  expect_error(response_aligned_mfa(Y = list(Y1 = Y1, X2 = Y2), X = list(X1 = X1, X2 = X2[1:9, , drop = FALSE])))
  expect_error(response_aligned_mfa(Y = list(X1 = Y1, X2 = matrix(rnorm(30), 10, 3)), X = list(X1 = X1, X2 = X2)))
  expect_error(response_aligned_mfa(
    Y = list(X1 = Y1, X2 = Y2),
    X = list(X1 = X1, X2 = X2),
    response_preproc = list(multivarious::center(), multivarious::center())
  ))
  expect_error(response_aligned_mfa(
    Y = list(X1 = Y1, X2 = Y2),
    X = list(X1 = X1, X2 = X2),
    anchor_response = matrix(rnorm(12), 6, 2)
  ), "anchor_map")
  colnames(X1) <- paste0("f", seq_len(ncol(X1)))
  colnames(X2) <- paste0("f", seq_len(ncol(X2)))
  expect_error(response_aligned_mfa(
    Y = list(X1 = Y1, X2 = Y2),
    X = list(X1 = X1, X2 = X2),
    feature_groups = "colnames",
    feature_graph = "colnames",
    feature_lambda = 1,
    graph_lambda = 1
  ), "cannot be combined")

  fit <- response_aligned_mfa(
    Y = list(X1 = Y1, X2 = Y2),
    X = list(X1 = X1, X2 = X2),
    ncomp = 2
  )

  expect_s3_class(fit, "response_aligned_mfa")
  expect_true(inherits(fit, "multiblock_biprojector"))
  expect_equal(fit$task, "response_prediction")
  expect_equal(fit$fit_spec$method, "response_aligned_mfa")
  expect_true(fit$fit_spec$refit_supported)
  expect_setequal(fit$oos_types, c("response", "scores", "reconstruction"))
  expect_true(all(c("X1", "X2") %in% names(fit$block_preproc)))
  expect_true(inherits(fit$response_preproc, "pre_processor"))
  expect_equal(names(fit$score_index), c("X1", "X2"))
  expect_equal(fit$score_representation, "stacked_block_scores")
  expect_equal(multivarious::scores(fit), do.call(rbind, fit$Z_list))
})

test_that("response_aligned_mfa recovers a shared response-informed latent space", {
  sim <- .sim_response_aligned_mfa(seed = 2)

  fit <- response_aligned_mfa(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  pooled_gram <- Reduce(
    "+",
    lapply(seq_along(fit$Z_list), function(k) {
      fit$score_weights[[k]] * crossprod(fit$Z_list[[k]])
    })
  )
  expect_equal(pooled_gram, diag(2), tolerance = 1e-5)

  mse_y <- mean(unlist(lapply(seq_along(sim$Y), function(k) {
    Y_hat <- fit$Z_list[[k]] %*% t(fit$B)
    (sim$Y[[k]] - Y_hat)^2
  })))
  mse_x <- mean(unlist(lapply(seq_along(sim$X), function(k) {
    X_hat <- fit$Z_list[[k]] %*% t(fit$V_list[[k]])
    (sim$X[[k]] - X_hat)^2
  })))

  expect_lt(mse_y, 0.05)
  expect_lt(mse_x, 0.05)
  expect_true(all(is.finite(fit$objective_trace)))
})

test_that("response_aligned_mfa drops anchor machinery when coupling is zero", {
  sim <- .sim_response_aligned_mfa(n = c(32, 30), p = c(6, 5), q = 3, K = 2, seed = 22)
  anchor_map <- list(
    X1 = sample.int(9, nrow(sim$X$X1), replace = TRUE),
    X2 = sample.int(9, nrow(sim$X$X2), replace = TRUE)
  )

  fit_base <- response_aligned_mfa(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  fit_zero <- response_aligned_mfa(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    anchor_map = anchor_map,
    coupling_lambda = 0,
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  expect_null(fit_zero$S)
  expect_null(fit_zero$anchor_map)
  expect_null(fit_zero$anchor_weight)
  expect_equal(fit_zero$B, fit_base$B, tolerance = 1e-8)
  expect_equal(fit_zero$V_list$X1, fit_base$V_list$X1, tolerance = 1e-8)
  expect_equal(fit_zero$V_list$X2, fit_base$V_list$X2, tolerance = 1e-8)
  expect_equal(fit_zero$Z_list$X1, fit_base$Z_list$X1, tolerance = 1e-8)
  expect_equal(fit_zero$Z_list$X2, fit_base$Z_list$X2, tolerance = 1e-8)
})

test_that("response_aligned_mfa supports graph-based feature coupling", {
  set.seed(12)
  K <- 2
  q <- 3
  Z1 <- scale(matrix(rnorm(24 * K), 24, K), center = TRUE, scale = FALSE)
  Z2 <- scale(matrix(rnorm(80 * K), 80, K), center = TRUE, scale = FALSE)
  B <- matrix(rnorm(q * K), q, K)
  V_shared <- matrix(
    c(1.0, 0.0,
      0.7, 0.3,
      -0.5, 0.9,
      0.4, -0.8),
    nrow = 4,
    byrow = TRUE
  )
  V1 <- rbind(
    V_shared,
    matrix(c(0.6, 0.2, -0.4, 0.5), nrow = 2, byrow = TRUE)
  )
  V2 <- rbind(
    V_shared,
    matrix(c(-0.3, 0.6, 0.8, -0.1), nrow = 2, byrow = TRUE)
  )

  X1 <- Z1 %*% t(V1) + matrix(rnorm(nrow(Z1) * nrow(V1), sd = 0.35), nrow(Z1), nrow(V1))
  X2 <- Z2 %*% t(V2) + matrix(rnorm(nrow(Z2) * nrow(V2), sd = 0.05), nrow(Z2), nrow(V2))
  Y1 <- Z1 %*% t(B) + matrix(rnorm(nrow(Z1) * q, sd = 0.04), nrow(Z1), q)
  Y2 <- Z2 %*% t(B) + matrix(rnorm(nrow(Z2) * q, sd = 0.04), nrow(Z2), q)
  colnames(X1) <- c(paste0("f", 1:4), "x1_u1", "x1_u2")
  colnames(X2) <- c(paste0("f", 1:4), "x2_u1", "x2_u2")

  fit_none <- response_aligned_mfa(
    Y = list(X1 = Y1, X2 = Y2),
    X = list(X1 = X1, X2 = X2),
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  fit_graph <- response_aligned_mfa(
    Y = list(X1 = Y1, X2 = Y2),
    X = list(X1 = X1, X2 = X2),
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    feature_graph = "colnames",
    graph_lambda = 8,
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  shared <- paste0("f", 1:4)
  gap_none <- mean(rowSums((fit_none$V_list$X1[shared, , drop = FALSE] - fit_none$V_list$X2[shared, , drop = FALSE])^2))
  gap_graph <- mean(rowSums((fit_graph$V_list$X1[shared, , drop = FALSE] - fit_graph$V_list$X2[shared, , drop = FALSE])^2))

  expect_lt(gap_graph, gap_none)
  expect_true(Matrix::isDiagonal(fit_graph$graph_laplacian) || inherits(fit_graph$graph_laplacian, "Matrix"))
  expect_equal(fit_graph$graph_lambda, 8)
})

test_that("response_aligned_mfa predict() supports X-only and response-refined projection", {
  sim <- .sim_response_aligned_mfa(seed = 3)

  fit <- response_aligned_mfa(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  set.seed(33)
  Z_new <- scale(matrix(rnorm(12 * 2), 12, 2), center = TRUE, scale = FALSE)
  X_new <- Z_new %*% t(sim$V_list$X1) + matrix(rnorm(12 * nrow(sim$V_list$X1), sd = 0.04), 12, nrow(sim$V_list$X1))
  Y_new <- Z_new %*% t(sim$B) + matrix(rnorm(12 * nrow(sim$B), sd = 0.04), 12, nrow(sim$B))

  score_x <- project(fit, X_new, block = "X1")
  resp_x <- predict(fit, X_new, block = "X1", type = "response")
  xhat <- predict(fit, X_new, block = "X1", type = "reconstruction")

  expect_equal(dim(score_x), c(nrow(X_new), multivarious::ncomp(fit)))
  expect_equal(dim(resp_x), dim(Y_new))
  expect_equal(dim(xhat), dim(X_new))
  expect_equal(score_x, predict(fit, X_new, block = "X1", type = "scores"))
  expect_error(
    project(fit, X_new, block = "X1", new_response = Y_new),
    "conditional = TRUE"
  )
  expect_error(
    predict(fit, X_new, block = "X1", type = "response", new_response = Y_new),
    "conditional = TRUE"
  )

  scores_refined <- project(
    fit,
    X_new,
    block = "X1",
    new_response = Y_new,
    conditional = TRUE
  )
  resp_refined <- predict(
    fit,
    X_new,
    block = "X1",
    type = "response",
    new_response = Y_new,
    conditional = TRUE
  )

  expect_equal(dim(scores_refined), dim(score_x))
  expect_equal(dim(resp_refined), dim(Y_new))

  mse_x <- mean((resp_x - Y_new)^2)
  mse_refined <- mean((resp_refined - Y_new)^2)
  expect_lte(mse_refined, mse_x + 1e-8)
  expect_error(project(fit, X_new, block = "missing"), "Unknown block")
})

test_that("response_aligned_mfa refit metadata can rebuild the fit", {
  sim <- .sim_response_aligned_mfa(seed = 4)
  weights <- list(
    X1 = seq(0.5, 1.5, length.out = nrow(sim$Y$X1)),
    X2 = seq(1.0, 0.6, length.out = nrow(sim$Y$X2))
  )

  fit <- response_aligned_mfa(
    Y = sim$Y,
    X = sim$X,
    response_weights = weights,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    max_iter = 40
  )

  refit <- fit$fit_spec$refit$fit_fn(fit$fit_spec$refit$data)

  expect_s3_class(refit, "response_aligned_mfa")
  expect_equal(names(refit$V_list), names(fit$V_list))
  expect_equal(vapply(refit$response_weights, length, integer(1)),
               vapply(fit$response_weights, length, integer(1)))
  expect_true(all(is.finite(refit$sdev)))
})

test_that("response_aligned_mfa stores parsed row weights and uses a deterministic sign convention", {
  sim <- .sim_response_aligned_mfa(n = c(26, 24), p = c(6, 5), q = 3, K = 2, seed = 14)
  N_anchor <- 11
  idx1 <- sample(seq_len(N_anchor), nrow(sim$X$X1), replace = TRUE)
  idx2 <- sample(seq_len(N_anchor), nrow(sim$X$X2), replace = TRUE)

  fit <- response_aligned_mfa(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    response_weights = c(X1 = 0.5, X2 = 1.2),
    anchor_map = list(X1 = idx1, X2 = idx2),
    anchor_weight = 0.75,
    coupling_lambda = 3,
    ridge = 1e-8,
    max_iter = 80,
    tol = 1e-10
  )

  refit_data <- fit$fit_spec$refit$data
  expect_true(is.list(refit_data$response_weights))
  expect_true(is.list(refit_data$anchor_weight))
  expect_equal(vapply(refit_data$response_weights, length, integer(1)), vapply(sim$Y, nrow, integer(1)))
  expect_equal(vapply(refit_data$anchor_weight, length, integer(1)), vapply(sim$X, nrow, integer(1)))

  dominant <- apply(abs(fit$B), 2, which.max)
  expect_true(all(fit$B[cbind(dominant, seq_len(ncol(fit$B)))] >= 0))

  permuted <- fit$fit_spec$refit$permutation_fn(refit_data)
  for (k in names(sim$Y)) {
    key_ref <- apply(signif(refit_data$Y[[k]], 14), 1, paste, collapse = "|")
    key_perm <- apply(signif(permuted$Y[[k]], 14), 1, paste, collapse = "|")
    idx <- match(key_perm, key_ref)
    expect_false(anyNA(idx))
    expect_equal(permuted$response_weights[[k]], refit_data$response_weights[[k]][idx])
  }
})

test_that("response_aligned_mfa can borrow supervision from anchor-level responses", {
  set.seed(41)
  K <- 2
  q <- 3
  N_anchor <- 12
  S_true <- scale(matrix(rnorm(N_anchor * K), N_anchor, K), center = TRUE, scale = FALSE)
  B_true <- matrix(rnorm(q * K), q, K)
  V1 <- matrix(rnorm(6 * K), 6, K)
  V2 <- matrix(rnorm(5 * K), 5, K)

  idx1 <- sample(seq_len(N_anchor), 38, replace = TRUE)
  idx2 <- sample(seq_len(N_anchor), 34, replace = TRUE)
  Z1 <- S_true[idx1, , drop = FALSE] + matrix(rnorm(length(idx1) * K, sd = 0.04), length(idx1), K)
  Z2 <- S_true[idx2, , drop = FALSE] + matrix(rnorm(length(idx2) * K, sd = 0.04), length(idx2), K)

  X <- list(
    X1 = Z1 %*% t(V1) + matrix(rnorm(length(idx1) * nrow(V1), sd = 0.06), length(idx1), nrow(V1)),
    X2 = Z2 %*% t(V2) + matrix(rnorm(length(idx2) * nrow(V2), sd = 0.06), length(idx2), nrow(V2))
  )
  Y_noise <- list(
    X1 = matrix(rnorm(length(idx1) * q), length(idx1), q),
    X2 = matrix(rnorm(length(idx2) * q), length(idx2), q)
  )
  Y_anchor <- S_true %*% t(B_true) + matrix(rnorm(N_anchor * q, sd = 0.04), N_anchor, q)

  fit <- response_aligned_mfa(
    Y = Y_noise,
    X = X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    response_alpha = 0,
    anchor_response = Y_anchor,
    anchor_map = list(X1 = idx1, X2 = idx2),
    coupling_lambda = 6,
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  expect_true(all(is.finite(fit$B)))
  expect_equal(dim(fit$S), c(N_anchor, 2))
  expect_false(is.null(fit$anchor_response_fit))
  expect_gt(fit$anchor_response_fit$r2, 0.9)
  expect_lt(mean((fit$S %*% t(fit$B) - Y_anchor)^2), 0.05)
})

test_that("response_aligned_mfa supports soft anchor coupling for mixed anchored and novel rows", {
  set.seed(5)
  K <- 2
  q <- 3
  N_anchor <- 14
  S_true <- scale(matrix(rnorm(N_anchor * K), N_anchor, K), center = TRUE, scale = FALSE)
  B_true <- matrix(rnorm(q * K), q, K)
  V1 <- matrix(rnorm(6 * K), 6, K)
  V2 <- matrix(rnorm(5 * K), 5, K)

  idx1 <- sample(c(seq_len(N_anchor), rep(NA_integer_, 8L)), 46, replace = TRUE)
  idx2 <- sample(c(seq_len(N_anchor), rep(NA_integer_, 8L)), 42, replace = TRUE)
  A1 <- matrix(0, nrow = length(idx1), ncol = N_anchor)
  A2 <- matrix(0, nrow = length(idx2), ncol = N_anchor)
  keep1 <- which(!is.na(idx1))
  keep2 <- which(!is.na(idx2))
  A1[cbind(keep1, idx1[keep1])] <- 1
  A2[cbind(keep2, idx2[keep2])] <- 1

  Z1 <- matrix(rnorm(length(idx1) * K, sd = 0.25), length(idx1), K)
  Z2 <- matrix(rnorm(length(idx2) * K, sd = 0.25), length(idx2), K)
  Z1[keep1, ] <- S_true[idx1[keep1], , drop = FALSE] + matrix(rnorm(length(keep1) * K, sd = 0.03), length(keep1), K)
  Z2[keep2, ] <- S_true[idx2[keep2], , drop = FALSE] + matrix(rnorm(length(keep2) * K, sd = 0.03), length(keep2), K)

  X <- list(
    X1 = Z1 %*% t(V1) + matrix(rnorm(length(idx1) * nrow(V1), sd = 0.06), length(idx1), nrow(V1)),
    X2 = Z2 %*% t(V2) + matrix(rnorm(length(idx2) * nrow(V2), sd = 0.06), length(idx2), nrow(V2))
  )
  Y <- list(
    X1 = Z1 %*% t(B_true) + matrix(rnorm(length(idx1) * q, sd = 0.06), length(idx1), q),
    X2 = Z2 %*% t(B_true) + matrix(rnorm(length(idx2) * q, sd = 0.06), length(idx2), q)
  )

  fit_free <- response_aligned_mfa(
    Y = Y,
    X = X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    coupling_lambda = 0,
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  fit_anchor <- response_aligned_mfa(
    Y = Y,
    X = X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    anchor_map = list(X1 = idx1, X2 = idx2),
    coupling_lambda = c(X1 = 5, X2 = 5),
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  S_free <- muscal:::.ramfa_update_anchor_scores(
    Z_list = fit_free$Z_list,
    anchor_map = list(X1 = A1, X2 = A2),
    anchor_weight = list(X1 = as.numeric(!is.na(idx1)), X2 = as.numeric(!is.na(idx2))),
    coupling_lambda = c(X1 = 1, X2 = 1),
    n_anchor = N_anchor,
    ridge = 1e-8
  )
  gap_free <- mean(c(
    rowSums((fit_free$Z_list$X1[keep1, , drop = FALSE] - A1[keep1, , drop = FALSE] %*% S_free)^2),
    rowSums((fit_free$Z_list$X2[keep2, , drop = FALSE] - A2[keep2, , drop = FALSE] %*% S_free)^2)
  ))
  gap_anchor <- mean(c(
    rowSums((fit_anchor$Z_list$X1[keep1, , drop = FALSE] - fit_anchor$anchor_map$X1[keep1, , drop = FALSE] %*% fit_anchor$S)^2),
    rowSums((fit_anchor$Z_list$X2[keep2, , drop = FALSE] - fit_anchor$anchor_map$X2[keep2, , drop = FALSE] %*% fit_anchor$S)^2)
  ))

  expect_lt(gap_anchor, gap_free)
  expect_equal(dim(fit_anchor$S), c(N_anchor, 2))
  expect_true(all(vapply(fit_anchor$anchor_weight, length, integer(1)) == c(length(idx1), length(idx2))))
})

test_that("response_aligned_mfa treats integer and one-hot anchor maps equivalently", {
  set.seed(51)
  K <- 2
  q <- 3
  N_anchor <- 9
  S_true <- scale(matrix(rnorm(N_anchor * K), N_anchor, K), center = TRUE, scale = FALSE)
  B_true <- matrix(rnorm(q * K), q, K)
  V1 <- matrix(rnorm(5 * K), 5, K)
  V2 <- matrix(rnorm(4 * K), 4, K)

  idx1 <- sample(seq_len(N_anchor), 28, replace = TRUE)
  idx2 <- sample(seq_len(N_anchor), 24, replace = TRUE)
  A1 <- matrix(0, nrow = length(idx1), ncol = N_anchor)
  A2 <- matrix(0, nrow = length(idx2), ncol = N_anchor)
  A1[cbind(seq_along(idx1), idx1)] <- 1
  A2[cbind(seq_along(idx2), idx2)] <- 1

  Z1 <- S_true[idx1, , drop = FALSE] + matrix(rnorm(length(idx1) * K, sd = 0.04), length(idx1), K)
  Z2 <- S_true[idx2, , drop = FALSE] + matrix(rnorm(length(idx2) * K, sd = 0.04), length(idx2), K)
  X <- list(
    X1 = Z1 %*% t(V1) + matrix(rnorm(length(idx1) * nrow(V1), sd = 0.05), length(idx1), nrow(V1)),
    X2 = Z2 %*% t(V2) + matrix(rnorm(length(idx2) * nrow(V2), sd = 0.05), length(idx2), nrow(V2))
  )
  Y <- list(
    X1 = Z1 %*% t(B_true) + matrix(rnorm(length(idx1) * q, sd = 0.05), length(idx1), q),
    X2 = Z2 %*% t(B_true) + matrix(rnorm(length(idx2) * q, sd = 0.05), length(idx2), q)
  )

  fit_idx <- response_aligned_mfa(
    Y = Y,
    X = X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    anchor_map = list(X1 = idx1, X2 = idx2),
    coupling_lambda = 4,
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )
  fit_mat <- response_aligned_mfa(
    Y = Y,
    X = X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    anchor_map = list(X1 = A1, X2 = A2),
    coupling_lambda = 4,
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  expect_equal(fit_idx$anchor_map$X1, A1, tolerance = 1e-12)
  expect_equal(fit_idx$anchor_map$X2, A2, tolerance = 1e-12)
  expect_equal(fit_idx$B, fit_mat$B, tolerance = 1e-10)
  expect_equal(fit_idx$S, fit_mat$S, tolerance = 1e-10)
  expect_equal(fit_idx$Z_list$X1, fit_mat$Z_list$X1, tolerance = 1e-10)
  expect_equal(fit_idx$V_list$X2, fit_mat$V_list$X2, tolerance = 1e-10)
  expect_equal(fit_idx$objective_trace, fit_mat$objective_trace, tolerance = 1e-10)
})

test_that("response_aligned_mfa can use test-time anchor_map to refine prediction", {
  set.seed(6)
  K <- 2
  q <- 3
  N_anchor <- 10
  S_true <- scale(matrix(rnorm(N_anchor * K), N_anchor, K), center = TRUE, scale = FALSE)
  B_true <- matrix(rnorm(q * K), q, K)
  V1 <- matrix(rnorm(6 * K), 6, K)
  V2 <- matrix(rnorm(5 * K), 5, K)
  idx1 <- sample(seq_len(N_anchor), 36, replace = TRUE)
  idx2 <- sample(seq_len(N_anchor), 34, replace = TRUE)

  Z1 <- S_true[idx1, , drop = FALSE] + matrix(rnorm(length(idx1) * K, sd = 0.04), length(idx1), K)
  Z2 <- S_true[idx2, , drop = FALSE] + matrix(rnorm(length(idx2) * K, sd = 0.04), length(idx2), K)

  X <- list(
    X1 = Z1 %*% t(V1) + matrix(rnorm(length(idx1) * nrow(V1), sd = 0.10), length(idx1), nrow(V1)),
    X2 = Z2 %*% t(V2) + matrix(rnorm(length(idx2) * nrow(V2), sd = 0.10), length(idx2), nrow(V2))
  )
  Y <- list(
    X1 = Z1 %*% t(B_true) + matrix(rnorm(length(idx1) * q, sd = 0.08), length(idx1), q),
    X2 = Z2 %*% t(B_true) + matrix(rnorm(length(idx2) * q, sd = 0.08), length(idx2), q)
  )

  fit <- response_aligned_mfa(
    Y = Y,
    X = X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    normalization = "None",
    anchor_map = list(X1 = idx1, X2 = idx2),
    coupling_lambda = 6,
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  new_idx <- sample(seq_len(N_anchor), 12, replace = TRUE)
  Z_new <- S_true[new_idx, , drop = FALSE] + matrix(rnorm(length(new_idx) * K, sd = 0.03), length(new_idx), K)
  X_new <- Z_new %*% t(V1) + matrix(rnorm(length(new_idx) * nrow(V1), sd = 0.18), length(new_idx), nrow(V1))
  Y_new <- Z_new %*% t(B_true) + matrix(rnorm(length(new_idx) * q, sd = 0.08), length(new_idx), q)

  pred_free <- predict(fit, X_new, block = "X1", type = "response")
  pred_anchor <- predict(
    fit,
    X_new,
    block = "X1",
    type = "response",
    new_anchor_map = new_idx
  )

  expect_lte(mean((pred_anchor - Y_new)^2), mean((pred_free - Y_new)^2) + 1e-8)
})
