library(testthat)
library(muscal)

test_that("orient_components makes anchored_mfa sign flips coherent and preserves fitted predictions", {
  set.seed(501)
  Y <- matrix(rnorm(72), nrow = 18)
  X1 <- matrix(rnorm(90), nrow = 18)
  X2 <- matrix(rnorm(72), nrow = 18)

  fit <- anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = 1:18, X2 = 1:18),
    ncomp = 2,
    ridge = 1e-6
  )

  fit_flip <- orient_components(fit, signs = c(-1, 1))

  expect_equal(fit_flip$s[, 1], -fit$s[, 1], tolerance = 1e-8)
  expect_equal(fit_flip$B[, 1], -fit$B[, 1], tolerance = 1e-8)
  expect_equal(fit_flip$V_list$X1[, 1], -fit$V_list$X1[, 1], tolerance = 1e-8)
  expect_equal(fit_flip$V_list$X2[, 1], -fit$V_list$X2[, 1], tolerance = 1e-8)

  new_rows <- X1[1:6, , drop = FALSE]
  score_ref <- predict(fit, new_rows, block = "X1", type = "scores")
  score_flip <- predict(fit_flip, new_rows, block = "X1", type = "scores")
  recon_ref <- predict(fit, new_rows, block = "X1", type = "reconstruction")
  recon_flip <- predict(fit_flip, new_rows, block = "X1", type = "reconstruction")
  resp_ref <- predict(fit, new_rows, block = "X1", type = "response")
  resp_flip <- predict(fit_flip, new_rows, block = "X1", type = "response")

  expect_equal(score_flip[, 1], -score_ref[, 1], tolerance = 1e-6)
  expect_equal(score_flip[, 2], score_ref[, 2], tolerance = 1e-6)
  expect_equal(recon_flip, recon_ref, tolerance = 1e-6)
  expect_equal(resp_flip, resp_ref, tolerance = 1e-6)

  fit_restored <- orient_components(fit_flip, reference = fit)
  expect_equal(fit_restored$s, fit$s, tolerance = 1e-8)
  expect_equal(fit_restored$B, fit$B, tolerance = 1e-8)
  expect_equal(fit_restored$V_list$X1, fit$V_list$X1, tolerance = 1e-8)
})

test_that("orient_components gives deterministic default orientation for anchored_mfa", {
  set.seed(502)
  Y <- matrix(rnorm(80), nrow = 20)
  X1 <- matrix(rnorm(100), nrow = 20)
  X2 <- matrix(rnorm(120), nrow = 20)

  fit <- anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = 1:20, X2 = 1:20),
    ncomp = 2,
    ridge = 1e-6
  )

  oriented_ref <- orient_components(fit)
  oriented_flip <- orient_components(orient_components(fit, signs = c(-1, -1)))

  expect_equal(oriented_flip$s, oriented_ref$s, tolerance = 1e-8)
  expect_equal(oriented_flip$B, oriented_ref$B, tolerance = 1e-8)
  expect_equal(oriented_flip$V_list$X2, oriented_ref$V_list$X2, tolerance = 1e-8)
})

test_that("orient_components updates mfa caches and keeps reconstruction invariant", {
  set.seed(503)
  X1 <- matrix(rnorm(72), nrow = 18)
  X2 <- matrix(rnorm(90), nrow = 18)

  fit <- mfa(list(X1 = X1, X2 = X2), ncomp = 2)
  fit_flip <- orient_components(fit, signs = c(-1, 1))

  expect_equal(fit_flip$s[, 1], -fit$s[, 1], tolerance = 1e-8)
  expect_equal(fit_flip$v[, 1], -fit$v[, 1], tolerance = 1e-8)
  expect_equal(fit_flip$ov[, 1], -fit$ov[, 1], tolerance = 1e-8)
  expect_equal(fit_flip$partial_scores[[1]][, 1], -fit$partial_scores[[1]][, 1], tolerance = 1e-8)

  new_data <- cbind(X1[1:5, , drop = FALSE], X2[1:5, , drop = FALSE])
  score_ref <- predict(fit, new_data, type = "scores")
  score_flip <- predict(fit_flip, new_data, type = "scores")
  recon_ref <- predict(fit, new_data, type = "reconstruction")
  recon_flip <- predict(fit_flip, new_data, type = "reconstruction")

  expect_equal(score_flip[, 1], -score_ref[, 1], tolerance = 1e-6)
  expect_equal(score_flip[, 2], score_ref[, 2], tolerance = 1e-6)
  expect_equal(recon_flip, recon_ref, tolerance = 1e-6)
})

test_that("orient_components updates mcca canonical weights coherently", {
  set.seed(504)
  X1 <- matrix(rnorm(96), nrow = 24)
  X2 <- matrix(rnorm(120), nrow = 24)

  fit <- mcca(list(X1 = X1, X2 = X2), ncomp = 2, ridge = 1e-6)
  fit_flip <- orient_components(fit, signs = c(-1, 1))

  expect_equal(fit_flip$canonical_weights[[1]][, 1], -fit$canonical_weights[[1]][, 1], tolerance = 1e-8)
  expect_equal(fit_flip$canonical_weights[[2]][, 1], -fit$canonical_weights[[2]][, 1], tolerance = 1e-8)
  expect_equal(fit_flip$partial_scores[[1]][, 1], -fit$partial_scores[[1]][, 1], tolerance = 1e-8)

  new_data <- cbind(X1[1:6, , drop = FALSE], X2[1:6, , drop = FALSE])
  score_ref <- predict(fit, new_data, type = "scores")
  score_flip <- predict(fit_flip, new_data, type = "scores")
  recon_ref <- predict(fit, new_data, type = "reconstruction")
  recon_flip <- predict(fit_flip, new_data, type = "reconstruction")

  expect_equal(score_flip[, 1], -score_ref[, 1], tolerance = 1e-6)
  expect_equal(score_flip[, 2], score_ref[, 2], tolerance = 1e-6)
  expect_equal(recon_flip, recon_ref, tolerance = 1e-6)
})

test_that("orient_components keeps ipca projection maps consistent", {
  set.seed(505)
  X1 <- matrix(rnorm(96), nrow = 24)
  X2 <- matrix(rnorm(72), nrow = 24)

  fit <- ipca(
    list(X1 = X1, X2 = X2),
    ncomp = 2,
    lambda = 1,
    method = "gram",
    max_iter = 40,
    tol = 1e-5,
    .return_state = TRUE
  )
  fit_flip <- orient_components(fit, signs = c(-1, 1))

  expect_equal(fit_flip$projection_map$coef_blocks[[1]][, 1], -fit$projection_map$coef_blocks[[1]][, 1], tolerance = 1e-8)
  expect_equal(fit_flip$Sigma_eigenvectors[, 1], -fit$Sigma_eigenvectors[, 1], tolerance = 1e-8)
  expect_equal(fit_flip$warm_state$U[, 1], -fit$warm_state$U[, 1], tolerance = 1e-8)

  new_data <- cbind(X1[1:5, , drop = FALSE], X2[1:5, , drop = FALSE])
  score_ref <- predict(fit, new_data, type = "scores")
  score_flip <- predict(fit_flip, new_data, type = "scores")
  recon_ref <- predict(fit, new_data, type = "reconstruction")
  recon_flip <- predict(fit_flip, new_data, type = "reconstruction")

  expect_equal(score_flip[, 1], -score_ref[, 1], tolerance = 1e-6)
  expect_equal(score_flip[, 2], score_ref[, 2], tolerance = 1e-6)
  expect_equal(recon_flip, recon_ref, tolerance = 1e-6)
})

test_that("orient_components updates penalized_mfa attributes coherently", {
  set.seed(506)
  X1 <- matrix(rnorm(80), nrow = 20)
  X2 <- matrix(rnorm(80), nrow = 20)

  fit <- penalized_mfa(
    list(X1 = X1, X2 = X2),
    ncomp = 2,
    lambda = 1,
    compute_consensus = TRUE,
    penalty_method = "projection",
    max_iter = 4,
    nsteps_inner = 2,
    verbose = FALSE
  )
  fit_flip <- orient_components(fit, signs = c(-1, 1))

  expect_equal(fit_flip$V_list[[1]][, 1], -fit$V_list[[1]][, 1], tolerance = 1e-8)
  expect_equal(fit_flip$scores_list[[2]][, 1], -fit$scores_list[[2]][, 1], tolerance = 1e-8)
  expect_equal(attr(fit_flip, "V_list")[[1]][, 1], -attr(fit, "V_list")[[1]][, 1], tolerance = 1e-8)
  expect_equal(attr(fit_flip, "consensus")[, 1], -attr(fit, "consensus")[, 1], tolerance = 1e-8)
})

test_that("orient_components keeps aligned_rrr predictions invariant", {
  set.seed(508)
  X1 <- matrix(rnorm(48 * 7), nrow = 48)
  X2 <- matrix(rnorm(42 * 6), nrow = 42)
  B_true <- qr.Q(qr(matrix(rnorm(3 * 2), 3, 2)))
  W1 <- matrix(rnorm(7 * 2), 7, 2)
  W2 <- matrix(rnorm(6 * 2), 6, 2)
  Y1 <- X1 %*% W1 %*% t(B_true) + matrix(rnorm(48 * 3, sd = 0.04), 48, 3)
  Y2 <- X2 %*% W2 %*% t(B_true) + matrix(rnorm(42 * 3, sd = 0.04), 42, 3)

  fit <- aligned_rrr(
    Y = list(X1 = Y1, X2 = Y2),
    X = list(X1 = X1, X2 = X2),
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )
  fit_flip <- orient_components(fit, signs = c(-1, 1))

  expect_equal(fit_flip$B[, 1], -fit$B[, 1], tolerance = 1e-8)
  expect_equal(fit_flip$W_list$X1[, 1], -fit$W_list$X1[, 1], tolerance = 1e-8)
  expect_equal(fit_flip$Z_list$X2[, 1], -fit$Z_list$X2[, 1], tolerance = 1e-8)

  new_rows <- X1[1:6, , drop = FALSE]
  score_ref <- predict(fit, new_rows, block = "X1", type = "scores")
  score_flip <- predict(fit_flip, new_rows, block = "X1", type = "scores")
  resp_ref <- predict(fit, new_rows, block = "X1", type = "response")
  resp_flip <- predict(fit_flip, new_rows, block = "X1", type = "response")

  expect_equal(score_flip[, 1], -score_ref[, 1], tolerance = 1e-8)
  expect_equal(score_flip[, 2], score_ref[, 2], tolerance = 1e-8)
  expect_equal(resp_flip, resp_ref, tolerance = 1e-8)

  fit_restored <- orient_components(fit_flip, reference = fit)
  expect_equal(fit_restored$B, fit$B, tolerance = 1e-8)
  expect_equal(fit_restored$W_list$X1, fit$W_list$X1, tolerance = 1e-8)
})

test_that("orient_components keeps response_aligned_mfa anchor predictions coherent", {
  set.seed(509)
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
    X1 = Z1 %*% t(V1) + matrix(rnorm(length(idx1) * nrow(V1), sd = 0.08), length(idx1), nrow(V1)),
    X2 = Z2 %*% t(V2) + matrix(rnorm(length(idx2) * nrow(V2), sd = 0.08), length(idx2), nrow(V2))
  )
  Y <- list(
    X1 = Z1 %*% t(B_true) + matrix(rnorm(length(idx1) * q, sd = 0.06), length(idx1), q),
    X2 = Z2 %*% t(B_true) + matrix(rnorm(length(idx2) * q, sd = 0.06), length(idx2), q)
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
  fit_flip <- orient_components(fit, signs = c(-1, 1))

  expect_equal(fit_flip$B[, 1], -fit$B[, 1], tolerance = 1e-8)
  expect_equal(fit_flip$V_list$X1[, 1], -fit$V_list$X1[, 1], tolerance = 1e-8)
  expect_equal(fit_flip$Z_list$X2[, 1], -fit$Z_list$X2[, 1], tolerance = 1e-8)
  expect_equal(fit_flip$S[, 1], -fit$S[, 1], tolerance = 1e-8)

  new_rows <- X$X1[1:6, , drop = FALSE]
  new_idx <- idx1[1:6]
  score_ref <- predict(fit, new_rows, block = "X1", type = "scores", new_anchor_map = new_idx)
  score_flip <- predict(fit_flip, new_rows, block = "X1", type = "scores", new_anchor_map = new_idx)
  recon_ref <- predict(fit, new_rows, block = "X1", type = "reconstruction", new_anchor_map = new_idx)
  recon_flip <- predict(fit_flip, new_rows, block = "X1", type = "reconstruction", new_anchor_map = new_idx)
  resp_ref <- predict(fit, new_rows, block = "X1", type = "response", new_anchor_map = new_idx)
  resp_flip <- predict(fit_flip, new_rows, block = "X1", type = "response", new_anchor_map = new_idx)

  expect_equal(score_flip[, 1], -score_ref[, 1], tolerance = 1e-8)
  expect_equal(score_flip[, 2], score_ref[, 2], tolerance = 1e-8)
  expect_equal(recon_flip, recon_ref, tolerance = 1e-8)
  expect_equal(resp_flip, resp_ref, tolerance = 1e-8)

  fit_restored <- orient_components(fit_flip, reference = fit)
  expect_equal(fit_restored$B, fit$B, tolerance = 1e-8)
  expect_equal(fit_restored$S, fit$S, tolerance = 1e-8)
  expect_equal(fit_restored$V_list$X1, fit$V_list$X1, tolerance = 1e-8)
})

test_that("loading_reliability aligns bootstrap loading signs to the reference", {
  set.seed(507)
  X1 <- matrix(rnorm(72), nrow = 18)
  X2 <- matrix(rnorm(72), nrow = 18)
  fit <- mfa(list(X1 = X1, X2 = X2), ncomp = 2)

  ref <- fit$v
  boot_loadings <- list(
    ref,
    sweep(ref, 2, c(-1, 1), `*`),
    sweep(ref, 2, c(1, -1), `*`)
  )

  res <- loading_reliability(fit, method = "bootstrap", boot_loadings = boot_loadings, V_ref = ref)

  expect_equal(unname(res$mean), unname(ref), tolerance = 1e-10)
  expect_true(all(res$sd < 1e-10))
  expect_true(all(res$sign_consistency == 1))
})
