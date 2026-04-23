library(testthat)
library(muscal)

.sim_aligned_rrr <- function(n = c(48, 42),
                             p = c(7, 6),
                             q = 3,
                             K = 2,
                             noise_y = 0.04,
                             seed = 1) {
  set.seed(seed)
  B <- qr.Q(qr(matrix(rnorm(q * K), q, K)))
  W_list <- lapply(p, function(pk) matrix(rnorm(pk * K), pk, K))
  X <- lapply(seq_along(n), function(k) matrix(rnorm(n[[k]] * p[[k]]), n[[k]], p[[k]]))
  Y <- lapply(seq_along(n), function(k) {
    X[[k]] %*% W_list[[k]] %*% t(B) +
      matrix(rnorm(n[[k]] * q, sd = noise_y), n[[k]], q)
  })

  names(X) <- paste0("X", seq_along(X))
  names(Y) <- names(X)
  names(W_list) <- names(X)

  list(X = X, Y = Y, W_list = W_list, B = B)
}

.aligned_rrr_weighted_gram <- function(fit) {
  Reduce(
    "+",
    lapply(seq_along(fit$Z_list), function(k) {
      wk <- sqrt(fit$block_weight[[k]] * fit$response_weights[[k]])
      crossprod(fit$Z_list[[k]] * wk)
    })
  )
}

.aligned_rrr_weighted_ls_coef <- function(X, Y, w) {
  Xw <- X * sqrt(w)
  Yw <- Y * sqrt(w)
  solve(crossprod(Xw), crossprod(Xw, Yw))
}

.aligned_rrr_expect_stable_descent <- function(trace,
                                               final_rel_tol = 5e-5,
                                               tail_rel_tol = 2e-4,
                                               max_tail_increase = 2e-4,
                                               tail_n = 5L) {
  expect_true(length(trace) >= 2L)
  expect_true(all(is.finite(trace)))
  scale <- max(1, abs(trace))
  expect_lte(trace[[length(trace)]] - min(trace), final_rel_tol * scale)
  tail_trace <- utils::tail(trace, min(length(trace), tail_n))
  tail_diff <- diff(tail_trace)
  expect_lte(max(abs(tail_diff)), tail_rel_tol * scale)
  expect_lte(max(tail_diff), max_tail_increase * scale)
}

test_that("aligned_rrr validates inputs and exposes the standard contract", {
  set.seed(1)
  X1 <- matrix(rnorm(30), nrow = 10)
  X2 <- matrix(rnorm(40), nrow = 10)
  Y1 <- matrix(rnorm(20), nrow = 10)
  Y2 <- matrix(rnorm(20), nrow = 10)

  expect_error(aligned_rrr(Y = list(Y1), X = list(X1, X2)))
  expect_error(aligned_rrr(Y = list(X1 = Y1, X2 = Y2[1:9, , drop = FALSE]), X = list(X1 = X1, X2 = X2)))
  expect_error(aligned_rrr(
    Y = list(X1 = Y1, X2 = matrix(rnorm(30), 10, 3)),
    X = list(X1 = X1, X2 = X2)
  ))
  expect_error(aligned_rrr(
    Y = list(X1 = Y1, X2 = Y2),
    X = list(X1 = X1, X2 = X2),
    response_preproc = list(multivarious::center(), multivarious::center())
  ))
  expect_error(aligned_rrr(
    Y = list(X1 = Y1, X2 = Y2),
    X = list(X1 = X1, X2 = X2),
    ncomp = 3
  ), "ncomp")

  fit <- aligned_rrr(
    Y = list(X1 = Y1, X2 = Y2),
    X = list(X1 = X1, X2 = X2),
    ncomp = 2
  )

  expect_s3_class(fit, "aligned_rrr")
  expect_true(inherits(fit, "multiblock_biprojector"))
  expect_equal(fit$task, "response_prediction")
  expect_equal(fit$fit_spec$method, "aligned_rrr")
  expect_true(fit$fit_spec$refit_supported)
  expect_setequal(fit$oos_types, c("response", "scores"))
  expect_equal(names(fit$score_index), c("X1", "X2"))
  expect_equal(fit$score_representation, "stacked_block_scores")
  expect_equal(multivarious::scores(fit), do.call(rbind, fit$Z_list))
})

test_that("aligned_rrr recovers a shared reduced-rank response model", {
  sim <- .sim_aligned_rrr(seed = 2)

  fit <- aligned_rrr(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    block_weight = c(X1 = 0.8, X2 = 1.4),
    response_weights = list(
      X1 = seq(0.6, 1.2, length.out = nrow(sim$X$X1)),
      X2 = seq(1.1, 0.7, length.out = nrow(sim$X$X2))
    ),
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  expect_equal(.aligned_rrr_weighted_gram(fit), diag(2), tolerance = 1e-5)

  mse_y <- mean(unlist(lapply(seq_along(sim$Y), function(k) {
    Y_hat <- sim$X[[k]] %*% fit$W_list[[k]] %*% t(fit$B)
    (sim$Y[[k]] - Y_hat)^2
  })))

  expect_lt(mse_y, 0.05)
  expect_true(all(is.finite(fit$objective_trace)))
})

test_that("aligned_rrr predict() and project() follow the shared score contract", {
  sim <- .sim_aligned_rrr(seed = 3)

  fit <- aligned_rrr(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  set.seed(33)
  X_new <- matrix(rnorm(12 * nrow(sim$W_list$X1)), 12, nrow(sim$W_list$X1))
  Y_new <- X_new %*% sim$W_list$X1 %*% t(sim$B) +
    matrix(rnorm(12 * nrow(sim$B), sd = 0.04), 12, nrow(sim$B))

  score_x <- project(fit, X_new, block = "X1")
  resp_x <- predict(fit, X_new, block = "X1", type = "response")

  expect_equal(dim(score_x), c(nrow(X_new), multivarious::ncomp(fit)))
  expect_equal(dim(resp_x), dim(Y_new))
  expect_equal(score_x, predict(fit, X_new, block = "X1", type = "scores"))
  expect_lt(mean((resp_x - Y_new)^2), 0.1)
})

test_that("aligned_rrr reduces to weighted least squares when ncomp spans the response", {
  set.seed(31)
  n1 <- 42
  n2 <- 38
  p1 <- 5
  p2 <- 4
  q <- 2

  X1 <- matrix(rnorm(n1 * p1), n1, p1)
  X2 <- matrix(rnorm(n2 * p2), n2, p2)
  C1 <- matrix(rnorm(p1 * q), p1, q)
  C2 <- matrix(rnorm(p2 * q), p2, q)
  Y1 <- X1 %*% C1 + matrix(rnorm(n1 * q, sd = 0.05), n1, q)
  Y2 <- X2 %*% C2 + matrix(rnorm(n2 * q, sd = 0.05), n2, q)
  w1 <- seq(0.5, 1.5, length.out = n1)
  w2 <- seq(1.4, 0.6, length.out = n2)

  fit <- aligned_rrr(
    Y = list(X1 = Y1, X2 = Y2),
    X = list(X1 = X1, X2 = X2),
    ncomp = q,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    response_weights = list(X1 = w1, X2 = w2),
    ridge = 1e-10,
    max_iter = 120,
    tol = 1e-10
  )

  coef1 <- .aligned_rrr_weighted_ls_coef(X1, Y1, w1)
  coef2 <- .aligned_rrr_weighted_ls_coef(X2, Y2, w2)
  yhat1_ref <- X1 %*% coef1
  yhat2_ref <- X2 %*% coef2
  yhat1_fit <- fit$Z_list$X1 %*% t(fit$B)
  yhat2_fit <- fit$Z_list$X2 %*% t(fit$B)

  expect_equal(yhat1_fit, yhat1_ref, tolerance = 1e-7)
  expect_equal(yhat2_fit, yhat2_ref, tolerance = 1e-7)
  .aligned_rrr_expect_stable_descent(fit$objective_trace)
})

test_that("aligned_rrr ignores blocks with zero effective response weight", {
  set.seed(32)
  n1 <- 34
  n2 <- 30
  p1 <- 5
  p2 <- 4
  q <- 2

  X1 <- matrix(rnorm(n1 * p1), n1, p1)
  Y1 <- matrix(rnorm(n1 * q), n1, q)
  X2 <- matrix(rnorm(n2 * p2), n2, p2)
  Y2 <- matrix(rnorm(n2 * q), n2, q)

  fit_ref <- suppressWarnings(aligned_rrr(
    Y = list(X1 = Y1, X2 = Y2),
    X = list(X1 = X1, X2 = X2),
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    response_weights = list(X1 = rep(1, n1), X2 = rep(0, n2)),
    ridge = 1e-8,
    max_iter = 100,
    tol = 1e-10
  ))

  fit_perturbed <- suppressWarnings(aligned_rrr(
    Y = list(
      X1 = Y1,
      X2 = matrix(rnorm(n2 * q, sd = 50), n2, q)
    ),
    X = list(
      X1 = X1,
      X2 = matrix(rnorm(n2 * p2, sd = 50), n2, p2)
    ),
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    response_weights = list(X1 = rep(1, n1), X2 = rep(0, n2)),
    ridge = 1e-8,
    max_iter = 100,
    tol = 1e-10
  ))

  expect_equal(fit_ref$B, fit_perturbed$B, tolerance = 1e-5)
  expect_equal(fit_ref$W_list$X1, fit_perturbed$W_list$X1, tolerance = 1e-6)
  expect_equal(fit_ref$Z_list$X1, fit_perturbed$Z_list$X1, tolerance = 1e-6)
})

test_that("aligned_rrr is invariant to block order and within-block row permutations", {
  sim <- .sim_aligned_rrr(n = c(34, 30), p = c(6, 5), q = 3, K = 2, seed = 41)
  block_weight <- c(X1 = 0.8, X2 = 1.3)
  response_weights <- list(
    X1 = seq(0.6, 1.2, length.out = nrow(sim$X$X1)),
    X2 = seq(1.1, 0.7, length.out = nrow(sim$X$X2))
  )

  fit <- aligned_rrr(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    block_weight = block_weight,
    response_weights = response_weights,
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  fit_swapped <- aligned_rrr(
    Y = sim$Y[c("X2", "X1")],
    X = sim$X[c("X2", "X1")],
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    block_weight = block_weight[c("X2", "X1")],
    response_weights = response_weights[c("X2", "X1")],
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  perm1 <- sample.int(nrow(sim$X$X1))
  perm2 <- sample.int(nrow(sim$X$X2))
  fit_perm <- aligned_rrr(
    Y = list(
      X1 = sim$Y$X1[perm1, , drop = FALSE],
      X2 = sim$Y$X2[perm2, , drop = FALSE]
    ),
    X = list(
      X1 = sim$X$X1[perm1, , drop = FALSE],
      X2 = sim$X$X2[perm2, , drop = FALSE]
    ),
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    block_weight = block_weight,
    response_weights = list(
      X1 = response_weights$X1[perm1],
      X2 = response_weights$X2[perm2]
    ),
    ridge = 1e-8,
    max_iter = 120,
    tol = 1e-10
  )

  expect_equal(fit$B, fit_swapped$B, tolerance = 5e-3)
  expect_equal(fit$W_list$X1, fit_swapped$W_list$X1, tolerance = 1e-4)
  expect_equal(fit$W_list$X2, fit_swapped$W_list$X2, tolerance = 1e-4)
  expect_equal(fit$Z_list$X1, fit_swapped$Z_list$X1, tolerance = 1e-4)
  expect_equal(fit$Z_list$X2, fit_swapped$Z_list$X2, tolerance = 1e-4)

  expect_equal(fit$B, fit_perm$B, tolerance = 1e-3)
  expect_equal(fit$W_list$X1, fit_perm$W_list$X1, tolerance = 1e-4)
  expect_equal(fit$W_list$X2, fit_perm$W_list$X2, tolerance = 1e-4)
  expect_equal(fit$Z_list$X1[perm1, , drop = FALSE], fit_perm$Z_list$X1, tolerance = 1e-4)
  expect_equal(fit$Z_list$X2[perm2, , drop = FALSE], fit_perm$Z_list$X2, tolerance = 1e-4)

  expect_equal(
    predict(fit, sim$X$X1, block = "X1", type = "response", preprocess = FALSE),
    predict(fit_swapped, sim$X$X1, block = "X1", type = "response", preprocess = FALSE),
    tolerance = 1e-4
  )
  expect_equal(
    predict(fit, sim$X$X2, block = "X2", type = "response", preprocess = FALSE),
    predict(fit_swapped, sim$X$X2, block = "X2", type = "response", preprocess = FALSE),
    tolerance = 1e-4
  )
})

test_that("aligned_rrr remains finite for nearly singular predictors and extreme row weights", {
  set.seed(42)
  n1 <- 40
  n2 <- 36
  q <- 3
  K <- 2

  base1 <- matrix(rnorm(n1 * 2), n1, 2)
  base2 <- matrix(rnorm(n2 * 2), n2, 2)
  X1 <- cbind(
    base1[, 1],
    base1[, 1] + 1e-8 * rnorm(n1),
    base1[, 2],
    base1[, 2] - 1e-8 * rnorm(n1),
    rnorm(n1)
  )
  X2 <- cbind(
    base2[, 1],
    base2[, 1] - 1e-8 * rnorm(n2),
    base2[, 2],
    base2[, 2] + 1e-8 * rnorm(n2)
  )
  B <- matrix(rnorm(q * K), q, K)
  W1 <- matrix(rnorm(ncol(X1) * K), ncol(X1), K)
  W2 <- matrix(rnorm(ncol(X2) * K), ncol(X2), K)
  Y1 <- X1 %*% W1 %*% t(B) + matrix(rnorm(n1 * q, sd = 0.05), n1, q)
  Y2 <- X2 %*% W2 %*% t(B) + matrix(rnorm(n2 * q, sd = 0.05), n2, q)
  rw1 <- exp(seq(log(1e-3), log(1e3), length.out = n1))
  rw2 <- exp(seq(log(1e3), log(1e-3), length.out = n2))

  fit <- aligned_rrr(
    Y = list(X1 = Y1, X2 = Y2),
    X = list(X1 = X1, X2 = X2),
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    response_weights = list(X1 = rw1, X2 = rw2),
    ridge = 1e-5,
    max_iter = 120,
    tol = 1e-10
  )

  pred <- predict(fit, X1[1:8, , drop = FALSE], block = "X1", type = "response")

  expect_true(all(is.finite(fit$B)))
  expect_true(all(vapply(fit$W_list, function(W) all(is.finite(W)), logical(1))))
  expect_true(all(vapply(fit$Z_list, function(Z) all(is.finite(Z)), logical(1))))
  expect_true(all(is.finite(pred)))
  .aligned_rrr_expect_stable_descent(
    fit$objective_trace,
    final_rel_tol = 1e-4,
    tail_rel_tol = 2e-4,
    max_tail_increase = 2e-4
  )
})

test_that("aligned_rrr stores parsed row weights and exposes deterministic refits", {
  sim <- .sim_aligned_rrr(n = c(28, 24), p = c(6, 5), q = 3, K = 2, seed = 14)

  fit <- aligned_rrr(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    block_weight = c(X1 = 0.9, X2 = 1.1),
    response_weights = c(X1 = 0.7, X2 = 1.3),
    ridge = 1e-8,
    max_iter = 100,
    tol = 1e-10
  )

  refit_data <- fit$fit_spec$refit$data
  expect_true(is.list(refit_data$response_weights))
  expect_equal(length(refit_data$response_weights$X1), nrow(sim$X$X1))
  expect_equal(length(refit_data$response_weights$X2), nrow(sim$X$X2))

  refit <- fit$fit_spec$refit$fit_fn(refit_data)
  expect_s3_class(refit, "aligned_rrr")
  expect_equal(refit$B, fit$B, tolerance = 1e-8)
  expect_equal(refit$W_list$X1, fit$W_list$X1, tolerance = 1e-8)
  expect_equal(refit$W_list$X2, fit$W_list$X2, tolerance = 1e-8)
})
