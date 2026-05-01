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

test_that("aligned_rrr prepares an opaque predictive multifer payload", {
  sim <- .sim_aligned_rrr(n = c(24, 22), p = c(5, 4), q = 3, K = 2, seed = 51)
  fit <- aligned_rrr(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    ridge = 1e-8,
    max_iter = 80,
    tol = 1e-9
  )

  payload <- muscal:::.arrr_multifer_payload(fit)
  prepared <- muscal:::.arrr_multifer_prepare_fit(fit, payload)

  expect_true(muscal:::.arrr_valid_multifer_payload(payload))
  expect_equal(names(payload$X), names(sim$X))
  expect_equal(dim(muscal:::.arrr_project_predictor_scores(prepared, payload)),
               c(sum(vapply(sim$X, nrow, integer(1))), 2))
  expect_equal(dim(muscal:::.arrr_project_response_scores(prepared, payload)),
               c(sum(vapply(sim$Y, nrow, integer(1))), 2))
  expect_true(all(muscal:::.arrr_predictive_roots(prepared, payload) >= 0))
})

test_that("aligned_rrr multifer adapter declares predictive and stability capabilities", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.arrr_multifer_contract_available())

  adapter <- muscal:::.arrr_multifer_adapter()

  for (target in c("component_significance", "variable_stability",
                   "score_stability", "subspace_stability")) {
    expect_true(multifer::adapter_supports(adapter, "adapter", "predictive", target))
  }
  expect_equal(adapter$domains(NULL, NULL), c("predictor", "response"))
  expect_equal(adapter$validity_level, "conditional")
  expect_true("block_rows_exchangeable_within_block" %in% adapter$declared_assumptions)
})

test_that("aligned_rrr multifer predictive roots match weighted SSE oracle", {
  sim <- .sim_aligned_rrr(n = c(24, 22), p = c(5, 4), q = 3, K = 2, seed = 54)
  fit <- aligned_rrr(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    block_weight = c(X1 = 0.8, X2 = 1.2),
    response_weights = list(
      X1 = seq(0.7, 1.1, length.out = nrow(sim$X$X1)),
      X2 = seq(1.2, 0.8, length.out = nrow(sim$X$X2))
    ),
    ridge = 1e-8,
    max_iter = 80,
    tol = 1e-9
  )
  payload <- muscal:::.arrr_multifer_payload(fit)
  prepared <- muscal:::.arrr_multifer_prepare_fit(fit, payload)

  zero_sse <- sum(vapply(names(payload$Y), function(nm) {
    payload$block_weight[[nm]] *
      sum(payload$response_weights[[nm]] * rowSums(payload$Y[[nm]]^2))
  }, numeric(1)))
  sse_for_k <- function(k) {
    Bk <- prepared$B[, seq_len(k), drop = FALSE]
    sum(vapply(names(payload$Y), function(nm) {
      pred <- prepared$Z_list[[nm]][, seq_len(k), drop = FALSE] %*% t(Bk)
      resid <- payload$Y[[nm]] - pred
      payload$block_weight[[nm]] *
        sum(payload$response_weights[[nm]] * rowSums(resid^2))
    }, numeric(1)))
  }
  oracle <- c(
    (zero_sse - sse_for_k(1L)) / zero_sse,
    (sse_for_k(1L) - sse_for_k(2L)) / zero_sse
  )
  oracle <- cummin(pmax(oracle, 0))

  expect_equal(muscal:::.arrr_predictive_roots(prepared, payload), oracle,
               tolerance = 1e-10)
  expect_true(all(diff(muscal:::.arrr_predictive_roots(prepared, payload)) <= 1e-12))
})

test_that("aligned_rrr multifer null/bootstrap/residual payloads preserve contracts", {
  sim <- .sim_aligned_rrr(n = c(22, 20), p = c(5, 4), q = 3, K = 2, seed = 55)
  fit <- aligned_rrr(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    ridge = 1e-8,
    max_iter = 80,
    tol = 1e-9
  )
  payload <- muscal:::.arrr_multifer_payload(fit)
  prepared <- muscal:::.arrr_multifer_prepare_fit(fit, payload)

  null_payload <- muscal:::.arrr_permute_multifer_payload(payload)
  boot <- muscal:::.arrr_bootstrap_multifer_payload(payload)
  residual <- muscal:::.arrr_residualize_multifer_payload(prepared, payload)
  prediction <- muscal:::.arrr_component_prediction(prepared, k = 1L)

  expect_true(muscal:::.arrr_valid_multifer_payload(null_payload))
  expect_true(muscal:::.arrr_valid_multifer_payload(boot$data))
  expect_true(muscal:::.arrr_valid_multifer_payload(residual))
  expect_named(boot$indices, names(payload$X))
  for (nm in names(payload$X)) {
    z <- prepared$Z_list[[nm]][, 1L, drop = FALSE]
    w <- payload$block_weight[[nm]] * payload$response_weights[[nm]]
    expect_equal(crossprod(z * w, residual$X[[nm]]),
                 matrix(0, 1L, ncol(residual$X[[nm]])),
                 tolerance = 1e-8)
    expect_equal(residual$Y[[nm]], payload$Y[[nm]] - prediction[[nm]],
                 tolerance = 1e-10)
  }
})

test_that("aligned_rrr multifer payload alignment is invariant to named block order", {
  sim <- .sim_aligned_rrr(n = c(24, 22), p = c(5, 4), q = 3, K = 2, seed = 56)
  fit <- aligned_rrr(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    ridge = 1e-8,
    max_iter = 80,
    tol = 1e-9
  )

  payload <- muscal:::.arrr_multifer_payload(fit)
  reversed_payload <- muscal:::.arrr_multifer_payload(
    fit,
    data = list(X = rev(sim$X), Y = rev(sim$Y))
  )

  expect_equal(names(reversed_payload$X), names(payload$X))
  expect_equal(
    muscal:::.arrr_predictive_roots(
      muscal:::.arrr_multifer_prepare_fit(fit, reversed_payload),
      reversed_payload
    ),
    muscal:::.arrr_predictive_roots(
      muscal:::.arrr_multifer_prepare_fit(fit, payload),
      payload
    ),
    tolerance = 1e-12
  )
})

test_that("aligned_rrr multifer payload validation rejects invalid weights and shapes", {
  sim <- .sim_aligned_rrr(n = c(20, 18), p = c(5, 4), q = 3, K = 2, seed = 57)
  fit <- aligned_rrr(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    ridge = 1e-8,
    max_iter = 80,
    tol = 1e-9
  )
  payload <- muscal:::.arrr_multifer_payload(fit)

  bad <- payload
  bad$response_weights$X1[1] <- NA_real_
  expect_false(muscal:::.arrr_valid_multifer_payload(bad))

  bad <- payload
  bad$Y$X2 <- bad$Y$X2[-1, , drop = FALSE]
  expect_false(muscal:::.arrr_valid_multifer_payload(bad))

  bad <- payload
  bad$ncomp <- ncol(payload$Y[[1]]) + 1L
  expect_false(muscal:::.arrr_valid_multifer_payload(bad))
})

test_that("infer_aligned_rrr delegates predictive component tests to multifer", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.arrr_multifer_contract_available())

  sim <- .sim_aligned_rrr(n = c(26, 24), p = c(5, 4), q = 3, K = 2, seed = 52)
  fit <- aligned_rrr(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    ridge = 1e-8,
    max_iter = 80,
    tol = 1e-9
  )

  res <- infer_aligned_rrr(
    fit,
    targets = "component_significance",
    B = 7L,
    R = 0L,
    seed = 52L
  )

  expect_s3_class(res, "infer_result")
  expect_match(res$provenance$capabilities,
               "adapter/predictive:component_significance",
               fixed = TRUE)
  expect_match(res$mc$estimand_label, "adapter predictive gain")
  expect_true(nrow(res$component_tests) >= 1L)
  expect_true(all(is.finite(res$component_tests$statistic)))
})

test_that("infer_aligned_rrr uses adapter-owned bootstrap projections for stability", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.arrr_multifer_contract_available())

  sim <- .sim_aligned_rrr(n = c(22, 20), p = c(5, 4), q = 3, K = 2, seed = 53)
  fit <- aligned_rrr(
    Y = sim$Y,
    X = sim$X,
    ncomp = 2,
    preproc = multivarious::pass(),
    response_preproc = multivarious::pass(),
    ridge = 1e-8,
    max_iter = 80,
    tol = 1e-9
  )

  bundle <- infer_aligned_rrr(
    fit,
    targets = c("variable_stability", "score_stability", "subspace_stability"),
    B = 0L,
    R = 4L,
    seed = 53L,
    return_bundle = TRUE
  )

  artifact <- bundle$artifacts$bootstrap
  expect_s3_class(bundle$result, "infer_result")
  expect_s3_class(artifact, "multifer_bootstrap_artifact")
  expect_true(artifact$used_bootstrap_action)
  expect_equal(artifact$score_source, "project_scores")
  expect_setequal(artifact$domains, c("predictor", "response"))
  expect_setequal(names(artifact$reps[[1]]$aligned_scores), c("predictor", "response"))
  expect_gt(nrow(bundle$result$variable_stability), 0L)
  expect_gt(nrow(bundle$result$score_stability), 0L)
  expect_gt(nrow(bundle$result$subspace_stability), 0L)
})
