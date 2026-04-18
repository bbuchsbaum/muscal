library(testthat)
library(muscal)

.projected_grad_norm <- function(S, G) {
  PG <- G - S %*% ((crossprod(S, G) + t(crossprod(S, G))) / 2)
  sqrt(sum(PG^2))
}

test_that("anchored_mfa is the primary interface (linked_mfa is an alias)", {
  set.seed(1)
  Y <- matrix(rnorm(30 * 4), 30, 4)
  X1 <- matrix(rnorm(10 * 6), 10, 6)
  X2 <- matrix(rnorm(12 * 5), 12, 5)
  idx1 <- sample.int(nrow(Y), nrow(X1), replace = FALSE)
  idx2 <- sample.int(nrow(Y), nrow(X2), replace = FALSE)

  fit <- anchored_mfa(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
  expect_s3_class(fit, "anchored_mfa")
  expect_s3_class(fit, "linked_mfa")

  fit_alias <- linked_mfa(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
  expect_s3_class(fit_alias, "linked_mfa")

  S1 <- multivarious::scores(fit)
  S2 <- multivarious::scores(fit_alias)
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
  expect_lt(rel, 1e-12)
})

test_that("linked_mfa validates row_index and returns expected structure", {
  set.seed(1)
  Y <- matrix(rnorm(30 * 4), 30, 4)
  X1 <- matrix(rnorm(10 * 6), 10, 6)
  X2 <- matrix(rnorm(12 * 5), 12, 5)
  idx1 <- sample.int(nrow(Y), nrow(X1), replace = FALSE)
  idx2 <- sample.int(nrow(Y), nrow(X2), replace = FALSE)

  expect_error(linked_mfa(Y, list(X1, X2), list(idx1)))
  expect_error(linked_mfa(Y, list(X1, X2), list(idx1, c(idx2, 1L))))
  expect_error(linked_mfa(Y, list(X1, X2), list(idx1, rep(999L, nrow(X2)))))
  expect_error(linked_mfa(Y, list(X1, X2), list(idx1, c(idx2[1], NA_integer_, idx2[-c(1, 2)]))))

  fit <- linked_mfa(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
  expect_s3_class(fit, "linked_mfa")
  expect_true(inherits(fit, "multiblock_biprojector"))
  expect_equal(length(fit$block_indices), 3)
  expect_equal(names(fit$block_indices), c("Y", "X1", "X2"))
  expect_equal(nrow(multivarious::scores(fit)), nrow(Y))
  expect_equal(multivarious::ncomp(fit), 2)
})

test_that("linked_mfa recovers a strong shared low-rank fit when all blocks share rows", {
  set.seed(2)
  N <- 40
  k <- 2
  S <- qr.Q(qr(matrix(rnorm(N * k), N, k)), complete = FALSE)
  B <- matrix(rnorm(6 * k), 6, k)
  V1 <- matrix(rnorm(8 * k), 8, k)
  V2 <- matrix(rnorm(5 * k), 5, k)
  Y <- S %*% t(B) + matrix(rnorm(N * 6, sd = 0.02), N, 6)
  X1 <- S %*% t(V1) + matrix(rnorm(N * 8, sd = 0.02), N, 8)
  X2 <- S %*% t(V2) + matrix(rnorm(N * 5, sd = 0.02), N, 5)

  # Full overlap mapping
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  fit_lmfa <- linked_mfa(
    Y,
    list(X1 = X1, X2 = X2),
    idx,
    ncomp = 2,
    normalization = "MFA",
    ridge = 1e-10,
    max_iter = 200,
    tol = 1e-10
  )

  expect_equal(crossprod(fit_lmfa$s), diag(k), tolerance = 1e-6)

  yhat <- fit_lmfa$s %*% t(fit_lmfa$B)
  x1hat <- fit_lmfa$s %*% t(fit_lmfa$V_list$X1)
  x2hat <- fit_lmfa$s %*% t(fit_lmfa$V_list$X2)

  expect_lt(mean((Y - yhat)^2), 0.02)
  expect_lt(mean((X1 - x1hat)^2), 0.02)
  expect_lt(mean((X2 - x2hat)^2), 0.02)
})

test_that("feature_groups='colnames' shrinks within-group loading differences", {
  set.seed(3)
  N <- 50
  Y <- matrix(rnorm(N * 3), N, 3)

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

  fit0 <- linked_mfa(Y, list(X1 = X1, X2 = X2), idx, ncomp = 2, feature_groups = "colnames", feature_lambda = 0)
  fit1 <- linked_mfa(Y, list(X1 = X1, X2 = X2), idx, ncomp = 2, feature_groups = "colnames", feature_lambda = 5)

  V0_1 <- fit0$V_list$X1[cn_shared, , drop = FALSE]
  V0_2 <- fit0$V_list$X2[cn_shared, , drop = FALSE]
  V1_1 <- fit1$V_list$X1[cn_shared, , drop = FALSE]
  V1_2 <- fit1$V_list$X2[cn_shared, , drop = FALSE]

  d0 <- rowSums((V0_1 - V0_2)^2)
  d1 <- rowSums((V1_1 - V1_2)^2)

  expect_lt(mean(d1), mean(d0))
})

test_that("anchored_mfa exposes a standard out-of-sample contract", {
  set.seed(5)
  Y <- matrix(rnorm(40 * 4), 40, 4)
  X1 <- matrix(rnorm(14 * 6), 14, 6)
  X2 <- matrix(rnorm(16 * 5), 16, 5)
  idx1 <- sample.int(nrow(Y), nrow(X1), replace = FALSE)
  idx2 <- sample.int(nrow(Y), nrow(X2), replace = FALSE)

  fit <- anchored_mfa(
    Y,
    list(X1 = X1, X2 = X2),
    list(X1 = idx1, X2 = idx2),
    ncomp = 2
  )

  expect_equal(fit$task, "response_prediction")
  expect_equal(fit$fit_spec$method, "anchored_mfa")
  expect_true(fit$fit_spec$refit_supported)
  expect_true(is.list(fit$fit_spec$refit))
  expect_setequal(fit$oos_types, c("response", "scores", "reconstruction"))
  expect_true(all(c("X1", "X2") %in% names(fit$block_preproc)))
  expect_true(inherits(fit$anchor_preproc, "pre_processor"))

  new_rows <- X1[1:3, , drop = FALSE] + matrix(rnorm(18, sd = 0.05), nrow = 3)
  scores_new <- project(fit, new_rows, block = "X1")
  yhat <- predict(fit, new_rows, block = "X1", type = "response")
  xhat <- predict(fit, new_rows, block = "X1", type = "reconstruction")

  expect_equal(scores_new, predict(fit, new_rows, block = "X1", type = "scores"))
  expect_equal(dim(scores_new), c(nrow(new_rows), multivarious::ncomp(fit)))
  expect_equal(dim(yhat), c(nrow(new_rows), ncol(Y)))
  expect_equal(dim(xhat), dim(new_rows))
  expect_error(predict(fit, new_rows, block = "missing"), "Unknown block")
})

test_that("anchored_mfa refit contract can rebuild from stored training data", {
  set.seed(42)
  Y <- matrix(rnorm(36 * 3), 36, 3)
  X1 <- matrix(rnorm(14 * 5), 14, 5)
  X2 <- matrix(rnorm(16 * 4), 16, 4)
  idx1 <- sample.int(nrow(Y), nrow(X1), replace = FALSE)
  idx2 <- sample.int(nrow(Y), nrow(X2), replace = FALSE)

  fit <- anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = idx1, X2 = idx2),
    ncomp = 2
  )

  refit <- fit$fit_spec$refit$fit_fn(fit$fit_spec$refit$data)

  expect_s3_class(refit, "anchored_mfa")
  expect_equal(names(refit$V_list), names(fit$V_list))
  expect_equal(length(refit$row_index), length(fit$row_index))
  expect_true(all(is.finite(refit$sdev)))
})

test_that("anchored_mfa objective trace is finite and non-increasing", {
  set.seed(43)
  Y <- matrix(rnorm(50 * 4), 50, 4)
  X1 <- matrix(rnorm(20 * 6), 20, 6)
  X2 <- matrix(rnorm(24 * 5), 24, 5)
  idx1 <- sample.int(nrow(Y), nrow(X1), replace = TRUE)
  idx2 <- sample.int(nrow(Y), nrow(X2), replace = TRUE)

  fit <- anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = idx1, X2 = idx2),
    ncomp = 2,
    max_iter = 30,
    tol = 0
  )

  obj <- fit$objective_trace
  expect_gt(length(obj), 0)
  expect_true(all(is.finite(obj)))
  expect_lte(max(diff(obj)), 1e-8)
})

test_that("anchored_mfa supports orthonormal score constraint", {
  set.seed(46)
  N <- 40
  Y <- matrix(rnorm(N * 4), N, 4)
  X1 <- matrix(rnorm(18 * 6), 18, 6)
  X2 <- matrix(rnorm(16 * 5), 16, 5)
  idx1 <- sample.int(N, nrow(X1), replace = TRUE)
  idx2 <- sample.int(N, nrow(X2), replace = TRUE)

  fit <- anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = idx1, X2 = idx2),
    ncomp = 2,
    score_constraint = "orthonormal",
    max_iter = 20
  )

  expect_equal(crossprod(fit$s), diag(2), tolerance = 1e-5)
  expect_true(all(is.finite(fit$objective_trace)))
})

test_that("orthonormal MM score step decreases the anchored subproblem objective", {
  set.seed(47)
  N <- 44
  Y <- matrix(rnorm(N * 4), N, 4)
  X1 <- matrix(rnorm(20 * 6), 20, 6)
  X2 <- matrix(rnorm(18 * 5), 18, 5)
  idx <- list(
    X1 = sample.int(N, nrow(X1), replace = TRUE),
    X2 = sample.int(N, nrow(X2), replace = TRUE)
  )

  init <- anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = idx,
    ncomp = 2,
    preproc = multivarious::pass(),
    normalization = "None",
    ridge = 1e-6,
    max_iter = 3,
    tol = 0
  )

  local <- muscal:::.lmfa_score_system(
    Y = Y,
    B = init$B,
    X_list = list(X1 = X1, X2 = X2),
    V_list = init$V_list,
    row_index = init$row_index,
    alpha_y = init$alpha_blocks[[1]],
    alpha_blocks = unname(init$alpha_blocks[-1])
  )
  S0 <- qr.Q(qr(matrix(rnorm(N * 2), N, 2)), complete = FALSE)
  obj0 <- muscal:::.muscal_score_objective_from_system(
    S0,
    local$A_list,
    local$rhs,
    ridge = init$ridge
  )
  grad0 <- muscal:::.lmfa_score_gradient(
    S = S0,
    Y = Y,
    B = init$B,
    X_list = list(X1 = X1, X2 = X2),
    V_list = init$V_list,
    row_index = init$row_index,
    alpha_y = init$alpha_blocks[[1]],
    alpha_blocks = unname(init$alpha_blocks[-1]),
    ridge = init$ridge
  )

  opt <- muscal:::.muscal_stiefel_mm(
    S = S0,
    A_list = local$A_list,
    rhs = local$rhs,
    ridge = init$ridge,
    max_iter = 50,
    tol = 1e-10
  )
  obj1 <- muscal:::.muscal_score_objective_from_system(
    opt$S,
    local$A_list,
    local$rhs,
    ridge = init$ridge
  )
  grad1 <- muscal:::.lmfa_score_gradient(
    S = opt$S,
    Y = Y,
    B = init$B,
    X_list = list(X1 = X1, X2 = X2),
    V_list = init$V_list,
    row_index = init$row_index,
    alpha_y = init$alpha_blocks[[1]],
    alpha_blocks = unname(init$alpha_blocks[-1]),
    ridge = init$ridge
  )
  pg0 <- .projected_grad_norm(S0, grad0)
  pg1 <- .projected_grad_norm(opt$S, grad1)

  expect_equal(crossprod(opt$S), diag(2), tolerance = 1e-6)
  expect_lte(obj1, obj0 + 1e-8)
  expect_lt(pg1, pg0)
})

test_that("anchored_mfa uses anchor information rather than arbitrary Y row assignments", {
  set.seed(44)
  N <- 60
  k <- 2
  S <- scale(matrix(rnorm(N * k), N, k), center = TRUE, scale = FALSE)
  B <- matrix(rnorm(4 * k), 4, k)
  V1 <- matrix(rnorm(8 * k), 8, k)
  V2 <- matrix(rnorm(7 * k), 7, k)

  Y <- S %*% t(B) + matrix(rnorm(N * 4, sd = 0.01), N, 4)
  X1 <- S %*% t(V1) + matrix(rnorm(N * 8, sd = 0.01), N, 8)
  X2 <- S %*% t(V2) + matrix(rnorm(N * 7, sd = 0.01), N, 7)
  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  fit_true <- anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = idx,
    ncomp = k,
    preproc = multivarious::pass(),
    normalization = "None",
    max_iter = 80,
    tol = 1e-8,
    ridge = 1e-10
  )
  fit_shuf <- anchored_mfa(
    Y = Y[sample.int(N), , drop = FALSE],
    X = list(X1 = X1, X2 = X2),
    row_index = idx,
    ncomp = k,
    preproc = multivarious::pass(),
    normalization = "None",
    max_iter = 80,
    tol = 1e-8,
    ridge = 1e-10
  )

  yhat_true <- predict(fit_true, X1, block = "X1", type = "response", preprocess = FALSE)
  yhat_shuf <- predict(fit_shuf, X1, block = "X1", type = "response", preprocess = FALSE)

  mse_true <- mean((yhat_true - Y)^2)
  mse_shuf <- mean((yhat_shuf - Y)^2)

  expect_lt(mse_true, 1e-3)
  expect_gt(mse_shuf, 10 * mse_true)
})
