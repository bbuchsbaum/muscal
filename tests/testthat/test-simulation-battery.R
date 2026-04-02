library(testthat)
library(muscal)

skip_if(
  tolower(Sys.getenv("MUSCAL_SKIP_SIM_BATTERY", unset = "false")) %in% c("1", "true", "yes"),
  "Simulation battery disabled via MUSCAL_SKIP_SIM_BATTERY."
)

.sim_battery_values <- function(default, strict = default) {
  strict_flag <- tolower(Sys.getenv("MUSCAL_SIM_BATTERY_STRICT", unset = "false"))
  if (strict_flag %in% c("1", "true", "yes")) strict else default
}

.sim_subspace_rel_diff <- function(S_fit, S_true) {
  P_fit <- S_fit %*% solve(crossprod(S_fit), t(S_fit))
  P_true <- S_true %*% solve(crossprod(S_true), t(S_true))
  norm(P_fit - P_true, type = "F") / (norm(P_true, type = "F") + 1e-12)
}

.sim_mfa_recovery_error <- function(seed, noise_sd, n = 70, k = 2, p_vec = c(9, 11, 10)) {
  set.seed(seed)
  Z <- matrix(rnorm(n * k), nrow = n, ncol = k)
  blocks <- lapply(p_vec, function(p) {
    A <- matrix(rnorm(p * k), nrow = p, ncol = k)
    Z %*% t(A) + matrix(rnorm(n * p, sd = noise_sd), nrow = n, ncol = p)
  })

  fit <- suppressMessages(
    mfa(
      blocks,
      ncomp = k,
      preproc = multivarious::pass(),
      normalization = "None"
    )
  )

  .sim_subspace_rel_diff(multivarious::scores(fit), Z)
}

.sim_anchor_prediction_mse <- function(seed,
                                       n_train,
                                       n_train_max = 60,
                                       n_test = 100,
                                       k = 2,
                                       noise_sd = 0.15,
                                       p1 = 5,
                                       p2 = 4,
                                       q = 3) {
  set.seed(seed)

  B <- matrix(rnorm(q * k), nrow = q, ncol = k)
  V1 <- matrix(rnorm(p1 * k), nrow = p1, ncol = k)
  V2 <- matrix(rnorm(p2 * k), nrow = p2, ncol = k)

  S_train_full <- scale(matrix(rnorm(n_train_max * k), n_train_max, k), center = TRUE, scale = FALSE)
  S_test <- scale(matrix(rnorm(n_test * k), n_test, k), center = TRUE, scale = FALSE)

  Y_train_full <- S_train_full %*% t(B) + matrix(rnorm(n_train_max * q, sd = noise_sd), n_train_max, q)
  X1_train_full <- S_train_full %*% t(V1) + matrix(rnorm(n_train_max * p1, sd = noise_sd), n_train_max, p1)
  X2_train_full <- S_train_full %*% t(V2) + matrix(rnorm(n_train_max * p2, sd = noise_sd), n_train_max, p2)

  Y_test <- S_test %*% t(B) + matrix(rnorm(n_test * q, sd = noise_sd), n_test, q)
  X1_test <- S_test %*% t(V1) + matrix(rnorm(n_test * p1, sd = noise_sd), n_test, p1)

  idx <- seq_len(n_train)
  fit <- suppressMessages(
    anchored_mfa(
      Y = Y_train_full[idx, , drop = FALSE],
      X = list(
        X1 = X1_train_full[idx, , drop = FALSE],
        X2 = X2_train_full[idx, , drop = FALSE]
      ),
      row_index = list(X1 = seq_len(n_train), X2 = seq_len(n_train)),
      ncomp = k,
      preproc = multivarious::pass(),
      normalization = "None",
      max_iter = 60,
      tol = 1e-8,
      ridge = 1e-8
    )
  )

  Y_hat <- predict(fit, X1_test, block = "X1", type = "response", preprocess = FALSE)
  mean((Y_hat - Y_test)^2)
}

.sim_aligned_mapping_errors <- function(seed,
                                        N = 80,
                                        n_each = 45,
                                        k = 2,
                                        noise_sd = 0.04) {
  set.seed(seed)

  S <- scale(matrix(rnorm(N * k), N, k), center = TRUE, scale = FALSE)
  idx1 <- sample.int(N, n_each, replace = TRUE)
  idx2 <- sample.int(N, n_each, replace = TRUE)

  V1 <- matrix(rnorm(8 * k), nrow = 8, ncol = k)
  V2 <- matrix(rnorm(7 * k), nrow = 7, ncol = k)
  X1 <- S[idx1, , drop = FALSE] %*% t(V1) + matrix(rnorm(length(idx1) * 8, sd = noise_sd), length(idx1), 8)
  X2 <- S[idx2, , drop = FALSE] %*% t(V2) + matrix(rnorm(length(idx2) * 7, sd = noise_sd), length(idx2), 7)

  fit_true <- suppressMessages(
    aligned_mfa(
      list(X1 = X1, X2 = X2),
      list(X1 = idx1, X2 = idx2),
      N = N,
      ncomp = k,
      preproc = multivarious::pass(),
      normalization = "None",
      max_iter = 50,
      tol = 1e-8,
      ridge = 1e-10
    )
  )

  fit_wrong <- suppressMessages(
    aligned_mfa(
      list(X1 = X1, X2 = X2),
      list(X1 = sample(idx1), X2 = sample(idx2)),
      N = N,
      ncomp = k,
      preproc = multivarious::pass(),
      normalization = "None",
      max_iter = 50,
      tol = 1e-8,
      ridge = 1e-10
    )
  )

  c(
    true = .sim_subspace_rel_diff(multivarious::scores(fit_true), S),
    wrong = .sim_subspace_rel_diff(multivarious::scores(fit_wrong), S)
  )
}

.sim_aligned_mcca_mapping_cc <- function(seed,
                                         N = 75,
                                         k = 2,
                                         p_vec = c(15, 17, 14),
                                         n_vec = c(55, 52, 58),
                                         noise_sd = 0.02) {
  set.seed(seed)

  Z <- matrix(rnorm(N * k), nrow = N, ncol = k)
  idx <- lapply(n_vec, function(nk) sample.int(N, nk, replace = TRUE))
  blocks <- Map(function(p, id) {
    A <- matrix(rnorm(p * k), nrow = p, ncol = k)
    Z[id, , drop = FALSE] %*% t(A) +
      matrix(rnorm(length(id) * p, sd = noise_sd), nrow = length(id), ncol = p)
  }, p_vec, idx)

  fit_true <- suppressMessages(
    aligned_mcca(blocks, idx, N = N, ncomp = k, ridge = 1e-6)
  )
  fit_wrong <- suppressMessages(
    aligned_mcca(blocks, lapply(idx, sample), N = N, ncomp = k, ridge = 1e-6)
  )

  cc_true <- cancor(
    scale(Z, center = TRUE, scale = FALSE),
    scale(multivarious::scores(fit_true), center = TRUE, scale = FALSE)
  )$cor
  cc_wrong <- cancor(
    scale(Z, center = TRUE, scale = FALSE),
    scale(multivarious::scores(fit_wrong), center = TRUE, scale = FALSE)
  )$cor

  c(true = mean(cc_true[seq_len(k)]), wrong = mean(cc_wrong[seq_len(k)]))
}

.sim_graph_anchor_prediction <- function(seed,
                                         mode = c("none", "correct", "wrong"),
                                         n_train = 16,
                                         n_test = 100,
                                         k = 2,
                                         q = 3,
                                         noise_sd = 0.18,
                                         lambda = 8) {
  mode <- match.arg(mode)
  set.seed(seed)

  shared <- paste0("f", seq_len(6))
  u1 <- paste0("u1_", seq_len(2))
  u2 <- paste0("u2_", seq_len(2))
  cn1 <- c(shared, u1)
  cn2 <- c(shared, u2)

  B <- matrix(rnorm(q * k), nrow = q, ncol = k)
  V_shared <- matrix(rnorm(length(shared) * k), nrow = length(shared), ncol = k)
  V1 <- rbind(V_shared, matrix(rnorm(length(u1) * k, sd = 0.2), nrow = length(u1), ncol = k))
  V2 <- rbind(V_shared, matrix(rnorm(length(u2) * k, sd = 0.2), nrow = length(u2), ncol = k))

  S_train <- scale(matrix(rnorm(n_train * k), n_train, k), center = TRUE, scale = FALSE)
  S_test <- scale(matrix(rnorm(n_test * k), n_test, k), center = TRUE, scale = FALSE)

  Y_train <- S_train %*% t(B) + matrix(rnorm(n_train * q, sd = noise_sd), n_train, q)
  X1_train <- S_train %*% t(V1) + matrix(rnorm(n_train * nrow(V1), sd = noise_sd), n_train, nrow(V1))
  X2_train <- S_train %*% t(V2) + matrix(rnorm(n_train * nrow(V2), sd = noise_sd), n_train, nrow(V2))
  Y_test <- S_test %*% t(B) + matrix(rnorm(n_test * q, sd = noise_sd), n_test, q)
  X1_test <- S_test %*% t(V1) + matrix(rnorm(n_test * nrow(V1), sd = noise_sd), n_test, nrow(V1))

  colnames(X1_train) <- colnames(X1_test) <- cn1
  colnames(X2_train) <- cn2

  feature_graph <- NULL
  graph_lambda <- 0
  if (identical(mode, "correct")) {
    feature_graph <- "colnames"
    graph_lambda <- lambda
  } else if (identical(mode, "wrong")) {
    feature_graph <- data.frame(
      block1 = "X1",
      feature1 = shared,
      block2 = "X2",
      feature2 = rev(shared),
      weight = 1,
      stringsAsFactors = FALSE
    )
    graph_lambda <- lambda
  }

  fit <- suppressMessages(
    graph_anchored_mfa(
      Y = Y_train,
      X = list(X1 = X1_train, X2 = X2_train),
      row_index = list(X1 = seq_len(n_train), X2 = seq_len(n_train)),
      ncomp = k,
      preproc = multivarious::pass(),
      normalization = "None",
      feature_graph = feature_graph,
      graph_lambda = graph_lambda,
      max_iter = 60,
      tol = 1e-8,
      ridge = 1e-8
    )
  )

  Y_hat <- predict(fit, X1_test, block = "X1", type = "response", preprocess = FALSE)
  mean((Y_hat - Y_test)^2)
}

test_that("simulation battery: MFA recovery worsens monotonically across a noise grid", {
  seeds <- .sim_battery_values(981:986, 981:990)
  noise_grid <- .sim_battery_values(c(0.02, 0.10, 0.30), c(0.02, 0.08, 0.18, 0.30))

  grid <- expand.grid(seed = seeds, noise_sd = noise_grid, KEEP.OUT.ATTRS = FALSE)
  grid$error <- mapply(.sim_mfa_recovery_error, grid$seed, grid$noise_sd)

  medians <- stats::aggregate(error ~ noise_sd, data = grid, FUN = stats::median)
  rho <- stats::cor(grid$noise_sd, grid$error, method = "spearman")

  expect_true(all(diff(medians$error) > 0))
  expect_gt(rho, 0.75)
})

test_that("simulation battery: anchored_mfa out-of-sample prediction improves with more training rows", {
  seeds <- .sim_battery_values(951:956, 951:960)

  mse_small <- vapply(seeds, .sim_anchor_prediction_mse, numeric(1), n_train = 12, n_train_max = 60)
  mse_large <- vapply(seeds, .sim_anchor_prediction_mse, numeric(1), n_train = 60, n_train_max = 60)

  expect_lt(stats::median(mse_large), stats::median(mse_small))
  expect_gte(sum(mse_large < mse_small), ceiling(0.75 * length(seeds)))
})

test_that("simulation battery: aligned_mfa with repeated indices benefits from the correct row mapping", {
  seeds <- .sim_battery_values(921:924, 921:928)
  errs <- t(vapply(seeds, .sim_aligned_mapping_errors, numeric(2)))

  expect_true(all(is.finite(errs)))
  expect_gte(sum(errs[, "true"] < errs[, "wrong"]), ceiling(0.75 * length(seeds)))
  expect_gt(stats::median(errs[, "wrong"]) - stats::median(errs[, "true"]), 0.3)
})

test_that("simulation battery: aligned_mcca depends on the correct repeated-row mapping", {
  seeds <- .sim_battery_values(1021:1024, 1021:1028)
  vals <- t(vapply(seeds, .sim_aligned_mcca_mapping_cc, numeric(2)))

  expect_true(all(is.finite(vals)))
  expect_gte(sum(vals[, "true"] > vals[, "wrong"]), ceiling(0.75 * length(seeds)))
  expect_gt(stats::median(vals[, "true"]) - stats::median(vals[, "wrong"]), 0.25)
})

test_that("simulation battery: graph_anchored_mfa benefits from a correct feature graph", {
  seeds <- .sim_battery_values(1011:1014, 1011:1018)
  mse_none <- vapply(seeds, .sim_graph_anchor_prediction, numeric(1), mode = "none")
  mse_correct <- vapply(seeds, .sim_graph_anchor_prediction, numeric(1), mode = "correct")
  mse_wrong <- vapply(seeds, .sim_graph_anchor_prediction, numeric(1), mode = "wrong")

  expect_lt(stats::median(mse_correct), stats::median(mse_none))
  expect_lt(stats::median(mse_correct), stats::median(mse_wrong))
  expect_gte(sum(mse_correct <= mse_none), ceiling(0.75 * length(seeds)))
  expect_gte(sum(mse_correct < mse_wrong), ceiling(0.75 * length(seeds)))
  expect_gt(stats::median(mse_wrong) - stats::median(mse_correct), 0.05)
})
