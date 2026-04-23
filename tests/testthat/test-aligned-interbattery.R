library(testthat)
library(muscal)

.aib_subspace_rel <- function(S1, S2) {
  P1 <- tcrossprod(scale(S1, center = TRUE, scale = FALSE))
  P2 <- tcrossprod(scale(S2, center = TRUE, scale = FALSE))
  sqrt(sum((P1 - P2)^2)) / max(1e-12, sqrt(sum(P1^2)))
}

.aib_offdiag_cov_norm <- function(score_list) {
  S <- do.call(rbind, score_list)
  C <- crossprod(S) / nrow(S)
  diag(C) <- 0
  sqrt(sum(C^2))
}

.aib_offdiag_cov_norm_matrix <- function(S) {
  C <- crossprod(S) / nrow(S)
  diag(C) <- 0
  sqrt(sum(C^2))
}

.aib_cov_identity_distance <- function(S) {
  C <- crossprod(S) / nrow(S)
  sqrt(sum((C - diag(ncol(C)))^2))
}

.sim_aib_subject <- function(N, P1, P2, Q1, noise = 0.05, include_d2 = TRUE) {
  K <- ncol(P1)
  T <- qr.Q(qr(matrix(rnorm(N * K), N, K)))
  U <- T + matrix(rnorm(N * K, sd = noise), N, K)
  V <- T + matrix(rnorm(N * K, sd = noise), N, K)

  X <- list(
    D1 = U %*% t(P1) + matrix(rnorm(N * nrow(P1), sd = noise), N, nrow(P1))
  )
  if (isTRUE(include_d2)) {
    X$D2 <- U %*% t(P2) + matrix(rnorm(N * nrow(P2), sd = noise), N, nrow(P2))
  }

  Y <- list(
    E1 = V %*% t(Q1) + matrix(rnorm(N * nrow(Q1), sd = noise), N, nrow(Q1))
  )

  list(X = X, Y = Y)
}

.sim_aib_subject_correlated <- function(N, P1, P2, Q1, noise = 0.05, include_d2 = TRUE) {
  K <- ncol(P1)
  Zx <- matrix(rnorm(N * K), N, K)
  Zy <- matrix(rnorm(N * K), N, K)
  mix_x <- matrix(c(1, 0.9, 0.25, 1), K, K)
  mix_y <- matrix(c(1, 0.75, -0.2, 1), K, K)
  U <- Zx %*% mix_x + matrix(rnorm(N * K, sd = noise), N, K)
  V <- Zy %*% mix_y + matrix(rnorm(N * K, sd = noise), N, K)

  X <- list(
    D1 = U %*% t(P1) + matrix(rnorm(N * nrow(P1), sd = noise), N, nrow(P1))
  )
  if (isTRUE(include_d2)) {
    X$D2 <- U %*% t(P2) + matrix(rnorm(N * nrow(P2), sd = noise), N, nrow(P2))
  }

  Y <- list(
    E1 = V %*% t(Q1) + matrix(rnorm(N * nrow(Q1), sd = noise), N, nrow(Q1))
  )

  list(X = X, Y = Y)
}

.sim_aib_data <- function(subject_sizes = c(S1 = 14, S2 = 13, S3 = 12),
                          noise = 0.05,
                          seed = 1,
                          missing_d2 = character()) {
  set.seed(seed)
  K <- 2L
  P1 <- matrix(rnorm(5 * K), 5, K)
  P2 <- matrix(rnorm(4 * K), 4, K)
  Q1 <- matrix(rnorm(3 * K), 3, K)

  X <- list()
  Y <- list()
  for (subj in names(subject_sizes)) {
    sim_subj <- .sim_aib_subject(
      N = as.integer(subject_sizes[[subj]]),
      P1 = P1,
      P2 = P2,
      Q1 = Q1,
      noise = noise,
      include_d2 = !(subj %in% missing_d2)
    )
    X[[subj]] <- sim_subj$X
    Y[[subj]] <- sim_subj$Y
  }

  list(
    X = X,
    Y = Y,
    params = list(P1 = P1, P2 = P2, Q1 = Q1, K = K)
  )
}

.sim_aib_data_correlated <- function(subject_sizes = c(S1 = 14, S2 = 13, S3 = 12),
                                     noise = 0.05,
                                     seed = 1) {
  set.seed(seed)
  K <- 2L
  P1 <- matrix(rnorm(5 * K), 5, K)
  P2 <- matrix(rnorm(4 * K), 4, K)
  Q1 <- matrix(rnorm(3 * K), 3, K)

  X <- list()
  Y <- list()
  for (subj in names(subject_sizes)) {
    sim_subj <- .sim_aib_subject_correlated(
      N = as.integer(subject_sizes[[subj]]),
      P1 = P1,
      P2 = P2,
      Q1 = Q1,
      noise = noise
    )
    X[[subj]] <- sim_subj$X
    Y[[subj]] <- sim_subj$Y
  }

  list(
    X = X,
    Y = Y,
    params = list(P1 = P1, P2 = P2, Q1 = Q1, K = K)
  )
}

.aib_flatten_nested <- function(blocks_nested) {
  blocks_flat <- list()
  info_rows <- vector("list", 0L)
  for (subj in names(blocks_nested)) {
    for (dom in names(blocks_nested[[subj]])) {
      block_name <- paste(subj, dom, sep = "__")
      blocks_flat[[block_name]] <- blocks_nested[[subj]][[dom]]
      info_rows[[length(info_rows) + 1L]] <- data.frame(
        block = block_name,
        subject = subj,
        domain = dom,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    blocks = blocks_flat,
    block_info = do.call(rbind, info_rows)
  )
}

.aib_chain_laplacian <- function(n) {
  A <- matrix(0, n, n)
  if (n >= 2L) {
    idx <- seq_len(n - 1L)
    A[cbind(idx, idx + 1L)] <- 1
    A[cbind(idx + 1L, idx)] <- 1
  }
  diag(rowSums(A)) - A
}

.aib_chain_adjacency <- function(n) {
  A <- matrix(0, n, n)
  if (n >= 2L) {
    idx <- seq_len(n - 1L)
    A[cbind(idx, idx + 1L)] <- 1
    A[cbind(idx + 1L, idx)] <- 1
  }
  A
}

.aib_normalized_laplacian <- function(A) {
  deg <- rowSums(A)
  L <- diag(deg) - A
  D_half <- diag(ifelse(deg > 0, 1 / sqrt(deg), 0), nrow = length(deg))
  D_half %*% L %*% D_half
}

test_that("aligned_interbattery supports nested subject-domain input with missing domains", {
  sim <- .sim_aib_data(seed = 11, missing_d2 = "S2")

  fit <- aligned_interbattery(
    X = sim$X,
    Y = sim$Y,
    ncomp = 2,
    max_iter = 8
  )

  expect_s3_class(fit, "aligned_interbattery")
  expect_equal(fit$score_representation, "common_scores")
  expect_equal(fit$task, "bidirectional_prediction")
  expect_equal(fit$prediction_layer, "native_coupling")
  expect_equal(nrow(multivarious::scores(fit)), sum(vapply(fit$subject_sizes, identity, integer(1))))
  expect_equal(sort(unique(fit$x_block_info$subject)), c("S1", "S2", "S3"))
  expect_equal(sort(unique(fit$y_block_info$domain)), "E1")
  expect_true(all(c("D1", "D2") %in% names(fit$x_loadings)))
  expect_equal(names(fit$y_loadings), "E1")

  common_proj <- project(fit, sim$X$S1$D1[1:4, , drop = FALSE], block = "D1", from = "x", side = "common")
  common_pred <- predict(fit, sim$X$S1$D1[1:4, , drop = FALSE], block = "D1", from = "x", type = "common_scores")
  expect_equal(common_proj, common_pred)
  expect_equal(dim(common_proj), c(4, 2))
})

test_that("aligned_interbattery gives equivalent results for nested and flat block layouts", {
  sim <- .sim_aib_data(seed = 12)
  x_flat <- .aib_flatten_nested(sim$X)
  y_flat <- .aib_flatten_nested(sim$Y)

  fit_nested <- aligned_interbattery(
    X = sim$X,
    Y = sim$Y,
    ncomp = 2,
    max_iter = 8
  )
  fit_flat <- aligned_interbattery(
    X = x_flat$blocks,
    Y = y_flat$blocks,
    x_block_info = x_flat$block_info,
    y_block_info = y_flat$block_info,
    ncomp = 2,
    max_iter = 8
  )

  expect_lt(.aib_subspace_rel(multivarious::scores(fit_nested), multivarious::scores(fit_flat)), 1e-8)

  yhat_nested <- predict(
    fit_nested,
    sim$X$S1$D1[1:5, , drop = FALSE],
    block = "D1",
    from = "x",
    type = "prediction",
    target_block = "E1"
  )
  yhat_flat <- predict(
    fit_flat,
    sim$X$S1$D1[1:5, , drop = FALSE],
    block = "D1",
    from = "x",
    type = "prediction",
    target_block = "E1"
  )
  expect_equal(yhat_nested, yhat_flat, tolerance = 1e-8)
})

test_that("aligned_interbattery predicts across sides on synthetic data", {
  train <- .sim_aib_data(seed = 13, subject_sizes = c(S1 = 15, S2 = 14, S3 = 13), noise = 0.04)
  fit <- aligned_interbattery(
    X = train$X,
    Y = train$Y,
    ncomp = 2,
    max_iter = 10
  )

  new_subj <- .sim_aib_subject(
    N = 9,
    P1 = train$params$P1,
    P2 = train$params$P2,
    Q1 = train$params$Q1,
    noise = 0.04
  )

  yhat <- predict(
    fit,
    new_subj$X$D1,
    block = "D1",
    from = "x",
    type = "prediction",
    target_block = "E1"
  )
  xhat <- predict(
    fit,
    new_subj$Y$E1,
    block = "E1",
    from = "y",
    type = "prediction",
    target_block = "D1"
  )

  expect_equal(dim(yhat), dim(new_subj$Y$E1))
  expect_equal(dim(xhat), dim(new_subj$X$D1))
  expect_lt(mean((yhat - new_subj$Y$E1)^2), 0.05)
  expect_lt(mean((xhat - new_subj$X$D1)^2), 0.05)
})

test_that("aligned_interbattery supports subject-wise row graphs and orientation-invariant prediction", {
  sim <- .sim_aib_data(seed = 14, subject_sizes = c(S1 = 12, S2 = 11))
  row_graph <- lapply(sim$X, function(blocks) .aib_chain_laplacian(nrow(blocks$D1)))

  fit <- aligned_interbattery(
    X = sim$X,
    Y = sim$Y,
    ncomp = 2,
    row_graph = row_graph,
    row_graph_lambda = 0.2,
    max_iter = 8
  )
  fit_flip <- orient_components(fit, signs = c(-1, 1))

  expect_true(all(vapply(fit$row_graph_laplacian, function(L) is.matrix(L) || methods::is(L, "Matrix"), logical(1))))
  expect_true(all(is.finite(fit$objective_trace)))

  new_rows <- sim$X$S1$D1[1:4, , drop = FALSE]
  pred_ref <- predict(fit, new_rows, block = "D1", from = "x", type = "prediction", target_block = "E1")
  pred_flip <- predict(fit_flip, new_rows, block = "D1", from = "x", type = "prediction", target_block = "E1")
  expect_equal(pred_ref, pred_flip, tolerance = 1e-10)
})

test_that("aligned_interbattery projects and predicts from multi-block side bundles", {
  train <- .sim_aib_data(seed = 15, subject_sizes = c(S1 = 15, S2 = 14, S3 = 13), noise = 0.04)
  fit <- aligned_interbattery(
    X = train$X,
    Y = train$Y,
    ncomp = 2,
    max_iter = 10
  )

  new_subj <- .sim_aib_subject(
    N = 8,
    P1 = train$params$P1,
    P2 = train$params$P2,
    Q1 = train$params$Q1,
    noise = 0.04
  )

  scores_x <- project(
    fit,
    new_subj$X,
    from = "x",
    side = "x"
  )
  common_x <- project(
    fit,
    new_subj$X,
    from = "x",
    side = "common"
  )
  yhat <- predict(
    fit,
    new_subj$X,
    from = "x",
    type = "prediction",
    target_block = "E1"
  )
  recon_x <- predict(
    fit,
    new_subj$X,
    from = "x",
    type = "reconstruction"
  )

  expect_equal(dim(scores_x), c(8, 2))
  expect_equal(dim(common_x), c(8, 2))
  expect_equal(dim(yhat), dim(new_subj$Y$E1))
  expect_true(is.list(recon_x))
  expect_equal(sort(names(recon_x)), c("D1", "D2"))
  expect_lt(mean((yhat - new_subj$Y$E1)^2), 0.05)
})

test_that("aligned_interbattery supports joint completion from both sides", {
  train <- .sim_aib_data(seed = 16, subject_sizes = c(S1 = 16, S2 = 15, S3 = 14), noise = 0.04)
  fit <- aligned_interbattery(
    X = train$X,
    Y = train$Y,
    ncomp = 2,
    max_iter = 10
  )

  new_subj <- .sim_aib_subject(
    N = 9,
    P1 = train$params$P1,
    P2 = train$params$P2,
    Q1 = train$params$Q1,
    noise = 0.04
  )
  observed <- list(
    x = list(D1 = new_subj$X$D1),
    y = list(E1 = new_subj$Y$E1)
  )

  proj_both <- project(fit, observed, from = "both", side = "both")
  pred_d2 <- predict(
    fit,
    observed,
    from = "both",
    type = "prediction",
    target_block = "x:D2"
  )
  pred_all <- predict(
    fit,
    observed,
    from = "both",
    type = "prediction"
  )
  recon_e1 <- predict(
    fit,
    observed,
    from = "both",
    type = "reconstruction",
    target_block = "y:E1"
  )

  expect_equal(sort(names(proj_both)), c("common", "x", "y"))
  expect_equal(dim(proj_both$common), c(9, 2))
  expect_equal(dim(proj_both$x), c(9, 2))
  expect_equal(dim(proj_both$y), c(9, 2))
  expect_equal(dim(pred_d2), dim(new_subj$X$D2))
  expect_true(is.list(pred_all))
  expect_equal(sort(names(pred_all)), c("x", "y"))
  expect_true(all(c("D1", "D2") %in% names(pred_all$x)))
  expect_equal(names(pred_all$y), "E1")
  expect_equal(dim(recon_e1), dim(new_subj$Y$E1))
  expect_lt(mean((pred_d2 - new_subj$X$D2)^2), 0.05)
  expect_lt(mean((recon_e1 - new_subj$Y$E1)^2), 0.05)
})

test_that("aligned_interbattery decorrelate penalty reduces within-side score covariance", {
  sim <- .sim_aib_data_correlated(seed = 17, subject_sizes = c(S1 = 18, S2 = 17, S3 = 16), noise = 0.04)

  fit_base <- aligned_interbattery(
    X = sim$X,
    Y = sim$Y,
    ncomp = 2,
    decorrelate = 0,
    max_iter = 10
  )
  fit_decor <- aligned_interbattery(
    X = sim$X,
    Y = sim$Y,
    ncomp = 2,
    decorrelate = 10,
    max_iter = 10
  )

  off_x_base <- .aib_offdiag_cov_norm(fit_base$x_scores)
  off_x_decor <- .aib_offdiag_cov_norm(fit_decor$x_scores)
  off_y_base <- .aib_offdiag_cov_norm(fit_base$y_scores)
  off_y_decor <- .aib_offdiag_cov_norm(fit_decor$y_scores)

  expect_lt(off_y_decor, off_y_base)
  expect_lt(off_x_decor + off_y_decor, off_x_base + off_y_base)
  expect_equal(fit_decor$decorrelate, 10)
  expect_equal(fit_decor$decorrelate_type, "penalty")
})

test_that("aligned_interbattery supports whitening-based decorrelation with coherent orientation", {
  sim <- .sim_aib_data_correlated(seed = 117, subject_sizes = c(S1 = 18, S2 = 17, S3 = 16), noise = 0.04)

  fit_white <- aligned_interbattery(
    X = sim$X,
    Y = sim$Y,
    ncomp = 2,
    decorrelate = 0.8,
    decorrelate_type = "whiten",
    max_iter = 10
  )
  fit_white_flip <- orient_components(fit_white, signs = c(-1, 1))

  x_raw <- do.call(rbind, fit_white$x_scores)
  y_raw <- do.call(rbind, fit_white$y_scores)
  x_white <- x_raw %*% fit_white$x_whitener
  y_white <- y_raw %*% fit_white$y_whitener

  expect_equal(fit_white$decorrelate_type, "whiten")
  expect_equal(fit_white$prediction_layer, "native_coupling")
  expect_gt(max(abs(fit_white$x_whitener - diag(ncol(fit_white$x_whitener)))), 1e-6)
  expect_gt(max(abs(fit_white$y_whitener - diag(ncol(fit_white$y_whitener)))), 1e-6)
  expect_lt(.aib_cov_identity_distance(x_white), .aib_cov_identity_distance(x_raw))
  expect_lt(.aib_cov_identity_distance(y_white), .aib_cov_identity_distance(y_raw))

  new_rows <- sim$X$S1$D1[1:5, , drop = FALSE]
  pred_ref <- predict(fit_white, new_rows, block = "D1", from = "x", type = "prediction", target_block = "E1")
  pred_flip <- predict(fit_white_flip, new_rows, block = "D1", from = "x", type = "prediction", target_block = "E1")
  expect_equal(pred_ref, pred_flip, tolerance = 1e-10)
})

test_that("aligned_interbattery supports out-of-sample row maps for partial bundles", {
  train <- .sim_aib_data(seed = 18, subject_sizes = c(S1 = 15, S2 = 14, S3 = 13), noise = 0.04)
  fit <- aligned_interbattery(
    X = train$X,
    Y = train$Y,
    ncomp = 2,
    max_iter = 10
  )

  new_subj <- .sim_aib_subject(
    N = 9,
    P1 = train$params$P1,
    P2 = train$params$P2,
    Q1 = train$params$Q1,
    noise = 0.04
  )
  x_obs <- list(
    D1 = new_subj$X$D1[c(1, 3, 5, 7, 9), , drop = FALSE],
    D2 = new_subj$X$D2[c(2, 3, 4, 8), , drop = FALSE]
  )
  x_row_map_new <- list(
    D1 = c(1L, 3L, 5L, 7L, 9L),
    D2 = c(2L, 3L, 4L, 8L)
  )

  common_x <- project(
    fit,
    x_obs,
    from = "x",
    side = "common",
    x_row_map_new = x_row_map_new
  )
  yhat <- predict(
    fit,
    x_obs,
    from = "x",
    type = "prediction",
    target_block = "E1",
    x_row_map_new = x_row_map_new
  )
  recon_x <- predict(
    fit,
    x_obs,
    from = "x",
    type = "reconstruction",
    x_row_map_new = x_row_map_new
  )

  expect_equal(dim(common_x), c(9, 2))
  expect_equal(dim(yhat), dim(new_subj$Y$E1))
  expect_true(is.list(recon_x))
  expect_equal(sort(names(recon_x)), c("D1", "D2"))
  expect_equal(dim(recon_x$D1), dim(x_obs$D1))
  expect_equal(dim(recon_x$D2), dim(x_obs$D2))
  expect_lt(mean((yhat - new_subj$Y$E1)^2), 0.08)
})

test_that("aligned_interbattery reconstruction from both returns only observed sides", {
  train <- .sim_aib_data(seed = 19, subject_sizes = c(S1 = 15, S2 = 14, S3 = 13), noise = 0.04)
  fit <- aligned_interbattery(
    X = train$X,
    Y = train$Y,
    ncomp = 2,
    max_iter = 10
  )

  new_subj <- .sim_aib_subject(
    N = 7,
    P1 = train$params$P1,
    P2 = train$params$P2,
    Q1 = train$params$Q1,
    noise = 0.04
  )
  observed <- list(x = list(D1 = new_subj$X$D1), y = NULL)

  recon <- predict(
    fit,
    observed,
    from = "both",
    type = "reconstruction"
  )

  expect_equal(names(recon), "x")
  expect_equal(dim(recon$x), dim(new_subj$X$D1))
  expect_error(
    predict(
      fit,
      observed,
      from = "both",
      type = "reconstruction",
      target_block = "y:E1"
    ),
    "requires observed y-side blocks"
  )
})

test_that("aligned_interbattery accepts precomputed normalized Laplacians", {
  sim <- .sim_aib_data(seed = 20, subject_sizes = c(S1 = 11, S2 = 10))
  row_graph_adj <- lapply(sim$X, function(blocks) .aib_chain_adjacency(nrow(blocks$D1)))
  row_graph_norm <- lapply(row_graph_adj, .aib_normalized_laplacian)

  fit_adj <- aligned_interbattery(
    X = sim$X,
    Y = sim$Y,
    ncomp = 2,
    row_graph = row_graph_adj,
    row_graph_form = "normalized_laplacian",
    row_graph_lambda = 0.2,
    max_iter = 6
  )
  fit_norm <- aligned_interbattery(
    X = sim$X,
    Y = sim$Y,
    ncomp = 2,
    row_graph = row_graph_norm,
    row_graph_form = "normalized_laplacian",
    row_graph_lambda = 0.2,
    max_iter = 6
  )

  expect_equal(as.matrix(fit_adj$row_graph_laplacian$S1), as.matrix(fit_norm$row_graph_laplacian$S1), tolerance = 1e-10)
  expect_equal(as.matrix(fit_adj$row_graph_laplacian$S2), as.matrix(fit_norm$row_graph_laplacian$S2), tolerance = 1e-10)
})
