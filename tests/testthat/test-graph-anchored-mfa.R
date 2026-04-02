library(testthat)
library(muscal)

.score_subspace_rel <- function(S1, S2) {
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
}

test_that("graph_anchored_mfa reduces to anchored_mfa when graph penalty is off", {
  set.seed(11)
  N <- 40
  Y <- matrix(rnorm(N * 4), N, 4)
  X1 <- matrix(rnorm(18 * 6), 18, 6)
  X2 <- matrix(rnorm(15 * 5), 15, 5)
  idx1 <- sample.int(N, nrow(X1), replace = FALSE)
  idx2 <- sample.int(N, nrow(X2), replace = FALSE)

  fit_anchor <- anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = idx1, X2 = idx2),
    ncomp = 2,
    ridge = 1e-8,
    max_iter = 100,
    tol = 1e-10
  )

  fit_graph <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = idx1, X2 = idx2),
    ncomp = 2,
    graph_lambda = 0,
    ridge = 1e-8,
    max_iter = 100,
    tol = 1e-10
  )

  expect_s3_class(fit_graph, "graph_anchored_mfa")
  expect_s3_class(fit_graph, "anchored_mfa")
  expect_s3_class(fit_graph, "linked_mfa")

  S1 <- multivarious::scores(fit_anchor)
  S2 <- multivarious::scores(fit_graph)
  rel <- .score_subspace_rel(S1, S2)
  expect_lt(rel, 1e-10)
})

test_that("graph_anchored_mfa supports nested subject-domain input with missing domains", {
  set.seed(12)
  N <- 50
  Y <- matrix(rnorm(N * 5), N, 5)

  X <- list(
    S1 = list(
      D1 = matrix(rnorm(16 * 6), 16, 6),
      D2 = matrix(rnorm(12 * 4), 12, 4)
    ),
    S2 = list(
      D1 = matrix(rnorm(14 * 6), 14, 6)
    ),
    S3 = list(
      D1 = matrix(rnorm(15 * 6), 15, 6),
      D2 = matrix(rnorm(11 * 4), 11, 4)
    )
  )

  row_index <- list(
    S1 = list(
      D1 = sample.int(N, nrow(X$S1$D1), replace = FALSE),
      D2 = sample.int(N, nrow(X$S1$D2), replace = FALSE)
    ),
    S2 = list(
      D1 = sample.int(N, nrow(X$S2$D1), replace = FALSE)
    ),
    S3 = list(
      D1 = sample.int(N, nrow(X$S3$D1), replace = FALSE),
      D2 = sample.int(N, nrow(X$S3$D2), replace = FALSE)
    )
  )

  fit <- graph_anchored_mfa(
    Y = Y,
    X = X,
    row_index = row_index,
    ncomp = 2,
    graph_lambda = 0
  )

  expect_s3_class(fit, "graph_anchored_mfa")
  expect_equal(nrow(multivarious::scores(fit)), N)
  expect_equal(nrow(fit$block_info), 5)
  expect_equal(sort(unique(fit$block_info$subject)), c("S1", "S2", "S3"))
  expect_equal(sort(na.omit(unique(fit$block_info$domain))), c("D1", "D2"))
  expect_true(all(c("S1__D1", "S1__D2", "S2__D1", "S3__D1", "S3__D2") %in% names(fit$V_list)))
})

test_that("graph penalty shrinks linked loadings across blocks", {
  set.seed(13)
  N <- 60
  Y <- matrix(rnorm(N * 4), N, 4)

  shared <- paste0("f", seq_len(4))
  cn1 <- c(shared, "u1")
  cn2 <- c(shared, "u2")

  X1 <- matrix(rnorm(N * length(cn1)), N, length(cn1))
  X2 <- matrix(rnorm(N * length(cn2)), N, length(cn2))
  colnames(X1) <- cn1
  colnames(X2) <- cn2

  idx <- list(X1 = seq_len(N), X2 = seq_len(N))

  fit0 <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = idx,
    ncomp = 2,
    feature_graph = "colnames",
    graph_lambda = 0
  )

  fit1 <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = idx,
    ncomp = 2,
    feature_graph = "colnames",
    graph_lambda = 10
  )

  V0_1 <- fit0$V_list$X1[shared, , drop = FALSE]
  V0_2 <- fit0$V_list$X2[shared, , drop = FALSE]
  V1_1 <- fit1$V_list$X1[shared, , drop = FALSE]
  V1_2 <- fit1$V_list$X2[shared, , drop = FALSE]

  d0 <- rowSums((V0_1 - V0_2)^2)
  d1 <- rowSums((V1_1 - V1_2)^2)

  expect_lt(mean(d1), mean(d0))
})

test_that("graph_anchored_mfa is invariant to equivalent graph representations", {
  set.seed(15)
  N <- 50
  Y <- matrix(rnorm(N * 4), N, 4)

  shared <- paste0("f", seq_len(4))
  cn1 <- c(shared, "u1")
  cn2 <- c(shared, "u2")

  X1 <- matrix(rnorm(N * length(cn1)), N, length(cn1))
  X2 <- matrix(rnorm(N * length(cn2)), N, length(cn2))
  colnames(X1) <- cn1
  colnames(X2) <- cn2

  idx <- list(X1 = seq_len(N), X2 = seq_len(N))
  edge_graph <- data.frame(
    block1 = rep("X1", length(shared)),
    feature1 = shared,
    block2 = rep("X2", length(shared)),
    feature2 = shared,
    weight = 1,
    stringsAsFactors = FALSE
  )

  adj <- Matrix::sparseMatrix(
    i = c(1:4, 6:9),
    j = c(6:9, 1:4),
    x = 1,
    dims = c(10, 10)
  )

  fit_names <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = idx,
    ncomp = 2,
    feature_graph = "colnames",
    graph_lambda = 5
  )

  fit_edges <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = idx,
    ncomp = 2,
    feature_graph = edge_graph,
    graph_lambda = 5
  )

  fit_adj <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = idx,
    ncomp = 2,
    feature_graph = adj,
    graph_form = "adjacency",
    graph_lambda = 5
  )

  expect_lt(.score_subspace_rel(multivarious::scores(fit_names), multivarious::scores(fit_edges)), 1e-10)
  expect_lt(.score_subspace_rel(multivarious::scores(fit_names), multivarious::scores(fit_adj)), 1e-10)
})

test_that("graph_anchored_mfa gives equivalent results for nested and flat block layouts", {
  set.seed(16)
  N <- 55
  Y <- matrix(rnorm(N * 4), N, 4)

  X_nested <- list(
    S1 = list(
      D1 = matrix(rnorm(14 * 5), 14, 5),
      D2 = matrix(rnorm(13 * 4), 13, 4)
    ),
    S2 = list(
      D1 = matrix(rnorm(15 * 5), 15, 5)
    ),
    S3 = list(
      D1 = matrix(rnorm(16 * 5), 16, 5),
      D2 = matrix(rnorm(12 * 4), 12, 4)
    )
  )

  row_nested <- list(
    S1 = list(
      D1 = sample.int(N, nrow(X_nested$S1$D1), replace = FALSE),
      D2 = sample.int(N, nrow(X_nested$S1$D2), replace = FALSE)
    ),
    S2 = list(
      D1 = sample.int(N, nrow(X_nested$S2$D1), replace = FALSE)
    ),
    S3 = list(
      D1 = sample.int(N, nrow(X_nested$S3$D1), replace = FALSE),
      D2 = sample.int(N, nrow(X_nested$S3$D2), replace = FALSE)
    )
  )

  X_flat <- list(
    S1__D1 = X_nested$S1$D1,
    S1__D2 = X_nested$S1$D2,
    S2__D1 = X_nested$S2$D1,
    S3__D1 = X_nested$S3$D1,
    S3__D2 = X_nested$S3$D2
  )
  row_flat <- list(
    S1__D1 = row_nested$S1$D1,
    S1__D2 = row_nested$S1$D2,
    S2__D1 = row_nested$S2$D1,
    S3__D1 = row_nested$S3$D1,
    S3__D2 = row_nested$S3$D2
  )

  fit_nested <- graph_anchored_mfa(Y = Y, X = X_nested, row_index = row_nested, ncomp = 2, graph_lambda = 0)
  fit_flat <- graph_anchored_mfa(Y = Y, X = X_flat, row_index = row_flat, ncomp = 2, graph_lambda = 0)

  expect_lt(.score_subspace_rel(multivarious::scores(fit_nested), multivarious::scores(fit_flat)), 1e-10)
})

test_that("graph_anchored_mfa objective trace is non-increasing", {
  set.seed(17)
  N <- 50
  Y <- matrix(rnorm(N * 4), N, 4)
  X1 <- matrix(rnorm(N * 5), N, 5)
  X2 <- matrix(rnorm(N * 5), N, 5)
  colnames(X1) <- c(paste0("f", 1:4), "u1")
  colnames(X2) <- c(paste0("f", 1:4), "u2")

  fit <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = seq_len(N), X2 = seq_len(N)),
    ncomp = 2,
    feature_graph = "colnames",
    graph_lambda = 3,
    max_iter = 20
  )

  diffs <- diff(fit$objective_trace)
  expect_true(all(diffs <= 1e-8))
})

test_that("graph_anchored_mfa returns finite outputs for p > n and duplicated-column blocks", {
  set.seed(18)
  N <- 30
  Y <- matrix(rnorm(N * 3), N, 3)

  X1_base <- matrix(rnorm(10 * 4), 10, 4)
  X2_base <- matrix(rnorm(12 * 4), 12, 4)
  X1 <- cbind(X1_base, X1_base, X1_base[, 1, drop = FALSE])
  X2 <- cbind(X2_base, X2_base, X2_base[, 1, drop = FALSE])
  colnames(X1) <- c(paste0("f", 1:4), paste0("f", 1:4), "f1_dup")
  colnames(X2) <- c(paste0("f", 1:4), paste0("f", 1:4), "f1_dup")

  fit <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(
      X1 = sample.int(N, nrow(X1), replace = TRUE),
      X2 = sample.int(N, nrow(X2), replace = TRUE)
    ),
    ncomp = 2,
    feature_graph = "colnames",
    graph_lambda = 2,
    ridge = 1e-6
  )

  expect_true(all(is.finite(multivarious::scores(fit))))
  expect_true(all(vapply(fit$V_list, function(v) all(is.finite(v)), logical(1))))

  new_rows <- X1[1:3, , drop = FALSE]
  expect_true(all(is.finite(project(fit, new_rows, block = "X1"))))
  expect_true(all(is.finite(predict(fit, new_rows, block = "X1"))))
})

test_that("graph_anchored_mfa validates malformed feature-graph inputs", {
  set.seed(19)
  N <- 20
  Y <- matrix(rnorm(N * 3), N, 3)
  X1 <- matrix(rnorm(N * 4), N, 4)
  X2 <- matrix(rnorm(N * 5), N, 5)

  bad_graph <- Matrix::Diagonal(3)
  bad_edges <- data.frame(
    block1 = "X1",
    feature1 = 1,
    block2 = "X2",
    feature2 = 1,
    weight = -1
  )

  expect_error(
    graph_anchored_mfa(
      Y = Y,
      X = list(X1 = X1, X2 = X2),
      row_index = list(X1 = seq_len(N), X2 = seq_len(N)),
      ncomp = 2,
      feature_graph = bad_graph,
      graph_form = "adjacency",
      graph_lambda = 1
    )
  )

  expect_error(
    graph_anchored_mfa(
      Y = Y,
      X = list(X1 = X1, X2 = X2),
      row_index = list(X1 = seq_len(N), X2 = seq_len(N)),
      ncomp = 2,
      feature_graph = bad_edges,
      graph_lambda = 1
    )
  )
})

test_that("graph_anchored_mfa projects scores and predicts Y for new rows", {
  set.seed(14)
  N <- 45
  Y <- matrix(rnorm(N * 3), N, 3)
  X1 <- matrix(rnorm(20 * 6), 20, 6)
  X2 <- matrix(rnorm(18 * 5), 18, 5)
  idx1 <- sample.int(N, nrow(X1), replace = FALSE)
  idx2 <- sample.int(N, nrow(X2), replace = FALSE)

  fit <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = idx1, X2 = idx2),
    ncomp = 2,
    graph_lambda = 0
  )

  new_rows <- X1[1:4, , drop = FALSE]
  scores_new <- project(fit, new_rows, block = "X1")
  yhat <- predict(fit, new_rows, block = "X1")
  scores_from_predict <- predict(fit, new_rows, block = "X1", type = "scores")
  xhat <- predict(fit, new_rows, block = "X1", type = "reconstruction")

  expect_equal(fit$fit_spec$method, "graph_anchored_mfa")
  expect_equal(fit$task, "response_prediction")
  expect_setequal(fit$oos_types, c("response", "scores", "reconstruction"))
  expect_equal(dim(scores_new), c(4, 2))
  expect_equal(dim(scores_from_predict), c(4, 2))
  expect_equal(dim(yhat), c(4, ncol(Y)))
  expect_equal(dim(xhat), dim(new_rows))
  expect_equal(unname(scores_from_predict), unname(scores_new), tolerance = 1e-10)
})
