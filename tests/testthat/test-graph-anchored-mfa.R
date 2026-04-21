library(testthat)
library(muscal)

.projected_grad_norm <- function(S, G) {
  PG <- G - S %*% ((crossprod(S, G) + t(crossprod(S, G))) / 2)
  sqrt(sum(PG^2))
}

.score_subspace_rel <- function(S1, S2) {
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
}

.graph_anchor_exact_fixture <- function() {
  set.seed(1201)

  N <- 72
  k <- 2
  S <- scale(matrix(rnorm(N * k), N, k), center = TRUE, scale = FALSE)
  B <- matrix(
    c(1.2, -0.4,
      0.8, 0.6,
      -0.7, 0.9),
    nrow = 3,
    byrow = TRUE
  )
  V_shared <- matrix(
    c(1.0, 0.0,
      0.8, 0.2,
      -0.6, 1.1,
      0.4, -0.9),
    nrow = 4,
    byrow = TRUE
  )
  V1_unique <- matrix(
    c(0.7, 0.3,
      -0.5, 0.6),
    nrow = 2,
    byrow = TRUE
  )
  V2_unique <- matrix(
    c(-0.4, 0.8,
      0.9, -0.2),
    nrow = 2,
    byrow = TRUE
  )

  Y <- S %*% t(B)
  X1 <- cbind(S %*% t(V_shared), S %*% t(V1_unique))
  X2 <- cbind(S %*% t(V_shared), S %*% t(V2_unique))
  colnames(X1) <- c(paste0("f", 1:4), "u1", "u2")
  colnames(X2) <- c(paste0("f", 1:4), "v1", "v2")

  edge_graph <- data.frame(
    block1 = rep("X1", 4),
    feature1 = paste0("f", 1:4),
    block2 = rep("X2", 4),
    feature2 = paste0("f", 1:4),
    weight = 1,
    stringsAsFactors = FALSE
  )

  list(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = seq_len(N), X2 = seq_len(N)),
    S = S,
    feature_graph = edge_graph,
    shared = paste0("f", 1:4),
    ncomp = k
  )
}

.graph_anchor_sparse_borrowing_fixture <- function() {
  set.seed(1202)

  N <- 90
  k <- 2
  S <- scale(matrix(rnorm(N * k), N, k), center = TRUE, scale = FALSE)
  B <- matrix(
    c(1.2, -0.4,
      0.8, 0.6,
      -0.7, 0.9),
    nrow = 3,
    byrow = TRUE
  )
  V_shared <- matrix(
    c(1.0, 0.0,
      0.8, 0.2,
      -0.6, 1.1,
      0.4, -0.9),
    nrow = 4,
    byrow = TRUE
  )
  V1_unique <- matrix(
    c(0.7, 0.3,
      -0.5, 0.6),
    nrow = 2,
    byrow = TRUE
  )
  V2_unique <- matrix(
    c(-0.4, 0.8,
      0.9, -0.2),
    nrow = 2,
    byrow = TRUE
  )

  Y <- S %*% t(B) + matrix(rnorm(N * 3, sd = 0.05), N, 3)
  X1_full <- cbind(S %*% t(V_shared), S %*% t(V1_unique)) +
    matrix(rnorm(N * 6, sd = 0.30), N, 6)
  X2 <- cbind(S %*% t(V_shared), S %*% t(V2_unique)) +
    matrix(rnorm(N * 6, sd = 0.05), N, 6)
  idx1 <- sort(sample.int(N, 10))
  X1 <- X1_full[idx1, , drop = FALSE]

  colnames(X1) <- c(paste0("f", 1:4), "u1", "u2")
  colnames(X2) <- c(paste0("f", 1:4), "v1", "v2")

  wrong_graph <- data.frame(
    block1 = rep("X1", 4),
    feature1 = paste0("f", 1:4),
    block2 = rep("X2", 4),
    feature2 = paste0("f", 4:1),
    weight = 1,
    stringsAsFactors = FALSE
  )

  list(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    X1_full = X1_full,
    row_index = list(X1 = idx1, X2 = seq_len(N)),
    S = S,
    shared = paste0("f", 1:4),
    wrong_graph = wrong_graph,
    ncomp = k
  )
}

.graph_anchor_score_graph_fixture <- function() {
  set.seed(1203)

  n_pairs <- 12
  k <- 2
  S_pair <- scale(matrix(rnorm(n_pairs * k), n_pairs, k), center = TRUE, scale = FALSE)
  S_true <- S_pair[rep(seq_len(n_pairs), each = 2), , drop = FALSE]
  N <- nrow(S_true)

  B <- matrix(
    c(1.1, -0.3,
      0.7, 0.8,
      -0.5, 1.0),
    nrow = 3,
    byrow = TRUE
  )
  V1 <- matrix(
    c(0.9, 0.1,
      -0.6, 0.8,
      0.5, -0.7,
      0.4, 0.6),
    nrow = 4,
    byrow = TRUE
  )
  V2 <- matrix(
    c(0.8, -0.2,
      -0.4, 0.9,
      0.7, -0.5,
      -0.3, 0.8),
    nrow = 4,
    byrow = TRUE
  )

  Y <- S_true %*% t(B) + matrix(rnorm(N * 3, sd = 0.60), N, 3)
  X1_full <- S_true %*% t(V1) + matrix(rnorm(N * 4, sd = 0.03), N, 4)
  X2_full <- S_true %*% t(V2) + matrix(rnorm(N * 4, sd = 0.03), N, 4)

  idx1 <- seq(1, N, by = 2)
  idx2 <- seq(2, N, by = 2)
  X1 <- X1_full[idx1, , drop = FALSE]
  X2 <- X2_full[idx2, , drop = FALSE]
  colnames(X1) <- paste0("x1_", seq_len(ncol(X1)))
  colnames(X2) <- paste0("x2_", seq_len(ncol(X2)))

  edge_graph <- data.frame(
    row1 = idx1,
    row2 = idx2,
    weight = 1,
    stringsAsFactors = FALSE
  )

  adj <- Matrix::sparseMatrix(
    i = c(idx1, idx2),
    j = c(idx2, idx1),
    x = 1,
    dims = c(N, N)
  )

  list(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = idx1, X2 = idx2),
    X1_full = X1_full,
    S_true = S_true,
    pair_index = split(seq_len(N), rep(seq_len(n_pairs), each = 2)),
    score_graph_edges = edge_graph,
    score_graph_adj = adj,
    ncomp = k
  )
}

test_that("graph_anchored_mfa exactly recovers a noiseless low-rank oracle with a consistent graph", {
  sim <- .graph_anchor_exact_fixture()

  fit <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    feature_graph = sim$feature_graph,
    graph_lambda = 6,
    max_iter = 100,
    tol = 1e-10,
    ridge = 1e-10
  )

  expect_lt(.score_subspace_rel(multivarious::scores(fit), sim$S), 1e-10)
  expect_lt(mean((predict(fit, sim$X$X1, block = "X1", type = "response", preprocess = FALSE) - sim$Y)^2), 1e-20)
  expect_lt(mean((predict(fit, sim$X$X1, block = "X1", type = "reconstruction", preprocess = FALSE) - sim$X$X1)^2), 1e-20)

  linked_gap <- rowSums((fit$V_list$X1[sim$shared, , drop = FALSE] - fit$V_list$X2[sim$shared, , drop = FALSE])^2)
  expect_lt(max(linked_gap), 1e-18)
})

test_that("graph_anchored_mfa borrows through the correct graph under sparse noisy coverage", {
  sim <- .graph_anchor_sparse_borrowing_fixture()

  fit_none <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    graph_lambda = 0,
    max_iter = 120,
    tol = 1e-8,
    ridge = 1e-8
  )

  fit_correct <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    feature_graph = "colnames",
    graph_lambda = 10,
    max_iter = 120,
    tol = 1e-8,
    ridge = 1e-8
  )

  fit_wrong <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    feature_graph = sim$wrong_graph,
    graph_lambda = 10,
    max_iter = 120,
    tol = 1e-8,
    ridge = 1e-8
  )

  mse_none <- mean((predict(fit_none, sim$X1_full, block = "X1", type = "response", preprocess = FALSE) - sim$Y)^2)
  mse_correct <- mean((predict(fit_correct, sim$X1_full, block = "X1", type = "response", preprocess = FALSE) - sim$Y)^2)
  mse_wrong <- mean((predict(fit_wrong, sim$X1_full, block = "X1", type = "response", preprocess = FALSE) - sim$Y)^2)

  recon_none <- mean((predict(fit_none, sim$X1_full, block = "X1", type = "reconstruction", preprocess = FALSE) - sim$X1_full)^2)
  recon_correct <- mean((predict(fit_correct, sim$X1_full, block = "X1", type = "reconstruction", preprocess = FALSE) - sim$X1_full)^2)

  linked_gap_none <- mean(rowSums((fit_none$V_list$X1[sim$shared, , drop = FALSE] - fit_none$V_list$X2[sim$shared, , drop = FALSE])^2))
  linked_gap_correct <- mean(rowSums((fit_correct$V_list$X1[sim$shared, , drop = FALSE] - fit_correct$V_list$X2[sim$shared, , drop = FALSE])^2))

  expect_lt(mse_correct, mse_none)
  expect_lt(mse_correct, 0.2 * mse_wrong)
  expect_lt(recon_correct, recon_none)
  expect_lt(linked_gap_correct, 1e-3 * linked_gap_none)
})

test_that("score-graph penalty is off-switch equivalent to the baseline graph_anchored_mfa fit", {
  sim <- .graph_anchor_score_graph_fixture()

  fit_base <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    graph_lambda = 0,
    score_graph = sim$score_graph_edges,
    score_graph_lambda = 0,
    max_iter = 80,
    tol = 1e-10,
    ridge = 1e-10
  )

  fit_plain <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    graph_lambda = 0,
    max_iter = 80,
    tol = 1e-10,
    ridge = 1e-10
  )

  expect_lt(.score_subspace_rel(multivarious::scores(fit_base), multivarious::scores(fit_plain)), 1e-10)
})

test_that("score graph pulls similar anchor rows together and improves score recovery under split coverage", {
  sim <- .graph_anchor_score_graph_fixture()

  fit_none <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    graph_lambda = 0,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  fit_score <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    graph_lambda = 0,
    score_graph = sim$score_graph_edges,
    score_graph_lambda = 8,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  pair_gap <- function(S) {
    mean(vapply(sim$pair_index, function(idx) sum((S[idx[1], ] - S[idx[2], ])^2), numeric(1)))
  }

  gap_none <- pair_gap(multivarious::scores(fit_none))
  gap_score <- pair_gap(multivarious::scores(fit_score))
  rel_none <- .score_subspace_rel(multivarious::scores(fit_none), sim$S_true)
  rel_score <- .score_subspace_rel(multivarious::scores(fit_score), sim$S_true)

  expect_lt(gap_score, 0.7 * gap_none)
  expect_lt(rel_score, rel_none)
})

test_that("larger sparse score-graph weights induce stronger score shrinkage", {
  sim <- .graph_anchor_score_graph_fixture()

  weighted_adj <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4),
    j = c(2, 1, 4, 3),
    x = c(10, 10, 1, 1),
    dims = c(nrow(sim$Y), nrow(sim$Y))
  )

  fit_none <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    score_graph = weighted_adj,
    score_graph_form = "adjacency",
    score_graph_lambda = 0,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  fit_weighted <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    score_graph = weighted_adj,
    score_graph_form = "adjacency",
    score_graph_lambda = 10,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  gap <- function(S, i, j) sum((S[i, ] - S[j, ])^2)

  high_ratio <- gap(fit_weighted$s, 1, 2) / gap(fit_none$s, 1, 2)
  low_ratio <- gap(fit_weighted$s, 3, 4) / gap(fit_none$s, 3, 4)

  expect_lt(high_ratio, low_ratio)
  expect_lt(high_ratio, 0.1)
})

test_that("score-graph representations are equivalent across edge and adjacency inputs", {
  sim <- .graph_anchor_score_graph_fixture()

  fit_edges <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    score_graph = sim$score_graph_edges,
    score_graph_lambda = 8,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  fit_adj <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    score_graph = sim$score_graph_adj,
    score_graph_form = "adjacency",
    score_graph_lambda = 8,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  expect_lt(.score_subspace_rel(multivarious::scores(fit_edges), multivarious::scores(fit_adj)), 1e-10)
})

test_that("graph_anchored_mfa accepts symmetric sparse adjacency matrices from Matrix::Matrix()", {
  sim <- .graph_anchor_score_graph_fixture()
  feature_edges <- data.frame(
    block1 = "X1",
    feature1 = colnames(sim$X$X1),
    block2 = "X2",
    feature2 = colnames(sim$X$X2),
    weight = 1,
    stringsAsFactors = FALSE
  )
  feature_adj <- muscal:::.gamfa_graph_from_edges(feature_edges, sim$X)
  feature_ds <- Matrix::Matrix(as.matrix(feature_adj), sparse = TRUE)
  score_ds <- Matrix::Matrix(as.matrix(sim$score_graph_adj), sparse = TRUE)

  expect_s4_class(feature_ds, "dsCMatrix")
  expect_s4_class(score_ds, "dsCMatrix")

  fit_dg <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    feature_graph = feature_adj,
    graph_form = "adjacency",
    graph_lambda = 3,
    score_graph = sim$score_graph_adj,
    score_graph_form = "adjacency",
    score_graph_lambda = 8,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  fit_ds <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    feature_graph = feature_ds,
    graph_form = "adjacency",
    graph_lambda = 3,
    score_graph = score_ds,
    score_graph_form = "adjacency",
    score_graph_lambda = 8,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  expect_lt(.score_subspace_rel(multivarious::scores(fit_dg), multivarious::scores(fit_ds)), 1e-10)
  expect_equal(unname(tail(fit_ds$objective_trace, 1)), unname(tail(fit_dg$objective_trace, 1)), tolerance = 1e-8)
})

test_that("precomputed graph Laplacians match equivalent adjacency-matrix inputs", {
  sim <- .graph_anchor_score_graph_fixture()
  feature_edges <- data.frame(
    block1 = "X1",
    feature1 = colnames(sim$X$X1),
    block2 = "X2",
    feature2 = colnames(sim$X$X2),
    weight = 1,
    stringsAsFactors = FALSE
  )
  feature_adj <- muscal:::.gamfa_graph_from_edges(feature_edges, sim$X)
  feature_lap <- muscal:::.gamfa_laplacian_from_adjacency(feature_adj)
  score_lap <- muscal:::.gamfa_laplacian_from_adjacency(sim$score_graph_adj)

  fit_adj <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    feature_graph = feature_adj,
    graph_form = "adjacency",
    graph_lambda = 3,
    score_graph = sim$score_graph_adj,
    score_graph_form = "adjacency",
    score_graph_lambda = 8,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  fit_lap <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    feature_graph = feature_lap,
    graph_form = "laplacian",
    graph_lambda = 3,
    score_graph = score_lap,
    score_graph_form = "laplacian",
    score_graph_lambda = 8,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  expect_lt(.score_subspace_rel(multivarious::scores(fit_adj), multivarious::scores(fit_lap)), 1e-10)
})

test_that("normalized Laplacian matrices match normalized graph construction from edge inputs", {
  sim <- .graph_anchor_score_graph_fixture()
  feature_edges <- data.frame(
    block1 = "X1",
    feature1 = colnames(sim$X$X1),
    block2 = "X2",
    feature2 = colnames(sim$X$X2),
    weight = 1,
    stringsAsFactors = FALSE
  )
  feature_norm_lap <- muscal:::.gamfa_laplacian_from_adjacency(
    muscal:::.gamfa_graph_from_edges(feature_edges, sim$X),
    normalized = TRUE
  )
  score_norm_lap <- muscal:::.gamfa_laplacian_from_adjacency(sim$score_graph_adj, normalized = TRUE)

  fit_edges <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    feature_graph = feature_edges,
    graph_form = "normalized_laplacian",
    graph_lambda = 3,
    score_graph = sim$score_graph_edges,
    score_graph_form = "normalized_laplacian",
    score_graph_lambda = 8,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  fit_norm <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    feature_graph = feature_norm_lap,
    graph_form = "normalized_laplacian",
    graph_lambda = 3,
    score_graph = score_norm_lap,
    score_graph_form = "normalized_laplacian",
    score_graph_lambda = 8,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  expect_lt(.score_subspace_rel(multivarious::scores(fit_edges), multivarious::scores(fit_norm)), 1e-10)
})

test_that("score_graph='knn' builds a weighted sparse similarity graph and returns finite outputs", {
  sim <- .graph_anchor_score_graph_fixture()

  fit_knn <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    score_graph = "knn",
    score_graph_k = 1,
    score_graph_weight_mode = "heat",
    score_graph_sigma = 0.75,
    score_graph_lambda = 8,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  pair_gap <- mean(vapply(sim$pair_index, function(idx) sum((fit_knn$s[idx[1], ] - fit_knn$s[idx[2], ])^2), numeric(1)))

  expect_true(all(is.finite(multivarious::scores(fit_knn))))
  expect_true(all(is.finite(fit_knn$B)))
  expect_true(all(is.finite(fit_knn$score_graph_laplacian)))
  expect_s4_class(fit_knn$score_graph_adjacency, "dgCMatrix")
  expect_gt(length(unique(round(fit_knn$score_graph_adjacency@x, 8))), 1)
  expect_lt(pair_gap, 0.1)
})

test_that("built-in weighted kNN score graph matches adjoin heat-kNN adjacency", {
  skip_if_not_installed("adjoin")
  sim <- .graph_anchor_score_graph_fixture()
  sigma <- 0.75

  fit_builtin <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    score_graph = "knn",
    score_graph_k = 3,
    score_graph_weight_mode = "heat",
    score_graph_sigma = sigma,
    score_graph_lambda = 8,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  ng <- adjoin::graph_weights(
    sim$Y,
    neighbor_mode = "knn",
    k = 3,
    weight_mode = "heat",
    sigma = sigma
  )
  A_adjoin <- adjoin::adjacency(ng)

  fit_adjoin <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    score_graph = A_adjoin,
    score_graph_form = "adjacency",
    score_graph_lambda = 8,
    max_iter = 100,
    tol = 1e-8,
    ridge = 1e-8
  )

  A_builtin <- fit_builtin$score_graph_adjacency
  adj_rel <- Matrix::norm(A_builtin - A_adjoin, "F") / (Matrix::norm(A_adjoin, "F") + 1e-12)

  expect_s4_class(A_builtin, "dgCMatrix")
  expect_lt(adj_rel, 1e-12)
  expect_lt(.score_subspace_rel(multivarious::scores(fit_builtin), multivarious::scores(fit_adjoin)), 1e-12)
})

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

test_that("graph_anchored_mfa supports orthonormal score constraint", {
  set.seed(15)
  Y <- matrix(rnorm(40 * 4), 40, 4)
  X1 <- matrix(rnorm(18 * 6), 18, 6)
  X2 <- matrix(rnorm(16 * 5), 16, 5)
  idx <- list(X1 = sample.int(nrow(Y), nrow(X1), replace = TRUE),
              X2 = sample.int(nrow(Y), nrow(X2), replace = TRUE))

  fit <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = idx,
    ncomp = 2,
    score_constraint = "orthonormal",
    score_graph = "knn",
    score_graph_lambda = 0.25,
    max_iter = 15
  )

  expect_equal(crossprod(fit$s), diag(2), tolerance = 1e-5)
  expect_true(all(is.finite(fit$objective_trace)))
})

test_that("graph_anchored_mfa orthonormal path has a small projected gradient residual", {
  sim <- .graph_anchor_score_graph_fixture()

  fit <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    score_constraint = "orthonormal",
    score_graph = sim$score_graph_adj,
    score_graph_form = "adjacency",
    score_graph_lambda = 2,
    ridge = 1e-6,
    max_iter = 40,
    tol = 1e-8
  )

  grad <- muscal:::.gamfa_score_gradient(
    S = fit$s,
    Y = sim$Y,
    B = fit$B,
    X_list = sim$X,
    V_list = fit$V_list,
    row_index = fit$row_index,
    alpha_y = fit$alpha_blocks[[1]],
    alpha_blocks = unname(fit$alpha_blocks[-1]),
    score_graph_laplacian = fit$score_graph_laplacian,
    score_graph_lambda = fit$score_graph_lambda,
    ridge = fit$ridge
  )
  pg_rel <- .projected_grad_norm(fit$s, grad) / (sqrt(sum(grad^2)) + 1e-12)

  expect_lt(pg_rel, 5e-3)
  expect_lte(max(diff(fit$objective_trace)), 1e-6)
})

test_that("graph_anchored_mfa orthonormal objective stays monotone on centered random data", {
  set.seed(4)
  N <- sample(18:30, 1)
  q <- sample(3:5, 1)
  p1 <- sample(4:7, 1)
  p2 <- sample(4:7, 1)
  Y <- matrix(rnorm(N * q), N, q)
  X1 <- matrix(rnorm(sample(12:20, 1) * p1), ncol = p1)
  X2 <- matrix(rnorm(sample(12:20, 1) * p2), ncol = p2)
  row_index <- list(
    X1 = sample.int(N, nrow(X1), replace = TRUE),
    X2 = sample.int(N, nrow(X2), replace = TRUE)
  )

  overlap <- sample.int(min(p1, p2), 1)
  nm1 <- paste0("x1_", seq_len(p1))
  nm2 <- paste0("x2_", seq_len(p2))
  nm1[seq_len(overlap)] <- paste0("shared", seq_len(overlap))
  nm2[seq_len(overlap)] <- paste0("shared", seq_len(overlap))
  colnames(X1) <- nm1
  colnames(X2) <- nm2

  fit <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = row_index,
    ncomp = min(2, q),
    preproc = multivarious::center(),
    normalization = "None",
    score_constraint = "orthonormal",
    graph_lambda = 0,
    score_graph_lambda = 0,
    ridge = 1e-6,
    max_iter = 25,
    tol = 1e-7
  )

  expect_lte(max(diff(fit$objective_trace)), 1e-8)
})

test_that("graph_anchored_mfa objective trace matches the fitted penalized objective", {
  sim <- .graph_anchor_score_graph_fixture()
  score_graph <- sim$score_graph_edges
  feature_graph <- data.frame(
    block1 = "X1",
    feature1 = colnames(sim$X$X1),
    block2 = "X2",
    feature2 = colnames(sim$X$X2),
    weight = 1,
    stringsAsFactors = FALSE
  )

  fit <- graph_anchored_mfa(
    Y = sim$Y,
    X = sim$X,
    row_index = sim$row_index,
    ncomp = sim$ncomp,
    preproc = multivarious::pass(),
    normalization = "None",
    feature_graph = feature_graph,
    graph_lambda = 3,
    score_graph = score_graph,
    score_graph_lambda = 5,
    max_iter = 80,
    tol = 1e-9,
    ridge = 1e-8
  )

  obj_fit <- muscal:::.gamfa_objective(
    Y = sim$Y,
    S = multivarious::scores(fit),
    B = fit$B,
    X_list = sim$X,
    V_list = fit$V_list,
    row_index = fit$row_index,
    alpha_y = fit$alpha_blocks[[1]],
    alpha_blocks = unname(fit$alpha_blocks[-1]),
    graph_laplacian = fit$graph_laplacian,
    graph_lambda = fit$graph_lambda,
    score_graph_laplacian = fit$score_graph_laplacian,
    score_graph_lambda = fit$score_graph_lambda,
    ridge = fit$ridge
  )

  expect_equal(unname(tail(fit$objective_trace, 1)), unname(obj_fit), tolerance = 1e-8)
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

test_that("graph_anchored_mfa matches nested row_index by names, not list position", {
  set.seed(1801)
  N <- 36
  Y <- matrix(rnorm(N * 3), N, 3)

  X_nested <- list(
    S1 = list(
      D1 = matrix(rnorm(12 * 4), 12, 4),
      D2 = matrix(rnorm(10 * 3), 10, 3)
    ),
    S2 = list(
      D1 = matrix(rnorm(11 * 4), 11, 4)
    )
  )
  row_nested <- list(
    S1 = list(
      D1 = sort(sample.int(N, 12)),
      D2 = sort(sample.int(N, 10))
    ),
    S2 = list(
      D1 = sort(sample.int(N, 11))
    )
  )

  fit_ref <- graph_anchored_mfa(
    Y = Y,
    X = X_nested,
    row_index = row_nested,
    ncomp = 2,
    preproc = multivarious::pass(),
    normalization = "None",
    graph_lambda = 0,
    max_iter = 60,
    tol = 1e-10
  )

  row_perm <- row_nested[c("S2", "S1")]
  row_perm$S1 <- row_perm$S1[c("D2", "D1")]

  fit_perm <- graph_anchored_mfa(
    Y = Y,
    X = X_nested,
    row_index = row_perm,
    ncomp = 2,
    preproc = multivarious::pass(),
    normalization = "None",
    graph_lambda = 0,
    max_iter = 60,
    tol = 1e-10
  )

  expect_lt(.score_subspace_rel(multivarious::scores(fit_ref), multivarious::scores(fit_perm)), 1e-10)
  expect_equal(names(fit_perm$row_index), names(fit_ref$row_index))
})

test_that("feature_graph='colnames' links across blocks but not within-block duplicates", {
  set.seed(1802)
  Y <- matrix(rnorm(24 * 3), 24, 3)
  X1 <- matrix(rnorm(24 * 3), 24, 3)
  X2 <- matrix(rnorm(24 * 2), 24, 2)
  colnames(X1) <- c("roi_a", "roi_a", "roi_b")
  colnames(X2) <- c("roi_a", "roi_c")

  fit <- graph_anchored_mfa(
    Y = Y,
    X = list(X1 = X1, X2 = X2),
    row_index = list(X1 = seq_len(nrow(Y)), X2 = seq_len(nrow(Y))),
    ncomp = 2,
    preproc = multivarious::pass(),
    normalization = "None",
    feature_graph = "colnames",
    graph_lambda = 1
  )

  A <- fit$graph_adjacency
  expect_equal(as.numeric(A[1, 2]), 0)
  expect_gt(as.numeric(A[1, 4]), 0)
  expect_gt(as.numeric(A[2, 4]), 0)
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

  bad_laplacian <- Matrix::Matrix(matrix(c(1, -2, -2, 1), 2, 2), sparse = TRUE)
  expect_error(
    graph_anchored_mfa(
      Y = Y,
      X = list(X1 = X1, X2 = X2),
      row_index = list(X1 = seq_len(N), X2 = seq_len(N)),
      ncomp = 2,
      feature_graph = as.matrix(Matrix::bdiag(bad_laplacian, Matrix::Diagonal(ncol(X1) + ncol(X2) - 2))),
      graph_form = "laplacian",
      graph_lambda = 1
    )
  )
})

test_that("graph_anchored_mfa validates malformed score-graph inputs", {
  set.seed(119)
  N <- 20
  Y <- matrix(rnorm(N * 3), N, 3)
  X1 <- matrix(rnorm(10 * 4), 10, 4)
  X2 <- matrix(rnorm(10 * 5), 10, 5)
  idx <- list(X1 = seq_len(10), X2 = 11:20)

  bad_graph <- Matrix::Diagonal(3)
  bad_edges <- data.frame(row1 = 1, row2 = 2, weight = -1)

  expect_error(
    graph_anchored_mfa(
      Y = Y,
      X = list(X1 = X1, X2 = X2),
      row_index = idx,
      ncomp = 2,
      score_graph = bad_graph,
      score_graph_form = "adjacency",
      score_graph_lambda = 1
    )
  )

  expect_error(
    graph_anchored_mfa(
      Y = Y,
      X = list(X1 = X1, X2 = X2),
      row_index = idx,
      ncomp = 2,
      score_graph = bad_edges,
      score_graph_lambda = 1
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
  expect_true(fit$fit_spec$refit_supported)
  expect_setequal(fit$oos_types, c("response", "scores", "reconstruction"))
  expect_equal(dim(scores_new), c(4, 2))
  expect_equal(dim(scores_from_predict), c(4, 2))
  expect_equal(dim(yhat), c(4, ncol(Y)))
  expect_equal(dim(xhat), dim(new_rows))
  expect_equal(unname(scores_from_predict), unname(scores_new), tolerance = 1e-10)
})
