.graph_line_coordinates <- function(p) cbind(x = seq_len(p), y = rep(0, p), z = rep(0, p))

test_that("clusterwise neighbors have explicit zero and singleton semantics", {
  for (p in c(1L, 2L, 4L)) {
    coords <- list(.graph_line_coordinates(p))
    zero <- muscal:::spatial_constraints(coords, k_nn = 0)
    expect_equal(as.matrix(zero), matrix(0, p, p))
    for (k in unique(c(1L, p - 1L, p, p + 5L))) {
      A <- muscal:::spatial_constraints(coords, k_nn = k)
      expect_equal(dim(A), c(p, p))
      expect_true(all(is.finite(A)))
      expect_equal(as.numeric(Matrix::diag(A)), rep(0, p))
      if (k >= p - 1L) expect_equal(as.matrix(A), matrix(1, p, p) - diag(p))
    }
  }
})

test_that("clusterwise graphs reject invalid neighbor counts", {
  coords <- list(.graph_line_coordinates(4))
  for (k in list(-1, 1.5, NA_real_, Inf, c(1, 2))) {
    expect_error(muscal:::spatial_constraints(coords, k_nn = k), "k_nn")
  }
})

test_that("clusterwise public defaults fit fewer than seven pooled features", {
  X <- matrix(c(1, 3, 2, 5, 2, 1, 4, 3), 4, 2)
  coords <- list(.graph_line_coordinates(2), .graph_line_coordinates(2) + 0.2)
  fit <- penalized_mfa_clusterwise(list(X, X + 3), coords, ncomp = 1,
                                   lambda = 0, max_iter = 1)
  expect_s3_class(fit, "penalized_mfa_clusterwise")
  expect_equal(dim(fit$Sadj), c(4L, 4L))
  expect_true(all(is.finite(fit$Sadj)))
  one_feature <- penalized_mfa_clusterwise(
    list(X[, 1, drop = FALSE], X[, 2, drop = FALSE]),
    list(.graph_line_coordinates(1), .graph_line_coordinates(1) + 0.2),
    ncomp = 1, lambda = 0, max_iter = 1
  )
  expect_equal(dim(one_feature$Sadj), c(2L, 2L))
  expect_error(penalized_mfa_clusterwise(list(X[FALSE, ], X), coords, ncomp = 1),
               "nonempty")
})

test_that("clusterwise union graph retains lower-only neighbors and the correct penalty", {
  coords <- cbind(x = c(0, 1, 3, 10), y = 0, z = 0)
  expected <- matrix(c(0, 1, 0, 0, 1, 0, 1, 0,
                       0, 1, 0, 1, 0, 0, 1, 0), 4, 4, byrow = TRUE)
  A <- muscal:::spatial_constraints(list(coords), k_nn = 1)
  expect_equal(as.matrix(A), expected)
  L <- muscal:::cluster_graph(list(coords), list(k_nn = 1))
  V <- matrix(c(1, 0, 0, 1, 2, 1, 1, -1), 4, 2, byrow = TRUE)
  expect_equal(sum(V * (L %*% V)), 11)
  permutation <- c(4, 2, 1, 3)
  permuted <- muscal:::spatial_constraints(list(coords[permutation, ]), k_nn = 1)
  expect_equal(as.matrix(permuted)[order(permutation), order(permutation)], expected)
  blocks <- list(a = coords[1:2, ], b = coords[3:4, ])
  reordered <- muscal:::spatial_constraints(blocks[2:1], k_nn = 1)
  expect_equal(as.matrix(reordered)[c(3, 4, 1, 2), c(3, 4, 1, 2)], expected)
})

test_that("clusterwise boundary ties include every tied neighbor without self edges", {
  coords <- cbind(x = c(0, 1, -1, 10), y = 0, z = 0)
  expected <- matrix(0, 4, 4)
  expected[1, c(2, 3)] <- expected[c(2, 3), 1] <- 1
  expected[2, 4] <- expected[4, 2] <- 1
  A <- muscal:::spatial_constraints(list(coords), k_nn = 1)
  expect_equal(as.matrix(A), expected)
  for (p in list(4:1, c(2, 4, 1, 3))) {
    Ap <- muscal:::spatial_constraints(list(coords[p, ]), k_nn = 1)
    expect_equal(as.matrix(Ap)[order(p), order(p)], expected)
  }
  duplicates <- cbind(x = c(0, 0, 1), y = 0, z = 0)
  A <- muscal:::spatial_constraints(list(duplicates), k_nn = 1)
  expect_equal(as.matrix(A), matrix(1, 3, 3) - diag(3))
})

test_that("clusterwise neighbor membership is invariant to coordinate units", {
  coords <- cbind(x = c(0, 1, 3, 10), y = 0, z = 0)
  reference <- muscal:::spatial_constraints(list(coords), k_nn = 1)
  for (factor in c(1e-200, 1e-20, 1e20, 1e200)) {
    A <- muscal:::spatial_constraints(list(coords * factor), k_nn = 1)
    expect_equal(as.matrix(A), as.matrix(reference))
  }
})

test_that("normalized Laplacians agree with the analytic path and isolated-vertex oracle", {
  A <- Matrix::Matrix(matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3), sparse = TRUE)
  expected <- matrix(c(1, -1 / sqrt(2), 0,
                       -1 / sqrt(2), 1, -1 / sqrt(2),
                       0, -1 / sqrt(2), 1), 3, 3, byrow = TRUE)
  L <- muscal:::.pmfa_graph_laplacian(A, normalized = TRUE)
  expect_equal(as.matrix(L), expected, tolerance = 1e-12)
  expect_equal(eigen(as.matrix(L), symmetric = TRUE)$values, c(2, 1, 0), tolerance = 1e-12)
  raw <- muscal:::.pmfa_graph_laplacian(A, normalized = FALSE)
  expect_equal(as.numeric(Matrix::diag(raw)), c(1, 2, 1))
  expect_equal(as.matrix(muscal:::.pmfa_graph_laplacian(3 * A, TRUE)), expected)
  expect_equal(as.matrix(muscal:::.pmfa_graph_laplacian(3 * A, FALSE)), 3 * as.matrix(raw))
  isolated <- Matrix::bdiag(A, Matrix::Matrix(0, 1, 1))
  Li <- muscal:::.pmfa_graph_laplacian(isolated, TRUE)
  expect_equal(as.numeric(Li[4, ]), rep(0, 4))
  expect_true(all(is.finite(Li)))
  expect_equal(as.numeric(Li %*% c(1, sqrt(2), 1, 0)), rep(0, 4), tolerance = 1e-12)
})

test_that("public clusterwise normalization uses retained adjacency degrees", {
  X <- matrix(c(1, 3, 2, 5, 2, 1, 4, 3), 4, 2)
  coords <- list(cbind(x = c(0, 1), y = 0, z = 0), matrix(c(3, 0, 0), 1, 3))
  fit <- penalized_mfa_clusterwise(list(X, matrix(c(3, 1, 2, 5), 4, 1)),
    coords, ncomp = 1, lambda = 0, max_iter = 1,
    adjacency_opts = list(k_nn = 1), normalized_laplacian = TRUE)
  expect_equal(as.numeric(Matrix::diag(fit$Sadj)), rep(1, 3))
  expect_equal(fit$graph_degree, c(1, 2, 1))
  expect_true(fit$normalized_laplacian)
  graph <- muscal:::cluster_graph(coords, list(k_nn = 1), return_graph = TRUE)
  expect_equal(graph$degree, c(1, 2, 1))
  expect_equal(as.numeric(Matrix::diag(graph$laplacian)), c(1, 2, 1))
})
