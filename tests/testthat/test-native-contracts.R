test_that("compiled row systems and graph operators match base-R oracles", {
  expect_true(muscal:::.muscal_has_native_symbol("muscal_apply_row_system_cpp"))
  S <- matrix(c(1, -2, 3, 4, 0.5, -1), 3, 2)
  A <- list(matrix(c(2, 1, 1, 3), 2), diag(c(4, 2)), matrix(c(3, -1, -1, 2), 2))
  expected <- t(vapply(seq_len(3), function(i) as.vector(A[[i]] %*% S[i, ]), numeric(2)))
  expect_equal(muscal:::muscal_apply_row_system_cpp(A, S), expected, tolerance = 1e-12)
  L <- matrix(c(1, -1, 0, -1, 2, -1, 0, -1, 1), 3)
  for (graph in list(L, as(Matrix::Matrix(L, sparse = TRUE), "generalMatrix"))) {
    expect_equal(muscal:::muscal_apply_row_operator_cpp(A, S, graph, 0.7, 0.2),
      expected + 0.7 * L %*% S + 0.2 * S, tolerance = 1e-12)
  }
  expect_equal(muscal:::muscal_apply_row_operator_cpp(A, S, NULL, 0, 0), expected)
})

test_that("compiled row kernels reject malformed data before Eigen indexing", {
  S <- matrix(seq_len(6) * 1.0, 3, 2)
  A <- rep(list(diag(2)), 3)
  expect_error(muscal:::muscal_apply_row_system_cpp(A[-1], S), "length")
  bad <- A; bad[[2]] <- matrix(1, 1, 2)
  expect_error(muscal:::muscal_apply_row_system_cpp(bad, S), "dimensions")
  bad[[2]] <- "bad"
  expect_error(muscal:::muscal_apply_row_system_cpp(bad, S), "numeric")
  bad <- A; bad[[1]][1] <- Inf
  expect_error(muscal:::muscal_apply_row_system_cpp(bad, S), "finite")
  bad_S <- S; bad_S[1] <- NA_real_
  expect_error(muscal:::muscal_apply_row_system_cpp(A, bad_S), "finite")
  expect_error(muscal:::muscal_apply_row_system_cpp(A, as.vector(S)))
  for (graph in list(diag(2), Matrix::Diagonal(2))) {
    expect_error(muscal:::muscal_apply_row_operator_cpp(A, S, graph, 0, 0), "dimensions|dgCMatrix|sparse")
  }
  L <- diag(3); L[1] <- Inf
  expect_error(muscal:::muscal_apply_row_operator_cpp(A, S, L, 1, 0), "finite")
  for (lambda in c(-1, NA_real_, Inf)) {
    expect_error(muscal:::muscal_apply_row_operator_cpp(A, S, NULL, lambda, 0), "nonnegative")
    expect_error(muscal:::muscal_apply_row_operator_cpp(A, S, NULL, 0, lambda), "nonnegative")
  }
})

test_that("compiled grouped rows agree with explicit sums and validate integer groups", {
  X <- matrix(c(1, 3, 2, -2, 5, 4, 7, 8), 4, 2)
  idx <- c(3, 1, 3, 1)
  out <- muscal:::muscal_rowsum_counts_cpp(X, idx, 4)
  expect_equal(out$sums, rbind(colSums(X[c(2,4), ]), c(0,0), colSums(X[c(1,3), ]), c(0,0)))
  expect_equal(out$counts, c(2L, 0L, 2L, 0L))
  for (N in list(-1, NA_real_, Inf, 1.5, c(3, 4), "4")) {
    expect_error(muscal:::muscal_rowsum_counts_cpp(X, idx, N), "N must")
  }
  for (bad in list(c(1, 2), c(0, 1, 1, 1), c(1, 1, 1, 5), c(1, NA, 1, 1), c(1, 1.5, 1, 1))) {
    expect_error(muscal:::muscal_rowsum_counts_cpp(X, bad, 4), "idx")
    expect_error(muscal:::.muscal_rowsum_counts(X, bad, 4), "idx")
  }
  bad <- X; bad[1] <- Inf
  expect_error(muscal:::muscal_rowsum_counts_cpp(bad, idx, 4), "finite")
  empty <- muscal:::muscal_rowsum_counts_cpp(matrix(numeric(0), 0, 2), integer(0), 0)
  expect_equal(dim(empty$sums), c(0L, 2L))
  expect_length(empty$counts, 0)
})
