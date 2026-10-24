library(testthat)
library(muscal)

.bs_sim <- function(seed = 11, N = 40, K = 3, n_k = c(20, 16, 25), p_k = c(90, 70, 120)) {
  set.seed(seed)
  Z <- matrix(rnorm(N * K), N, K)
  idx <- lapply(n_k, function(n) sort(sample.int(N, n)))
  X <- Map(function(n, p, ix) {
    Z[ix, , drop = FALSE] %*% matrix(rnorm(K * p), K, p) +
      matrix(rnorm(n * p, sd = 1.5), n, p) + rep(rnorm(p, sd = 5), each = n)
  }, n_k, p_k, idx)
  names(X) <- names(idx) <- paste0("S", seq_along(X))
  list(X = X, idx = idx, N = N, K = K)
}

# A loader-backed source that records every read.
.bs_counting_source <- function(X, chunk_size) {
  log <- new.env()
  log$reads <- list()
  src <- block_source(
    loader = function(k, cols) {
      log$reads[[length(log$reads) + 1L]] <- list(k = k, cols = cols)
      X[[k]][, cols, drop = FALSE]
    },
    nrow = vapply(X, nrow, integer(1L)),
    ncol = vapply(X, ncol, integer(1L)),
    names = names(X),
    chunk_size = chunk_size
  )
  list(src = src, log = log)
}

.bs_fit <- function(X, idx, N, K, ...) {
  suppressMessages(aligned_mfa(X, idx, N = N, ncomp = K, max_iter = 300, tol = 1e-12, ...))
}

test_that("block_source validates loaders and reads chunks", {
  X <- list(A = matrix(rnorm(40), 10), B = matrix(rnorm(70), 10))
  src <- as_block_source(X, chunk_size = 3)

  expect_s3_class(src, "muscal_block_source")
  expect_equal(names(src), c("A", "B"))
  expect_equal(length(src), 2L)
  expect_equal(read_block(src, "B"), X$B)
  expect_equal(read_block(src, 1, cols = 2:3), X$A[, 2:3])
  chunks <- muscal:::.muscal_source_chunks(src, 2)
  expect_equal(unlist(chunks), seq_len(7))
  expect_true(all(lengths(chunks) <= 3))

  bad_dims <- block_source(function(k, cols) matrix(0, 3, length(cols)), nrow = c(4, 4), ncol = c(2, 2))
  expect_error(read_block(bad_dims, 1), "expected 4 x 2")
  with_na <- block_source(function(k, cols) matrix(NA_real_, 4, length(cols)), nrow = c(4, 4), ncol = c(2, 2))
  expect_error(read_block(with_na, 1), "missing or infinite")
  expect_error(block_source(function(k, cols) NULL, nrow = c(4, 4), ncol = 2), "equal length")
  expect_error(as_block_source(1:3), "Cannot convert")
})

test_that("streamed Gram matrices and preprocessing statistics match in-memory values", {
  sim <- .bs_sim()
  src <- as_block_source(sim$X, chunk_size = 17)
  for (kind in c("none", "center", "standardize")) {
    gram <- muscal:::.muscal_stream_gram(src, rep(kind, 3))
    for (k in seq_along(sim$X)) {
      Xp <- switch(kind,
        none = sim$X[[k]],
        center = scale(sim$X[[k]], scale = FALSE),
        standardize = scale(sim$X[[k]])
      )
      expect_equal(gram$G[[k]], tcrossprod(unname(Xp)), tolerance = 1e-10)
      expect_equal(gram$prep[[k]]$col_ss, unname(colSums(scale(Xp, scale = FALSE)^2)), tolerance = 1e-10)
    }
  }
})

test_that("column-local preprocessors are recognised and others rejected", {
  kinds <- muscal:::.muscal_stream_preproc_kinds
  expect_equal(kinds(NULL, 2), c("none", "none"))
  expect_equal(kinds(multivarious::pass(), 1), "none")
  expect_equal(kinds(multivarious::center(), 1), "center")
  expect_equal(kinds(multivarious::standardize(), 1), "standardize")
  expect_equal(kinds(list(multivarious::center(), multivarious::standardize()), 2),
               c("center", "standardize"))
  expect_error(kinds(multivarious::colscale(type = "unit"), 1), "column-local")
  expect_error(kinds(list(multivarious::center()), 2), "must match")
})

test_that("aligned_mfa Gram path reproduces the in-memory fit", {
  sim <- .bs_sim()
  counted <- .bs_counting_source(sim$X, chunk_size = 25)

  fit_mem <- .bs_fit(sim$X, sim$idx, sim$N, sim$K, score_constraint = "orthonormal")
  fit_src <- .bs_fit(counted$src, sim$idx, sim$N, sim$K, score_constraint = "orthonormal")

  expect_s3_class(fit_src, "aligned_mfa")
  expect_equal(fit_src$strategy, "gram")
  expect_equal(multivarious::scores(fit_src), multivarious::scores(fit_mem), tolerance = 1e-8)
  expect_equal(fit_src$objective_trace, fit_mem$objective_trace, tolerance = 1e-8)
  expect_equal(fit_src$alpha_blocks, fit_mem$alpha_blocks, tolerance = 1e-10)

  # Loadings are the ridge regressions on the final scores.
  S <- multivarious::scores(fit_mem)
  for (k in seq_along(sim$X)) {
    Sk <- S[sim$idx[[k]], , drop = FALSE]
    a <- fit_mem$alpha_blocks[[k]]
    Xc <- scale(sim$X[[k]], scale = FALSE)
    V_ref <- t(solve(a * crossprod(Sk) + fit_mem$ridge * diag(sim$K), a * crossprod(Sk, Xc)))
    expect_equal(unname(fit_src$V_list[[k]]), unname(V_ref), tolerance = 1e-8)
    expect_equal(unname(fit_src$V_list[[k]]), unname(fit_mem$V_list[[k]]), tolerance = 1e-4)
  }
  expect_equal(fit_src$block_fit$r2, fit_mem$block_fit$r2, tolerance = 1e-6)
  expect_equal(unname(fit_src$cor_loadings), unname(fit_mem$cor_loadings), tolerance = 1e-6)
  expect_equal(unname(fit_src$partial_scores), unname(fit_mem$partial_scores), tolerance = 1e-4)

  # Two passes over the data, never more than chunk_size columns per read.
  n_chunks <- sum(ceiling(vapply(sim$X, ncol, integer(1L)) / 25))
  expect_equal(length(counted$log$reads), 2L * n_chunks)
  expect_true(all(vapply(counted$log$reads, function(r) length(r$cols), integer(1L)) <= 25))
})

test_that("aligned_mfa Gram path with score_constraint = 'none' finds the same solution", {
  sim <- .bs_sim(seed = 12)
  set.seed(1)
  fit_mem <- .bs_fit(sim$X, sim$idx, sim$N, sim$K)
  set.seed(1)
  fit_src <- .bs_fit(as_block_source(sim$X, chunk_size = 30), sim$idx, sim$N, sim$K)

  Q1 <- qr.Q(qr(multivarious::scores(fit_mem)))
  Q2 <- qr.Q(qr(multivarious::scores(fit_src)))
  expect_gt(min(svd(crossprod(Q1, Q2))$d), 1 - 1e-6)
  expect_equal(tail(fit_src$objective_trace, 1), tail(fit_mem$objective_trace, 1), tolerance = 1e-6)
})

test_that("aligned_mfa Gram path supports standardize and per-block preprocessing", {
  sim <- .bs_sim(seed = 13)
  pre <- list(multivarious::standardize(), multivarious::center(), multivarious::pass())
  fit_mem <- .bs_fit(sim$X, sim$idx, sim$N, sim$K, preproc = pre, score_constraint = "orthonormal")
  fit_src <- .bs_fit(as_block_source(sim$X, chunk_size = 40), sim$idx, sim$N, sim$K,
                     preproc = pre, score_constraint = "orthonormal")
  expect_equal(multivarious::scores(fit_src), multivarious::scores(fit_mem), tolerance = 1e-8)
})

test_that("Gram-path fits project new rows like the in-memory fit", {
  sim <- .bs_sim(seed = 14)
  fit_mem <- .bs_fit(sim$X, sim$idx, sim$N, sim$K, score_constraint = "orthonormal")
  fit_src <- .bs_fit(as_block_source(sim$X, chunk_size = 20), sim$idx, sim$N, sim$K,
                     score_constraint = "orthonormal")
  new_rows <- sim$X$S2[1:5, , drop = FALSE] + 0.1
  expect_equal(
    project(fit_src, new_rows, block = "S2"),
    project(fit_mem, new_rows, block = "S2"),
    tolerance = 1e-4
  )
  expect_equal(
    predict(fit_src, new_rows, block = "S2", type = "reconstruction"),
    predict(fit_mem, new_rows, block = "S2", type = "reconstruction"),
    tolerance = 1e-4
  )
})

test_that("aligned_mfa strategies interoperate and refit reuses the source", {
  sim <- .bs_sim(seed = 15)
  src <- as_block_source(sim$X, chunk_size = 50)
  fit_list_gram <- .bs_fit(sim$X, sim$idx, sim$N, sim$K, strategy = "gram", score_constraint = "orthonormal")
  fit_src_mem <- .bs_fit(src, sim$idx, sim$N, sim$K, strategy = "in_memory", score_constraint = "orthonormal")
  fit_mem <- .bs_fit(sim$X, sim$idx, sim$N, sim$K, score_constraint = "orthonormal")

  expect_equal(fit_list_gram$strategy, "gram")
  expect_equal(multivarious::scores(fit_src_mem), multivarious::scores(fit_mem))
  expect_equal(multivarious::scores(fit_list_gram), multivarious::scores(fit_mem), tolerance = 1e-8)

  refit <- suppressMessages(fit_list_gram$fit_spec$refit$fit_fn(fit_list_gram$fit_spec$refit$data))
  expect_equal(multivarious::scores(refit), multivarious::scores(fit_list_gram), tolerance = 1e-10)
  expect_error(fit_list_gram$fit_spec$refit$bootstrap_fn(fit_list_gram$fit_spec$refit$data),
               "not yet available")
})

test_that("aligned_mfa Gram path rejects unsupported settings", {
  sim <- .bs_sim(seed = 16)
  src <- as_block_source(sim$X)
  expect_error(
    aligned_mfa(src, sim$idx, N = sim$N, ncomp = 2, feature_groups = "colnames", feature_lambda = 1),
    "feature prior"
  )
  expect_error(
    aligned_mfa(src, sim$idx, N = sim$N, ncomp = 2, preproc = multivarious::colscale(type = "unit")),
    "column-local"
  )
  bad_idx <- sim$idx
  bad_idx$S1 <- bad_idx$S1[-1]
  expect_error(aligned_mfa(src, bad_idx, N = sim$N, ncomp = 2), "has length")
})
