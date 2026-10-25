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

# Component signs are arbitrary and can differ between fitting paths (the
# in-memory path's block weights come from a randomized SVD). Flip `fit`'s
# components to agree with `ref` before comparing.
.bs_align <- function(fit, ref) {
  sg <- sign(colSums(multivarious::scores(fit) * multivarious::scores(ref)))
  sg[sg == 0] <- 1
  flip <- function(M) sweep(M, 2L, sg, `*`)
  fit$s <- flip(fit$s)
  fit$v <- flip(fit$v)
  fit$V_list <- lapply(fit$V_list, flip)
  fit$partial_scores <- lapply(fit$partial_scores, flip)
  if (!is.null(fit$cor_loadings)) fit$cor_loadings <- flip(fit$cor_loadings)
  fit
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
  fit_src <- .bs_align(fit_src, fit_mem)

  expect_s3_class(fit_src, "aligned_mfa")
  expect_equal(fit_src$strategy, "gram")
  expect_equal(multivarious::scores(fit_src), multivarious::scores(fit_mem), tolerance = 1e-8)
  expect_equal(fit_src$objective_trace, fit_mem$objective_trace, tolerance = 1e-8)
  expect_equal(fit_src$alpha_blocks, fit_mem$alpha_blocks, tolerance = 1e-10)

  # Loadings are the fit's own loadings mapped back to feature space.
  for (k in seq_along(sim$X)) {
    expect_equal(unname(fit_src$V_list[[k]]), unname(fit_mem$V_list[[k]]), tolerance = 1e-8)
  }
  expect_equal(unname(fit_src$v), unname(fit_mem$v), tolerance = 1e-8)
  expect_equal(fit_src$block_fit$r2, fit_mem$block_fit$r2, tolerance = 1e-8)
  expect_equal(unname(fit_src$cor_loadings), unname(fit_mem$cor_loadings), tolerance = 1e-8)
  expect_equal(unname(fit_src$partial_scores), unname(fit_mem$partial_scores), tolerance = 1e-8)

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
  fit_src <- .bs_align(fit_src, fit_mem)
  expect_equal(multivarious::scores(fit_src), multivarious::scores(fit_mem), tolerance = 1e-8)
})

test_that("Gram-path fits project new rows like the in-memory fit", {
  sim <- .bs_sim(seed = 14)
  fit_mem <- .bs_fit(sim$X, sim$idx, sim$N, sim$K, score_constraint = "orthonormal")
  fit_src <- .bs_fit(as_block_source(sim$X, chunk_size = 20), sim$idx, sim$N, sim$K,
                     score_constraint = "orthonormal")
  fit_src <- .bs_align(fit_src, fit_mem)
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
  expect_equal(multivarious::scores(.bs_align(fit_list_gram, fit_mem)), multivarious::scores(fit_mem),
               tolerance = 1e-8)

  refit <- suppressMessages(fit_list_gram$fit_spec$refit$fit_fn(fit_list_gram$fit_spec$refit$data))
  expect_equal(multivarious::scores(refit), multivarious::scores(fit_list_gram), tolerance = 1e-10)
  expect_s3_class(fit_list_gram$fit_spec$refit$bootstrap_fn(fit_list_gram$fit_spec$refit$data)$X,
                  "muscal_block_source")
})

test_that("aligned_mfa Gram path rejects unsupported settings", {
  sim <- .bs_sim(seed = 16)
  src <- as_block_source(sim$X)
  expect_error(
    aligned_mfa(src, sim$idx, N = sim$N, ncomp = 2, feature_groups = "colnames",
                feature_lambda = 1, strategy = "gram"),
    "cannot use the feature prior"
  )
  expect_error(
    aligned_mfa(src, sim$idx, N = sim$N, ncomp = 2, preproc = multivarious::colscale(type = "unit")),
    "column-local"
  )
  bad_idx <- sim$idx
  bad_idx$S1 <- bad_idx$S1[-1]
  expect_error(aligned_mfa(src, bad_idx, N = sim$N, ncomp = 2), "has length")
})

# Blocks with shared feature names, so feature groups link features across blocks.
.bs_sim_named <- function(seed = 21, N = 36, K = 2, n_k = c(18, 22, 15), p = 60) {
  sim <- .bs_sim(seed = seed, N = N, K = K, n_k = n_k, p_k = rep(p, length(n_k)))
  sim$X <- lapply(sim$X, function(x) { colnames(x) <- paste0("f", seq_len(ncol(x))); x })
  sim
}

test_that("streaming path reproduces the in-memory fit with a feature prior", {
  sim <- .bs_sim_named()
  counted <- .bs_counting_source(sim$X, chunk_size = 23)
  counted$src$colnames <- lapply(sim$X, colnames)
  args <- list(feature_groups = "colnames", feature_lambda = 2, score_constraint = "orthonormal")

  fit_mem <- do.call(.bs_fit, c(list(sim$X, sim$idx, sim$N, sim$K), args))
  fit_str <- do.call(.bs_fit, c(list(counted$src, sim$idx, sim$N, sim$K), args))
  fit_str <- .bs_align(fit_str, fit_mem)

  expect_equal(fit_str$strategy, "stream")
  expect_equal(multivarious::scores(fit_str), multivarious::scores(fit_mem), tolerance = 1e-8)
  expect_equal(fit_str$objective_trace, unname(fit_mem$objective_trace), tolerance = 1e-8)
  expect_equal(fit_str$alpha_blocks, fit_mem$alpha_blocks, tolerance = 1e-10)
  expect_equal(unname(fit_str$V_list), unname(fit_mem$V_list), tolerance = 1e-8)
  expect_equal(fit_str$block_fit$r2, fit_mem$block_fit$r2, tolerance = 1e-8)
  expect_equal(unname(fit_str$cor_loadings), unname(fit_mem$cor_loadings), tolerance = 1e-8)
  expect_equal(unname(fit_str$partial_scores), unname(fit_mem$partial_scores), tolerance = 1e-8)

  # The prior matters: it changes the fit relative to no prior.
  fit_none <- .bs_fit(sim$X, sim$idx, sim$N, sim$K, score_constraint = "orthonormal")
  expect_gt(max(abs(unlist(fit_none$V_list) - unlist(fit_mem$V_list))), 1e-3)

  # One read of the data per iteration plus three fixed passes.
  n_chunks <- sum(ceiling(vapply(sim$X, ncol, integer(1L)) / 23))
  expect_equal(fit_str$passes, length(fit_str$objective_trace) + 3L)
  expect_equal(length(counted$log$reads), fit_str$passes * n_chunks)
  expect_true(all(vapply(counted$log$reads, function(r) length(r$cols), integer(1L)) <= 23))
})

test_that("streaming path matches in memory without a prior, with explicit groups, and at max_iter", {
  sim <- .bs_sim_named(seed = 22)
  src <- as_block_source(sim$X, chunk_size = 17)

  fit_mem <- .bs_fit(sim$X, sim$idx, sim$N, sim$K, score_constraint = "orthonormal")
  fit_str <- .bs_fit(src, sim$idx, sim$N, sim$K, score_constraint = "orthonormal", strategy = "stream")
  fit_str <- .bs_align(fit_str, fit_mem)
  expect_equal(multivarious::scores(fit_str), multivarious::scores(fit_mem), tolerance = 1e-8)

  groups <- data.frame(
    block = rep(c("S1", "S2", "S3"), each = 10),
    feature = rep(paste0("f", 1:10), 3),
    group = rep(paste0("g", 1:10), 3),
    weight = rep(c(1, 0.5, 2), each = 10)
  )
  capped <- list(feature_groups = groups, feature_lambda = 1, max_iter = 4, tol = 0,
                 score_constraint = "orthonormal")
  fit_mem <- suppressMessages(do.call(aligned_mfa, c(list(sim$X, sim$idx, N = sim$N, ncomp = sim$K), capped)))
  fit_str <- suppressMessages(do.call(aligned_mfa, c(list(src, sim$idx, N = sim$N, ncomp = sim$K), capped)))
  expect_length(fit_str$objective_trace, 4L)
  expect_equal(fit_str$objective_trace, unname(fit_mem$objective_trace), tolerance = 1e-8)
})

test_that("streaming path with score_constraint = 'none' lands where in-memory fits do", {
  # With a feature prior the "none" path converges slowly and its end point
  # depends on the random start (in memory too), so compare against the spread
  # of in-memory fits rather than one of them.
  sim <- .bs_sim_named(seed = 23)
  src <- as_block_source(sim$X, chunk_size = 25)
  args <- list(N = sim$N, ncomp = sim$K, feature_groups = "colnames", feature_lambda = 1,
               max_iter = 60)
  mem_obj <- vapply(1:4, function(seed) {
    set.seed(seed)
    tail(suppressMessages(do.call(aligned_mfa, c(list(sim$X, sim$idx), args)))$objective_trace, 1)
  }, numeric(1L))
  set.seed(1)
  fit_str <- suppressMessages(do.call(aligned_mfa, c(list(src, sim$idx), args)))

  expect_equal(fit_str$strategy, "stream")
  expect_true(all(is.finite(fit_str$objective_trace)))
  expect_lt(tail(fit_str$objective_trace, 1), max(mem_obj) * 1.01)
  expect_gt(tail(fit_str$objective_trace, 1), min(mem_obj) * 0.99)
})

test_that("streaming fits project, predict and refit", {
  sim <- .bs_sim_named(seed = 24)
  src <- as_block_source(sim$X, chunk_size = 30)
  fit <- .bs_fit(src, sim$idx, sim$N, sim$K, feature_groups = "colnames", feature_lambda = 1,
                 score_constraint = "orthonormal")
  fit_mem <- .bs_fit(sim$X, sim$idx, sim$N, sim$K, feature_groups = "colnames", feature_lambda = 1,
                     score_constraint = "orthonormal")
  new_rows <- sim$X$S1[1:4, , drop = FALSE]
  expect_equal(project(.bs_align(fit, fit_mem), new_rows, block = "S1"), project(fit_mem, new_rows, block = "S1"),
               tolerance = 1e-8)

  refit <- suppressMessages(fit$fit_spec$refit$fit_fn(fit$fit_spec$refit$data))
  expect_equal(refit$strategy, "stream")
  expect_equal(multivarious::scores(refit), multivarious::scores(fit), tolerance = 1e-10)
})

test_that("block_source column names are validated and needed for named groups", {
  X <- list(A = matrix(rnorm(20), 5), B = matrix(rnorm(30), 5))
  expect_error(
    block_source(function(k, cols) X[[k]][, cols, drop = FALSE], nrow = c(5, 5), ncol = c(4, 6),
                 colnames = list(letters[1:3], NULL)),
    "must have length"
  )
  src <- block_source(function(k, cols) X[[k]][, cols, drop = FALSE], nrow = c(5, 5), ncol = c(4, 6))
  expect_error(
    suppressMessages(aligned_mfa(src, list(1:5, 1:5), ncomp = 1, feature_groups = "colnames",
                                 feature_lambda = 1)),
    "colnames"
  )
})

test_that("block-source bootstrap draws the same replicate data as in memory", {
  sim <- .bs_sim(seed = 31)
  data_mem <- list(X = sim$X, row_index = sim$idx, N = sim$N)
  data_src <- list(X = as_block_source(sim$X, chunk_size = 30), row_index = sim$idx, N = sim$N)
  for (seed in 1:3) {
    set.seed(seed)
    boot_mem <- muscal:::.muscal_bootstrap_aligned_data(data_mem)
    set.seed(seed)
    boot_src <- muscal:::.muscal_source_resampler(muscal:::.muscal_bootstrap_aligned_data)(data_src)
    expect_equal(unname(boot_src$row_index), unname(boot_mem$row_index))
    for (k in seq_along(sim$X)) {
      expect_equal(read_block(boot_src$X, k), unname(boot_mem$X[[k]]), ignore_attr = TRUE)
    }
  }
})

test_that("infer_muscal resamples block-source fits like in-memory fits", {
  sim <- .bs_sim(seed = 32)
  src <- as_block_source(sim$X, chunk_size = 40)
  fit_mem <- .bs_fit(sim$X, sim$idx, sim$N, sim$K, score_constraint = "orthonormal")
  fit_gram <- .bs_fit(src, sim$idx, sim$N, sim$K, score_constraint = "orthonormal")
  boot_mem <- suppressMessages(infer_muscal(fit_mem, method = "bootstrap", nrep = 3, seed = 7))
  boot_gram <- suppressMessages(infer_muscal(fit_gram, method = "bootstrap", nrep = 3, seed = 7))
  expect_equal(boot_gram$n_success, 3L)
  # The first replicate's resample is drawn before any fit consumes random
  # numbers, so it is identical across paths; later draws can differ because
  # the in-memory path's randomized SVD consumes the stream.
  expect_equal(boot_gram$replicates[[1]], boot_mem$replicates[[1]], tolerance = 1e-6)
  expect_true(all(vapply(boot_gram$replicates, function(r) all(is.finite(unlist(r))), logical(1L))))

  perm_mem <- suppressMessages(infer_muscal(fit_mem, method = "permutation", nrep = 3, seed = 8))
  perm_gram <- suppressMessages(infer_muscal(fit_gram, method = "permutation", nrep = 3, seed = 8))
  expect_equal(perm_gram$replicates[[1]], perm_mem$replicates[[1]], tolerance = 1e-6)
  expect_equal(perm_gram$n_success, 3L)
})

test_that("gram_rank truncates compression and reports discarded variance", {
  sim <- .bs_sim(seed = 33)
  src <- as_block_source(sim$X, chunk_size = 35)
  fit_exact <- .bs_fit(src, sim$idx, sim$N, sim$K, score_constraint = "orthonormal")
  fit_full_rank <- .bs_fit(src, sim$idx, sim$N, sim$K, score_constraint = "orthonormal", gram_rank = 100)
  fit_trunc <- .bs_fit(src, sim$idx, sim$N, sim$K, score_constraint = "orthonormal", gram_rank = 6)

  expect_equal(multivarious::scores(fit_full_rank), multivarious::scores(fit_exact), tolerance = 1e-10)
  expect_true(all(fit_exact$gram_discarded < 1e-10))
  expect_true(all(fit_trunc$gram_rank <= 6))
  expect_true(all(fit_trunc$gram_discarded > 0 & fit_trunc$gram_discarded < 1))

  # Truncation keeps the leading structure: the same latent subspace.
  Q1 <- qr.Q(qr(multivarious::scores(fit_exact)))
  Q2 <- qr.Q(qr(multivarious::scores(fit_trunc)))
  expect_gt(min(svd(crossprod(Q1, Q2))$d), 0.95)

  # Truncated loadings lie in the kept right-singular subspace, and the reported
  # R^2 is the actual R^2 of reconstructing the data with them.
  S <- multivarious::scores(fit_trunc)
  for (k in seq_along(sim$X)) {
    Sk <- S[sim$idx[[k]], , drop = FALSE]
    Xc <- scale(sim$X[[k]], scale = FALSE)
    Vk <- unname(fit_trunc$V_list[[k]])
    W <- svd(Xc)$v[, seq_len(fit_trunc$gram_rank[[k]]), drop = FALSE]
    expect_equal(W %*% crossprod(W, Vk), Vk, tolerance = 1e-8)
    r2_ref <- 1 - sum((Xc - Sk %*% t(Vk))^2) / sum(Xc^2)
    expect_equal(fit_trunc$block_fit$r2[[k]], r2_ref, tolerance = 1e-8)
  }
  expect_error(.bs_fit(src, sim$idx, sim$N, sim$K, gram_rank = 0), "gram_rank")
})

test_that("feature_lambda_scale = 'block' is scale-free and streams exactly", {
  sim <- .bs_sim_named(seed = 34)
  rescaled <- sim$X
  rescaled$S2 <- rescaled$S2 * 1000
  args <- list(feature_groups = "colnames", feature_lambda = 1, score_constraint = "orthonormal",
               ridge = 0)
  fit_block <- do.call(.bs_fit, c(list(sim$X, sim$idx, sim$N, sim$K, feature_lambda_scale = "block"), args))
  fit_block2 <- do.call(.bs_fit, c(list(rescaled, sim$idx, sim$N, sim$K, feature_lambda_scale = "block"), args))
  expect_equal(multivarious::scores(fit_block2), multivarious::scores(fit_block), tolerance = 1e-8)
  expect_equal(fit_block2$V_list$S2 / 1000, fit_block$V_list$S2, tolerance = 1e-8)

  fit_raw <- do.call(.bs_fit, c(list(sim$X, sim$idx, sim$N, sim$K), args))
  fit_raw2 <- do.call(.bs_fit, c(list(rescaled, sim$idx, sim$N, sim$K), args))
  expect_gt(max(abs(multivarious::scores(fit_raw2) - multivarious::scores(fit_raw))), 1e-3)

  fit_stream <- do.call(.bs_fit, c(list(as_block_source(sim$X, chunk_size = 19), sim$idx, sim$N, sim$K,
                                        feature_lambda_scale = "block"), args))
  expect_equal(fit_stream$strategy, "stream")
  expect_equal(multivarious::scores(fit_stream), multivarious::scores(fit_block), tolerance = 1e-8)
  expect_equal(fit_stream$objective_trace, unname(fit_block$objective_trace), tolerance = 1e-8)

  refit <- suppressMessages(fit_stream$fit_spec$refit$fit_fn(fit_stream$fit_spec$refit$data))
  expect_equal(multivarious::scores(refit), multivarious::scores(fit_stream), tolerance = 1e-10)
})
