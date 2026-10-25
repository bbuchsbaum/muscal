library(testthat)
library(muscal)

# Out-of-core (block source) fits through the shared Gram engine must match the
# in-memory fits. Component signs are arbitrary, so comparisons flip each
# component of `b` to agree with `a` first.
.gm_expect_same <- function(a, b, fields = character(), tolerance = 1e-7) {
  S_a <- multivarious::scores(a)
  S_b <- multivarious::scores(b)
  sg <- sign(colSums(S_a * S_b))
  sg[sg == 0] <- 1
  flip <- function(M) sweep(as.matrix(M), 2L, sg, `*`)
  expect_equal(unname(S_b %*% diag(sg, length(sg))), unname(S_a), tolerance = tolerance)
  for (f in fields) {
    x <- a[[f]]
    y <- b[[f]]
    if (is.data.frame(x)) {
      expect_equal(y$r2, x$r2, tolerance = tolerance, info = f)
    } else if (is.list(x)) {
      for (nm in names(x)) {
        expect_equal(unname(flip(y[[nm]])), unname(as.matrix(x[[nm]])), tolerance = tolerance,
                     info = paste(f, nm))
      }
    } else {
      expect_equal(unname(flip(y)), unname(as.matrix(x)), tolerance = tolerance, info = f)
    }
  }
}

.gm_shared_rows <- function(seed = 1, n = 30, K = 3, p = c(40, 25, 55)) {
  set.seed(seed)
  Z <- matrix(rnorm(n * K), n)
  X <- lapply(p, function(pk) Z %*% matrix(rnorm(K * pk), K) + matrix(rnorm(n * pk), n))
  names(X) <- paste0("B", seq_along(X))
  X
}

.gm_linked_rows <- function(seed = 2, N = 40, K = 2, n_k = c(25, 30, 22), p = c(45, 35, 60)) {
  set.seed(seed)
  Z <- matrix(rnorm(N * K), N)
  idx <- lapply(n_k, function(m) sort(sample.int(N, m)))
  X <- Map(function(ix, pk) {
    Z[ix, , drop = FALSE] %*% matrix(rnorm(K * pk), K) + matrix(rnorm(length(ix) * pk), length(ix))
  }, idx, p)
  names(X) <- names(idx) <- paste0("B", seq_along(X))
  Y <- Z %*% matrix(rnorm(K * 6), K) + matrix(rnorm(N * 6), N)
  list(X = X, idx = idx, Y = Y, N = N, K = K)
}

test_that("mfa on a block source matches the in-memory fit", {
  X <- .gm_shared_rows()
  src <- as_block_source(X, chunk_size = 17)
  for (normalization in c("MFA", "RV", "None")) {
    fit_mem <- suppressMessages(mfa(X, ncomp = 3, normalization = normalization))
    fit_src <- suppressMessages(mfa(src, ncomp = 3, normalization = normalization))
    expect_equal(fit_src$strategy, "gram")
    .gm_expect_same(fit_mem, fit_src, c("v", "ov", "cor_loadings", "partial_scores"))
  }
  Xc <- do.call(cbind, X)
  expect_equal(
    predict(fit_src, Xc, type = "reconstruction"),
    predict(fit_mem, Xc, type = "reconstruction"),
    tolerance = 1e-8
  )
  expect_equal(Matrix::diag(fit_src$A_genpca), Matrix::diag(fit_mem$A_genpca), tolerance = 1e-10)

  refit <- suppressMessages(fit_src$fit_spec$refit$fit_fn(fit_src$fit_spec$refit$data))
  .gm_expect_same(fit_src, refit)
  boot <- suppressMessages(infer_muscal(fit_src, method = "bootstrap", nrep = 2, seed = 1))
  expect_equal(boot$n_success, 2L)

  expect_error(mfa(src, ncomp = 2, normalization = "custom", A = rep(1, 120)), "column metric")
})

test_that("gram_rank truncates a block-source mfa and reports what was dropped", {
  X <- .gm_shared_rows(seed = 3)
  fit <- suppressMessages(mfa(as_block_source(X), ncomp = 2, gram_rank = 5))
  expect_true(all(fit$gram_rank == 5L))
  expect_true(all(fit$gram_discarded > 0 & fit$gram_discarded < 1))
  expect_equal(nrow(fit$v), sum(vapply(X, ncol, integer(1L))))
})

test_that("mcca family on block sources matches the in-memory fits", {
  X <- .gm_shared_rows(seed = 4)
  fields <- c("v", "canonical_weights", "cor_loadings", "scaled_loadings_by_block", "partial_scores")
  fit_mem <- suppressMessages(mcca(X, ncomp = 3))
  fit_src <- suppressMessages(mcca(as_block_source(X, chunk_size = 9), ncomp = 3))
  # RSpectra's iterative eigensolver agrees to about 1e-6 across inputs.
  .gm_expect_same(fit_mem, fit_src, fields, tolerance = 1e-5)

  sim <- .gm_linked_rows(seed = 5)
  src <- as_block_source(sim$X, chunk_size = 13)
  for (normalization in c("MFA", "balanced")) {
    fit_mem <- suppressMessages(aligned_mcca(sim$X, sim$idx, N = sim$N, ncomp = 3,
                                             normalization = normalization))
    fit_src <- suppressMessages(aligned_mcca(src, sim$idx, N = sim$N, ncomp = 3,
                                             normalization = normalization))
    # The balanced weights come from an iteration that stops at tol = 1e-6.
    .gm_expect_same(fit_mem, fit_src, fields, tolerance = 1e-4)
  }

  fit_mem <- suppressMessages(anchored_mcca(sim$Y, sim$X, sim$idx, ncomp = 3))
  fit_src <- suppressMessages(anchored_mcca(sim$Y, src, sim$idx, ncomp = 3))
  .gm_expect_same(fit_mem, fit_src, fields, tolerance = 1e-6)
  refit <- suppressMessages(fit_src$fit_spec$refit$fit_fn(fit_src$fit_spec$refit$data))
  .gm_expect_same(fit_src, refit, tolerance = 1e-8)
})

test_that("anchored MFA variants on block sources match the in-memory fits", {
  sim <- .gm_linked_rows(seed = 6)
  src <- as_block_source(sim$X, chunk_size = 13)
  fields <- c("v", "V_list", "cor_loadings", "partial_scores", "block_fit")
  ctrl <- list(ncomp = sim$K, max_iter = 200, tol = 1e-12, score_constraint = "orthonormal")

  fit_mem <- suppressMessages(do.call(anchored_mfa, c(list(sim$Y, sim$X, sim$idx), ctrl)))
  fit_src <- suppressMessages(do.call(anchored_mfa, c(list(sim$Y, src, sim$idx), ctrl)))
  .gm_expect_same(fit_mem, fit_src, fields)
  expect_equal(
    predict(fit_src, sim$X$B2[1:4, ], block = "B2"),
    predict(fit_mem, sim$X$B2[1:4, ], block = "B2"),
    tolerance = 1e-8
  )
  expect_equal(suppressMessages(linked_mfa(sim$Y, src, sim$idx, ncomp = sim$K))$strategy, "gram")
  boot <- suppressMessages(infer_muscal(fit_src, method = "bootstrap", nrep = 2, seed = 2))
  expect_equal(boot$n_success, 2L)

  graph_ctrl <- c(ctrl, list(score_graph_lambda = 1, score_graph_k = 5))
  for (fn in list(graph_anchored_mfa, coupled_graph_anchored_mfa)) {
    fit_mem <- suppressMessages(do.call(fn, c(list(sim$Y, sim$X, sim$idx), graph_ctrl)))
    fit_src <- suppressMessages(do.call(fn, c(list(sim$Y, src, sim$idx), graph_ctrl)))
    .gm_expect_same(fit_mem, fit_src, fields)
  }

  expect_error(
    anchored_mfa(sim$Y, src, sim$idx, ncomp = 2, feature_groups = "colnames", feature_lambda = 1),
    "not supported for block sources"
  )
})

test_that("response_aligned_mfa on a block source matches the in-memory fit", {
  sim <- .gm_linked_rows(seed = 7)
  Y <- lapply(sim$X, function(x) x[, 1:3] %*% matrix(rnorm(9), 3) + rnorm(nrow(x) * 3))
  src <- as_block_source(sim$X, chunk_size = 11)
  fit_mem <- suppressMessages(response_aligned_mfa(Y, sim$X, ncomp = sim$K, max_iter = 200, tol = 1e-12))
  fit_src <- suppressMessages(response_aligned_mfa(Y, src, ncomp = sim$K, max_iter = 200, tol = 1e-12))
  .gm_expect_same(fit_mem, fit_src, c("v", "V_list", "cor_loadings", "partial_scores", "B"),
                  tolerance = 1e-5)
  expect_error(
    response_aligned_mfa(Y, src, ncomp = 2, feature_graph = "colnames", graph_lambda = 1),
    "not supported for block sources"
  )
})

test_that("aligned_interbattery compresses shared-loading domains as stacked blocks", {
  set.seed(11)
  K <- 2
  n <- 40
  mk <- function(Z, W) Z %*% W + matrix(rnorm(nrow(Z) * ncol(W), sd = 0.1), nrow(Z))
  Zs <- list(S1 = matrix(rnorm(n * K), n), S2 = matrix(rnorm(n * K), n))
  Wx <- list(fmri = matrix(rnorm(K * 60), K), eeg = matrix(rnorm(K * 30), K))
  Wy <- matrix(rnorm(K * 8), K)
  X <- list(S1_fmri = mk(Zs$S1, Wx$fmri), S1_eeg = mk(Zs$S1, Wx$eeg),
            S2_fmri = mk(Zs$S2, Wx$fmri), S2_eeg = mk(Zs$S2, Wx$eeg))
  Y <- list(S1_beh = mk(Zs$S1, Wy), S2_beh = mk(Zs$S2, Wy))
  xi <- data.frame(block = names(X), subject = c("S1", "S1", "S2", "S2"),
                   domain = c("fmri", "eeg", "fmri", "eeg"))
  yi <- data.frame(block = names(Y), subject = c("S1", "S2"), domain = "beh")
  fit <- function(XX, YY) {
    suppressMessages(aligned_interbattery(XX, YY, x_block_info = xi, y_block_info = yi,
                                          ncomp = K, max_iter = 500, tol = 1e-10))
  }
  fit_mem <- fit(X, Y)
  fit_x <- fit(as_block_source(X, chunk_size = 11), Y)
  fit_xy <- fit(as_block_source(X, chunk_size = 7), as_block_source(Y, chunk_size = 3))

  for (f in list(fit_x, fit_xy)) {
    .gm_expect_same(fit_mem, f, c("v", "x_loadings", "y_loadings", "block_fit"), tolerance = 1e-8)
    expect_equal(names(f$x_loadings), c("fmri", "eeg"))
    expect_equal(f$block_fit$p, fit_mem$block_fit$p)
    expect_equal(f$x_block_info, fit_mem$x_block_info)
  }
  # The stacked fMRI domain has 2n rows (centered: rank 2n - 1) and 60 features.
  expect_equal(fit_xy$gram_rank[["x::fmri"]], min(2L * n - 1L, 60L))
  expect_error(
    aligned_interbattery(as_block_source(X), Y, x_block_info = xi, y_block_info = yi,
                         ncomp = K, gram_rank = 5),
    "not supported"
  )
})
