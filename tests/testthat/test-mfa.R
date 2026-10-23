library(testthat)
library(muscal)

# Basic MFA run should produce a valid object

.sim_mfa_blocks <- function(n, k, p_vec, noise_sd = 0) {
  Z <- scale(matrix(rnorm(n * k), nrow = n, ncol = k), center = TRUE, scale = FALSE)
  blocks <- lapply(p_vec, function(p) {
    A <- matrix(rnorm(p * k), nrow = p, ncol = k)
    signal <- Z %*% t(A)
    signal + matrix(rnorm(n * p, sd = noise_sd), nrow = n, ncol = p)
  })
  list(blocks = blocks, Z = Z)
}

.subspace_rel_diff <- function(S1, S2) {
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
}

.mfa_mean_init <- function(blocks) {
  lapply(blocks, function(block) {
    out <- block
    cm <- colMeans(out, na.rm = TRUE)
    miss <- is.na(out)
    if (any(miss)) {
      out[miss] <- cm[col(out)[miss]]
    }
    out
  })
}

.mfa_mask_random <- function(block, prop) {
  mask <- matrix(FALSE, nrow(block), ncol(block))
  nmiss <- max(1L, floor(prop * length(block)))
  mask[sample.int(length(block), nmiss)] <- TRUE
  all_missing_cols <- which(colSums(!mask) == 0L)
  if (length(all_missing_cols) > 0L) {
    mask[1L, all_missing_cols] <- FALSE
  }
  mask
}

.mfa_missing_values <- function(blocks, masks) {
  unlist(Map(function(x, m) x[m], blocks, masks), use.names = FALSE)
}

.mfa_assert_observed_preserved <- function(original, imputed) {
  for (b in seq_along(original)) {
    obs <- !is.na(original[[b]])
    expect_equal(imputed[[b]][obs], original[[b]][obs])
  }
}

.mfa_diag_values <- function(x) {
  if (inherits(x, "Matrix")) {
    return(Matrix::diag(x))
  }
  diag(x)
}

test_that("mfa.list produces expected structure", {
  set.seed(1)
  blocks <- replicate(3, matrix(rnorm(20), nrow = 5), simplify = FALSE)
  res <- mfa(blocks, ncomp = 2)
  expect_s3_class(res, "mfa")
  expect_equal(length(res$block_indices), 3)
  ncols <- sapply(blocks, ncol)
  expect_equal(res$block_indices[[1]], 1:ncols[1])
  expect_equal(res$block_indices[[2]], (ncols[1]+1):(ncols[1]+ncols[2]))
  expect_equal(res$block_indices[[3]], (ncols[1]+ncols[2]+1):sum(ncols))
  expect_length(res$alpha, 3)
  expect_true(all(res$alpha > 0))
  expect_equal(multivarious::ncomp(res), 2)
})

# Error conditions in mfa.multiblock

test_that("mfa.multiblock validates input", {
  x1 <- matrix(rnorm(10), nrow = 5)
  expect_error(mfa(list(x1)))

  x2 <- matrix(rnorm(12), nrow = 6)
  expect_error(mfa(list(x1, x2)))

  expect_error(mfa(list(x1, x1), normalization = "custom"))
})

# Normalization factor computation for MFA mode

test_that("normalization_factors MFA matches svd", {
  X1 <- matrix(1:4, nrow = 2)
  X2 <- matrix(2:5, nrow = 2)
  nf <- muscal:::normalization_factors(list(X1, X2), type = "MFA")
  expect_equal(nf[1], 1/(svd(X1)$d[1]^2))
  expect_equal(nf[2], 1/(svd(X2)$d[1]^2))
})

test_that("mfa exposes a standard out-of-sample contract", {
  set.seed(4)
  blocks <- list(
    X1 = matrix(rnorm(60), nrow = 10),
    X2 = matrix(rnorm(50), nrow = 10)
  )
  fit <- mfa(blocks, ncomp = 2)

  expect_equal(fit$task, "reconstruction")
  expect_equal(fit$fit_spec$method, "mfa")
  expect_true(fit$fit_spec$refit_supported)
  expect_true(is.list(fit$fit_spec$refit))
  expect_setequal(fit$oos_types, c("scores", "reconstruction"))

  Xnew <- do.call(cbind, blocks) + matrix(rnorm(110, sd = 0.05), nrow = 10)
  scores_new <- predict(fit, Xnew, type = "scores")
  recon_new <- predict(fit, Xnew, type = "reconstruction")

  expect_equal(scores_new, multivarious::project(fit, Xnew))
  expect_equal(dim(scores_new), c(nrow(Xnew), multivarious::ncomp(fit)))
  expect_equal(dim(recon_new), dim(Xnew))
  expect_equal(predict(fit, type = "scores"), multivarious::scores(fit))
  expect_error(predict(fit, type = "reconstruction"), "`new_data` must be supplied")
})

test_that("mfa refit contract can refit the stored training data", {
  set.seed(41)
  blocks <- list(
    X1 = matrix(rnorm(48), nrow = 12),
    X2 = matrix(rnorm(36), nrow = 12)
  )
  fit <- mfa(blocks, ncomp = 2)

  refit <- fit$fit_spec$refit$fit_fn(fit$fit_spec$refit$data)

  expect_s3_class(refit, "mfa")
  expect_equal(dim(refit$v), dim(fit$v))
  expect_equal(length(refit$block_indices), length(fit$block_indices))
  expect_true(all(is.finite(refit$sdev)))
})

test_that("mfa refuses missing values by default", {
  set.seed(42)
  blocks <- list(
    X1 = matrix(rnorm(48), nrow = 12),
    X2 = matrix(rnorm(36), nrow = 12)
  )
  blocks$X1[1, 1] <- NA_real_

  expect_error(
    mfa(blocks, ncomp = 2),
    "does not accept NA values"
  )
})

test_that("mfa missing path is identical when there are no NAs", {
  set.seed(421)
  blocks <- list(
    X1 = matrix(rnorm(72), nrow = 12),
    X2 = matrix(rnorm(60), nrow = 12),
    X3 = matrix(rnorm(48), nrow = 12)
  )

  fit0 <- suppressMessages(mfa(blocks, ncomp = 3, normalization = "MFA"))
  fit1 <- suppressMessages(mfa(blocks, ncomp = 3, normalization = "MFA", missing = "regularized"))

  expect_equal(fit1$sdev, fit0$sdev, tolerance = 1e-8)
  signs <- sign(colSums(fit1$v[, seq_len(3), drop = FALSE] * fit0$v[, seq_len(3), drop = FALSE]))
  signs[signs == 0] <- 1
  expect_equal(
    sweep(fit1$v[, seq_len(3), drop = FALSE], 2, signs, "*"),
    fit0$v[, seq_len(3), drop = FALSE],
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
})

test_that("mfa regularized missing path imputes and preserves observed cells", {
  set.seed(422)
  blocks <- list(
    X1 = matrix(rnorm(80), nrow = 16),
    X2 = matrix(rnorm(64), nrow = 16)
  )
  observed <- blocks
  blocks$X1[cbind(c(1, 3, 7), c(1, 2, 4))] <- NA_real_
  blocks$X2[5, ] <- NA_real_

  fit <- suppressMessages(
    mfa(
      blocks,
      ncomp = 2,
      missing = "regularized",
      ncp_impute = 2,
      missing_maxiter = 100,
      return_imputed = TRUE
    )
  )

  expect_s3_class(fit, "mfa")
  expect_true(is.list(fit$missing))
  expect_equal(fit$missing$method, "regularized")
  expect_lte(fit$missing$iterations, 100)
  expect_true(is.finite(fit$missing$delta))
  expect_false(anyNA(do.call(cbind, fit$imputed_data)))
  expect_true(anyNA(do.call(cbind, fit$fit_spec$refit$data)))

  for (b in seq_along(blocks)) {
    obs <- !is.na(blocks[[b]])
    expect_equal(fit$imputed_data[[b]][obs], observed[[b]][obs])
  }

  refit <- suppressMessages(fit$fit_spec$refit$fit_fn(fit$fit_spec$refit$data))
  expect_s3_class(refit, "mfa")
  expect_equal(refit$missing$method, "regularized")
})

test_that("mfa missing path rejects all-missing variables", {
  set.seed(423)
  blocks <- list(
    X1 = matrix(rnorm(48), nrow = 12),
    X2 = matrix(rnorm(36), nrow = 12)
  )
  blocks$X1[, 2] <- NA_real_

  expect_error(
    mfa(blocks, ncomp = 2, missing = "em"),
    "all values missing"
  )
})

test_that("mfa missing controls validate the imputation contract", {
  set.seed(424)
  blocks <- list(
    X1 = matrix(rnorm(48), nrow = 12),
    X2 = matrix(rnorm(36), nrow = 12)
  )
  blocks$X1[1, 1] <- NA_real_

  expect_error(
    mfa(blocks, ncomp = 2, missing = "regularized", ncp_impute = 0),
    "`ncp_impute`"
  )
  expect_error(
    mfa(blocks, ncomp = 2, missing = "regularized", missing_tol = 0),
    "`missing_tol`"
  )
  expect_error(
    mfa(blocks, ncomp = 2, missing = "regularized", missing_maxiter = 0),
    "`missing_maxiter`"
  )
  expect_error(
    mfa(blocks, ncomp = 2, missing = "regularized", return_imputed = NA),
    "`return_imputed`"
  )
})

test_that("mfa EM imputation recovers noiseless low-rank held-out entries", {
  set.seed(425)
  sim <- .sim_mfa_blocks(n = 50, k = 2, p_vec = c(6, 7, 5), noise_sd = 0)
  blocks <- sim$blocks
  names(blocks) <- paste0("X", seq_along(blocks))

  set.seed(426)
  masks <- lapply(blocks, .mfa_mask_random, prop = 0.12)
  missing_blocks <- blocks
  for (b in seq_along(missing_blocks)) {
    missing_blocks[[b]][masks[[b]]] <- NA_real_
  }

  mean_init <- .mfa_mean_init(missing_blocks)
  fit <- suppressMessages(
    mfa(
      missing_blocks,
      ncomp = 2,
      preproc = multivarious::pass(),
      normalization = "None",
      missing = "em",
      ncp_impute = 2,
      missing_tol = 1e-8,
      missing_maxiter = 200,
      return_imputed = TRUE
    )
  )

  truth <- .mfa_missing_values(blocks, masks)
  estimated <- .mfa_missing_values(fit$imputed_data, masks)
  baseline <- .mfa_missing_values(mean_init, masks)
  rmse_est <- sqrt(mean((estimated - truth)^2))
  rmse_baseline <- sqrt(mean((baseline - truth)^2))

  expect_true(fit$missing$converged)
  expect_lt(rmse_est, 1e-5)
  expect_lt(rmse_est, rmse_baseline * 1e-4)
  .mfa_assert_observed_preserved(missing_blocks, fit$imputed_data)
})

test_that("mfa missing imputation is row-permutation equivariant", {
  set.seed(427)
  blocks <- list(
    X1 = matrix(rnorm(120), nrow = 20),
    X2 = matrix(rnorm(100), nrow = 20)
  )
  missing_blocks <- blocks
  missing_blocks$X1[cbind(c(2, 5, 9, 13), c(1, 3, 5, 2))] <- NA_real_
  missing_blocks$X2[cbind(c(1, 7, 12), c(2, 4, 1))] <- NA_real_

  perm <- sample.int(nrow(blocks$X1))
  inv_perm <- order(perm)
  permuted <- lapply(missing_blocks, function(block) block[perm, , drop = FALSE])

  fit <- suppressMessages(
    mfa(
      missing_blocks,
      ncomp = 2,
      missing = "em",
      ncp_impute = 2,
      missing_tol = 1e-8,
      missing_maxiter = 100,
      return_imputed = TRUE
    )
  )
  fit_perm <- suppressMessages(
    mfa(
      permuted,
      ncomp = 2,
      missing = "em",
      ncp_impute = 2,
      missing_tol = 1e-8,
      missing_maxiter = 100,
      return_imputed = TRUE
    )
  )

  for (b in seq_along(blocks)) {
    expect_equal(
      fit_perm$imputed_data[[b]][inv_perm, , drop = FALSE],
      fit$imputed_data[[b]],
      tolerance = 1e-8
    )
  }
  expect_equal(fit_perm$missing$iterations, fit$missing$iterations)
})

test_that("mfa regularized missing path is finite across randomized masks", {
  for (seed in 428:432) {
    set.seed(seed)
    blocks <- list(
      X1 = matrix(rnorm(72), nrow = 18),
      X2 = matrix(rnorm(90), nrow = 18),
      X3 = matrix(rnorm(54), nrow = 18)
    )
    masks <- lapply(blocks, .mfa_mask_random, prop = 0.10)
    missing_blocks <- blocks
    for (b in seq_along(missing_blocks)) {
      missing_blocks[[b]][masks[[b]]] <- NA_real_
    }

    fit <- suppressMessages(
      mfa(
        missing_blocks,
        ncomp = 2,
        missing = "regularized",
        ncp_impute = 2,
        missing_tol = 1e-5,
        missing_maxiter = 100,
        return_imputed = TRUE
      )
    )

    expect_equal(fit$missing$mask, masks)
    expect_lte(fit$missing$iterations, 100)
    expect_true(is.finite(fit$missing$delta))
    expect_true(all(is.finite(fit$sdev)))
    expect_true(all(is.finite(fit$s)))
    expect_true(all(is.finite(fit$v)))
    expect_false(anyNA(do.call(cbind, fit$imputed_data)))
    .mfa_assert_observed_preserved(missing_blocks, fit$imputed_data)

    recon <- predict(fit, do.call(cbind, fit$imputed_data), type = "reconstruction")
    expect_true(all(is.finite(recon)))
  }
})

test_that("mfa missing paths handle custom positive row weights", {
  set.seed(433)
  blocks <- list(
    X1 = matrix(rnorm(72), nrow = 12),
    X2 = matrix(rnorm(60), nrow = 12)
  )
  missing_blocks <- blocks
  missing_blocks$X1[cbind(c(1, 3, 8), c(1, 4, 2))] <- NA_real_
  missing_blocks$X2[cbind(c(4, 7, 10), c(3, 5, 1))] <- NA_real_

  row_weights <- seq(0.6, 1.8, length.out = nrow(blocks$X1))

  for (method in c("regularized", "em")) {
    fit_vec <- suppressMessages(
      mfa(
        missing_blocks,
        ncomp = 2,
        normalization = "custom",
        M = row_weights,
        missing = method,
        ncp_impute = 2,
        missing_tol = 1e-6,
        missing_maxiter = 100,
        return_imputed = TRUE
      )
    )
    fit_diag <- suppressMessages(
      mfa(
        missing_blocks,
        ncomp = 2,
        normalization = "custom",
        M = diag(row_weights),
        missing = method,
        ncp_impute = 2,
        missing_tol = 1e-6,
        missing_maxiter = 100,
        return_imputed = TRUE
      )
    )

    expect_s3_class(fit_vec, "mfa")
    expect_equal(fit_vec$missing$method, method)
    expect_equal(dim(fit_vec$M_genpca), c(length(row_weights), length(row_weights)))
    expect_equal(.mfa_diag_values(fit_vec$M_genpca), row_weights, tolerance = 1e-12)
    expect_equal(fit_vec$sdev, fit_diag$sdev, tolerance = 1e-10)
    expect_equal(fit_vec$imputed_data, fit_diag$imputed_data, tolerance = 1e-10)
    expect_false(anyNA(do.call(cbind, fit_vec$imputed_data)))
    expect_true(all(is.finite(fit_vec$sdev)))
    expect_true(all(is.finite(fit_vec$s)))
    expect_true(all(is.finite(fit_vec$v)))
    .mfa_assert_observed_preserved(missing_blocks, fit_vec$imputed_data)

    refit <- suppressMessages(fit_vec$fit_spec$refit$fit_fn(fit_vec$fit_spec$refit$data))
    expect_equal(refit$missing$method, method)
    expect_equal(.mfa_diag_values(refit$M_genpca), row_weights, tolerance = 1e-12)
  }
})

test_that("mfa row weights affect missing-data fits instead of being ignored", {
  set.seed(434)
  blocks <- list(
    X1 = matrix(rnorm(72), nrow = 12),
    X2 = matrix(rnorm(60), nrow = 12)
  )
  missing_blocks <- blocks
  missing_blocks$X1[cbind(c(2, 6, 11), c(2, 5, 1))] <- NA_real_
  missing_blocks$X2[cbind(c(1, 5, 9), c(4, 2, 5))] <- NA_real_

  fit_uniform <- suppressMessages(
    mfa(
      missing_blocks,
      ncomp = 2,
      normalization = "custom",
      M = rep(1, nrow(blocks$X1)),
      missing = "em",
      ncp_impute = 2,
      missing_tol = 1e-6,
      missing_maxiter = 30,
      return_imputed = TRUE
    )
  )
  fit_weighted <- suppressMessages(
    mfa(
      missing_blocks,
      ncomp = 2,
      normalization = "custom",
      M = seq(0.5, 2.0, length.out = nrow(blocks$X1)),
      missing = "em",
      ncp_impute = 2,
      missing_tol = 1e-6,
      missing_maxiter = 30,
      return_imputed = TRUE
    )
  )

  expect_gt(max(abs(fit_weighted$sdev - fit_uniform$sdev)), 1e-4)
  expect_gt(
    max(abs(do.call(cbind, fit_weighted$imputed_data) - do.call(cbind, fit_uniform$imputed_data))),
    1e-4
  )
  .mfa_assert_observed_preserved(missing_blocks, fit_weighted$imputed_data)
})

test_that("mfa missing path tolerates zero row weights in custom M", {
  set.seed(435)
  blocks <- list(
    X1 = matrix(rnorm(80), nrow = 16),
    X2 = matrix(rnorm(64), nrow = 16)
  )
  missing_blocks <- blocks
  missing_blocks$X1[2, ] <- NA_real_
  missing_blocks$X2[cbind(c(3, 9), c(1, 4))] <- NA_real_

  row_weights <- rep(1, nrow(blocks$X1))
  row_weights[2] <- 0

  fit <- suppressMessages(
    mfa(
      missing_blocks,
      ncomp = 2,
      normalization = "custom",
      M = row_weights,
      missing = "regularized",
      ncp_impute = 2,
      missing_tol = 1e-6,
      missing_maxiter = 30,
      return_imputed = TRUE
    )
  )

  expect_s3_class(fit, "mfa")
  expect_equal(.mfa_diag_values(fit$M_genpca), row_weights, tolerance = 1e-12)
  expect_true(all(is.finite(fit$sdev)))
  expect_true(all(is.finite(fit$s)))
  expect_true(all(is.finite(fit$v)))
  expect_false(anyNA(do.call(cbind, fit$imputed_data)))
  .mfa_assert_observed_preserved(missing_blocks, fit$imputed_data)
})

test_that("mfa exactly recovers noiseless low-rank structure and reconstruction improves with ncomp", {
  set.seed(43)
  sim <- .sim_mfa_blocks(n = 60, k = 3, p_vec = c(10, 12, 8), noise_sd = 0)
  X_concat <- do.call(cbind, sim$blocks)

  fits <- lapply(1:3, function(k) {
    suppressMessages(
      mfa(
        sim$blocks,
        ncomp = k,
        preproc = multivarious::pass(),
        normalization = "None"
      )
    )
  })

  mse <- vapply(fits, function(fit) {
    Xhat <- predict(fit, X_concat, type = "reconstruction")
    mean((Xhat - X_concat)^2)
  }, numeric(1))

  expect_true(all(diff(mse) <= 1e-10))
  expect_lt(mse[[3]], 1e-10)

  S_fit <- multivarious::scores(fits[[3]])
  P_fit <- S_fit %*% solve(crossprod(S_fit), t(S_fit))
  P_true <- sim$Z %*% solve(crossprod(sim$Z), t(sim$Z))
  rel <- norm(P_fit - P_true, type = "F") / (norm(P_true, type = "F") + 1e-12)
  expect_lt(rel, 1e-8)
})

test_that("mfa latent recovery degrades as noise increases", {
  set.seed(44)
  low_noise <- .sim_mfa_blocks(n = 70, k = 2, p_vec = c(9, 11, 10), noise_sd = 0.02)
  high_noise <- .sim_mfa_blocks(n = 70, k = 2, p_vec = c(9, 11, 10), noise_sd = 0.35)

  fit_low <- suppressMessages(
    mfa(low_noise$blocks, ncomp = 2, preproc = multivarious::pass(), normalization = "None")
  )
  fit_high <- suppressMessages(
    mfa(high_noise$blocks, ncomp = 2, preproc = multivarious::pass(), normalization = "None")
  )

  rel_low <- .subspace_rel_diff(multivarious::scores(fit_low), low_noise$Z)
  rel_high <- .subspace_rel_diff(multivarious::scores(fit_high), high_noise$Z)

  expect_lt(rel_low, rel_high)
  expect_lt(rel_low, 0.15)
})

# Classical regularized iterative MFA on centered blocks (Josse & Husson, 2012),
# written directly with base::svd as an independent reference.
.ref_impute_mfa <- function(blocks, ncp, method, tol = 1e-12, maxiter = 5000) {
  X <- do.call(cbind, blocks)
  miss <- is.na(X)
  grp <- rep(seq_along(blocks), vapply(blocks, ncol, integer(1L)))
  X[miss] <- colMeans(X, na.rm = TRUE)[col(X)[miss]]
  n <- nrow(X)
  p <- ncol(X)
  for (it in seq_len(maxiter)) {
    mu <- colMeans(X)
    Xc <- sweep(X, 2L, mu)
    w <- vapply(seq_along(blocks), function(b) {
      svd(Xc[, grp == b, drop = FALSE], nu = 0, nv = 0)$d[[1L]]
    }, numeric(1L))
    s <- svd(sweep(Xc, 2L, w[grp], `/`))
    d <- s$d
    sigma2 <- if (identical(method, "em")) {
      0
    } else {
      min(
        n * p / min(p, n - 1) * sum(d[-seq_len(ncp)]^2) / ((n - 1 - ncp) * (p - ncp)),
        d[[ncp + 1L]]^2
      )
    }
    keep <- seq_len(ncp)
    fitted <- s$u[, keep, drop = FALSE] %*%
      diag((d[keep]^2 - sigma2) / d[keep], nrow = ncp) %*%
      t(s$v[, keep, drop = FALSE])
    fitted <- sweep(sweep(fitted, 2L, w[grp], `*`), 2L, mu, `+`)
    delta <- max(abs(fitted[miss] - X[miss]))
    X[miss] <- fitted[miss]
    if (delta < tol) break
  }
  X
}

.sim_mfa_scaled_blocks <- function(n, k, p_vec, scales, shifts, noise_sd = 0) {
  Z <- matrix(rnorm(n * k), nrow = n, ncol = k)
  blocks <- Map(function(p, sc, sh) {
    signal <- Z %*% matrix(rnorm(k * p), nrow = k, ncol = p)
    (signal + matrix(rnorm(n * p, sd = noise_sd), nrow = n, ncol = p)) * sc + sh
  }, p_vec, scales, shifts)
  names(blocks) <- paste0("X", seq_along(blocks))
  blocks
}

.mfa_std_rmse <- function(estimate, truth_blocks, masks) {
  truth <- .mfa_missing_values(truth_blocks, masks)
  scale <- unlist(Map(function(b, m) apply(b, 2L, stats::sd)[col(m)[m]], truth_blocks, masks),
                  use.names = FALSE)
  sqrt(mean(((estimate - truth) / scale)^2))
}

test_that("mfa reconstruction back-projects in the data metric under MFA block weights", {
  set.seed(601)
  blocks <- list(
    X1 = matrix(rnorm(100, sd = 3), nrow = 20),
    X2 = matrix(rnorm(100, sd = 30), nrow = 20),
    X3 = matrix(rnorm(100, sd = 0.3), nrow = 20)
  )
  Xc <- do.call(cbind, blocks)
  fit <- suppressMessages(mfa(blocks, ncomp = 2))
  gsvd_rank2 <- sweep(
    as.matrix(fit$ou) %*% diag(fit$sdev, nrow = 2) %*% t(as.matrix(fit$ov)),
    2L, colMeans(Xc), `+`
  )

  expect_equal(multivarious::reconstruct(fit), gsvd_rank2, tolerance = 1e-10)
  expect_equal(
    multivarious::reconstruct(fit, rowind = 3:7, colind = 4:12),
    gsvd_rank2[3:7, 4:12],
    tolerance = 1e-10
  )
  expect_equal(predict(fit, Xc, type = "reconstruction"), gsvd_rank2, tolerance = 1e-10)

  fit_full <- suppressMessages(mfa(blocks, ncomp = 15))
  expect_equal(predict(fit_full, Xc, type = "reconstruction"), Xc, tolerance = 1e-8)
})

test_that("mfa EM imputation recovers noiseless low-rank entries under MFA block weighting", {
  set.seed(602)
  blocks <- .sim_mfa_scaled_blocks(
    n = 50, k = 2, p_vec = c(6, 7, 5),
    scales = c(1, 100, 0.01), shifts = c(0, 500, -2)
  )
  masks <- lapply(blocks, .mfa_mask_random, prop = 0.12)
  masks$X2[3:5, ] <- TRUE
  missing_blocks <- Map(function(b, m) { b[m] <- NA_real_; b }, blocks, masks)

  fit <- suppressMessages(
    mfa(
      missing_blocks,
      ncomp = 2,
      missing = "em",
      ncp_impute = 2,
      missing_tol = 1e-10,
      missing_maxiter = 2000,
      return_imputed = TRUE
    )
  )
  estimated <- .mfa_missing_values(fit$imputed_data, masks)
  baseline <- .mfa_missing_values(.mfa_mean_init(missing_blocks), masks)

  expect_true(fit$missing$converged)
  expect_lt(.mfa_std_rmse(estimated, blocks, masks), 1e-5)
  expect_lt(.mfa_std_rmse(estimated, blocks, masks), 1e-4 * .mfa_std_rmse(baseline, blocks, masks))
  .mfa_assert_observed_preserved(missing_blocks, fit$imputed_data)
})

test_that("mfa iterative imputation matches a reference regularized iterative MFA", {
  set.seed(603)
  blocks <- .sim_mfa_scaled_blocks(
    n = 40, k = 2, p_vec = c(5, 6, 4),
    scales = c(1, 20, 0.1), shifts = c(0, 10, 3), noise_sd = 0.5
  )
  masks <- lapply(blocks, .mfa_mask_random, prop = 0.15)
  masks$X3[7:9, ] <- TRUE
  missing_blocks <- Map(function(b, m) { b[m] <- NA_real_; b }, blocks, masks)
  all_masks <- do.call(cbind, masks)

  for (method in c("regularized", "em")) {
    fit <- suppressMessages(
      mfa(
        missing_blocks,
        ncomp = 2,
        missing = method,
        ncp_impute = 2,
        missing_tol = 1e-11,
        missing_maxiter = 5000,
        return_imputed = TRUE
      )
    )
    reference <- .ref_impute_mfa(missing_blocks, ncp = 2, method = method)
    scale <- apply(do.call(cbind, blocks), 2L, stats::sd)[col(all_masks)[all_masks]]

    expect_true(fit$missing$converged)
    expect_lt(
      max(abs(do.call(cbind, fit$imputed_data)[all_masks] - reference[all_masks]) / scale),
      1e-6
    )
  }
})

test_that("mfa regularized imputation beats mean imputation on noisy low-rank data", {
  set.seed(604)
  blocks <- .sim_mfa_scaled_blocks(
    n = 60, k = 2, p_vec = c(8, 6, 10),
    scales = c(1, 50, 0.01), shifts = c(0, 100, -3), noise_sd = 0.5
  )
  masks <- lapply(blocks, .mfa_mask_random, prop = 0.2)
  missing_blocks <- Map(function(b, m) { b[m] <- NA_real_; b }, blocks, masks)

  fit <- suppressMessages(
    mfa(missing_blocks, ncomp = 2, missing = "regularized", return_imputed = TRUE)
  )
  rmse_reg <- .mfa_std_rmse(.mfa_missing_values(fit$imputed_data, masks), blocks, masks)
  rmse_mean <- .mfa_std_rmse(
    .mfa_missing_values(.mfa_mean_init(missing_blocks), masks), blocks, masks
  )

  expect_true(fit$missing$converged)
  expect_lt(rmse_reg, 0.7 * rmse_mean)
})

test_that("mfa imputation is equivariant to rescaling and shifting a block", {
  set.seed(605)
  blocks <- .sim_mfa_scaled_blocks(
    n = 30, k = 2, p_vec = c(5, 5), scales = c(1, 1), shifts = c(0, 0), noise_sd = 0.3
  )
  masks <- lapply(blocks, .mfa_mask_random, prop = 0.15)
  missing_blocks <- Map(function(b, m) { b[m] <- NA_real_; b }, blocks, masks)
  shifted <- missing_blocks
  shifted$X2 <- shifted$X2 * 1000 + 250

  for (method in c("regularized", "em")) {
    fit <- suppressMessages(mfa(missing_blocks, ncomp = 2, missing = method, return_imputed = TRUE))
    fit_shifted <- suppressMessages(mfa(shifted, ncomp = 2, missing = method, return_imputed = TRUE))

    expect_equal(fit_shifted$missing$iterations, fit$missing$iterations)
    expect_equal(fit_shifted$imputed_data$X1, fit$imputed_data$X1, tolerance = 1e-8)
    expect_equal(fit_shifted$imputed_data$X2, fit$imputed_data$X2 * 1000 + 250, tolerance = 1e-8)
  }
})

test_that("mfa regularized shrinkage follows the residual-variance formula", {
  d <- c(10, 6, 2, 1, 0.5)
  n <- 20
  p <- 5
  total <- sum(d^2)
  shrink <- muscal:::.mfa_regularized_shrinkage(d, ncp = 2, n = n, p = p, total_inertia = total)
  sigma2 <- n * p / min(p, n - 1) * sum(d[3:5]^2) / ((n - 1 - 2) * (p - 2))
  sigma2 <- min(sigma2, d[[3]]^2)

  expect_equal(shrink, (d[1:2]^2 - sigma2) / d[1:2]^2)
  expect_true(all(shrink > 0 & shrink <= 1))
  # A large residual is capped at the next singular value.
  capped <- muscal:::.mfa_regularized_shrinkage(d, ncp = 2, n = n, p = p, total_inertia = 1e4)
  expect_equal(capped, (d[1:2]^2 - d[[3]]^2) / d[1:2]^2)
})

test_that("mfa missing path warns when imputation does not converge", {
  set.seed(606)
  blocks <- list(X1 = matrix(rnorm(60), nrow = 12), X2 = matrix(rnorm(48), nrow = 12))
  blocks$X1[cbind(c(1, 4, 9), c(1, 3, 5))] <- NA_real_

  expect_warning(
    fit <- suppressMessages(mfa(blocks, ncomp = 2, missing = "em", missing_maxiter = 1)),
    "did not converge"
  )
  expect_false(fit$missing$converged)
  expect_equal(fit$missing$iterations, 1L)
})

test_that("mfa missing path caps ncp_impute at the identifiable rank", {
  set.seed(607)
  blocks <- list(X1 = matrix(rnorm(20), nrow = 5), X2 = matrix(rnorm(20), nrow = 5))
  blocks$X1[2, 1] <- NA_real_

  expect_warning(
    fit <- suppressMessages(
      mfa(blocks, ncomp = 2, missing = "regularized", ncp_impute = 10, missing_maxiter = 500)
    ),
    "too large"
  )
  expect_true(fit$missing$converged)
  expect_equal(fit$missing$ncp_impute, 3L)
  # The refit contract stores the capped value, so refits do not warn again.
  expect_no_warning(
    refit <- suppressMessages(fit$fit_spec$refit$fit_fn(fit$fit_spec$refit$data))
  )
  expect_equal(refit$missing$ncp_impute, 3L)
})

test_that("mfa missing path flags rows with no observed data and rejects bad inputs", {
  set.seed(608)
  blocks <- list(X1 = matrix(rnorm(60), nrow = 12), X2 = matrix(rnorm(48), nrow = 12))
  empty <- blocks
  empty$X1[3, ] <- NA_real_
  empty$X2[3, ] <- NA_real_
  expect_warning(
    suppressMessages(mfa(empty, ncomp = 2, missing = "em", missing_maxiter = 500)),
    "no observed values"
  )

  infinite <- blocks
  infinite$X1[1, 1] <- NA_real_
  infinite$X2[2, 2] <- Inf
  expect_error(
    suppressMessages(mfa(infinite, ncomp = 2, missing = "em")),
    "infinite values"
  )

  ragged <- list(X1 = matrix(rnorm(60), nrow = 12), X2 = matrix(rnorm(44), nrow = 11))
  ragged$X1[1, 1] <- NA_real_
  expect_error(
    suppressMessages(muscal:::.mfa_impute_iterative(
      ragged, preproc = multivarious::center(), ncp = 2, normalization = "MFA",
      A = NULL, M = NULL, method = "em"
    )),
    "same number of rows"
  )
})
