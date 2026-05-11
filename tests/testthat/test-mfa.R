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
      missing_maxiter = 25,
      return_imputed = TRUE
    )
  )

  expect_s3_class(fit, "mfa")
  expect_true(is.list(fit$missing))
  expect_equal(fit$missing$method, "regularized")
  expect_lte(fit$missing$iterations, 25)
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
        missing_maxiter = 20,
        return_imputed = TRUE
      )
    )

    expect_equal(fit$missing$mask, masks)
    expect_lte(fit$missing$iterations, 20)
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
        missing_maxiter = 30,
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
        missing_maxiter = 30,
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
