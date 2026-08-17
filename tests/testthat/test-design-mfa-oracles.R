test_that("shared-grid weighted fit agrees with independent lm.wfit oracles", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = c(31, 37, 29, 43), p = 9, seed = 20)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    nuisance = ~ run,
    contrasts = list(saliency = "poly", probe = "sum", run = "sum"),
    store = "fitted"
  )

  X <- fit$design$X
  Y <- do.call(rbind, sim$Y)
  w <- unlist(fit$weights, use.names = FALSE)
  fitted_muscal <- X %*% fit$full_model$coefficient
  fitted_lm <- vapply(seq_len(ncol(Y)), function(j) {
    stats::lm.wfit(X, Y[, j], w = w)$fitted.values
  }, numeric(nrow(Y)))

  expect_equal(unname(fitted_muscal), unname(fitted_lm), tolerance = 2e-8)
  sse_lm <- sum((Y - fitted_lm)^2 * w)
  expect_equal(fit$full_model$sse, sse_lm, tolerance = 2e-8)

  saliency <- effects(fit, "saliency")
  reduced_names <- saliency$model_columns$reduced
  reduced <- X[, match(reduced_names, colnames(X)), drop = FALSE]
  reduced_fit <- vapply(seq_len(ncol(Y)), function(j) {
    stats::lm.wfit(reduced, Y[, j], w = w)$fitted.values
  }, numeric(nrow(Y)))
  partial_ss_lm <- sum((Y - reduced_fit)^2 * w) - sse_lm
  expect_equal(saliency$partial_ss, partial_ss_lm, tolerance = 2e-8)
  full_coef <- vapply(seq_len(ncol(Y)), function(j) {
    stats::lm.wfit(X, Y[, j], w = w)$coefficients
  }, numeric(ncol(X)))
  L_full <- fit$design$X_reference
  removed_names <- setdiff(saliency$model_columns$full, reduced_names)
  removed_columns <- match(removed_names, colnames(L_full))
  H_oracle <- L_full[, removed_columns, drop = FALSE] %*%
    full_coef[removed_columns, , drop = FALSE]
  H_oracle <- sweep(
    H_oracle,
    2L,
    colSums(H_oracle * saliency$reference_weights),
    "-"
  )
  expect_equal(unname(saliency$effect_matrix), unname(H_oracle), tolerance = 2e-8)
  expect_identical(fit$computation, "subject_streamed_sufficient")
})

test_that("subject-streamed sufficient statistics equal the stacked solver", {
  sim <- .sim_design_mfa(n_subjects = 5, trials = c(17, 23, 19, 29, 31), p = 11, seed = 201)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    nuisance = ~ run,
    contrasts = list(saliency = "poly", probe = "sum", run = "sum")
  )
  X <- fit$design$X
  rows <- split(seq_len(nrow(X)), rep(names(sim$Y), fit$n_trials))
  streamed <- muscal:::.dmfa_wls_list_sufficient(
    X,
    sim$Y,
    fit$weights,
    rows,
    fit$control$rank_tol
  )
  chunked <- muscal:::.dmfa_wls_list_sufficient(
    X,
    sim$Y,
    fit$weights,
    rows,
    fit$control$rank_tol,
    feature_chunk_size = 3L
  )
  stacked <- muscal:::.dmfa_wls_sufficient(
    X,
    do.call(rbind, sim$Y),
    unlist(fit$weights, use.names = FALSE),
    fit$control$rank_tol,
    subject_rows = rows
  )

  expect_equal(streamed$XtX, stacked$XtX, tolerance = 1e-13)
  expect_equal(streamed$XtY, stacked$XtY, tolerance = 1e-12)
  expect_equal(streamed$response_ss, stacked$response_ss, tolerance = 1e-12)
  expect_equal(streamed$coefficient, stacked$coefficient, tolerance = 1e-11)
  expect_equal(streamed$sse, stacked$sse, tolerance = 1e-10)
  expect_equal(chunked$coefficient, streamed$coefficient, tolerance = 1e-12)
  expect_equal(chunked$response_ss, streamed$response_ss, tolerance = 1e-12)
})

test_that("residual SSE is stable under large subject offsets", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = c(37, 41, 43, 47), p = 8, seed = 202)
  args <- list(
    formula = ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum")
  )
  baseline <- do.call(design_mfa, c(list(Y = sim$Y, design = sim$design), args))
  shifted_Y <- Map(function(y, offset) y + offset, sim$Y, 1e5 + seq_along(sim$Y))
  names(shifted_Y) <- names(sim$Y)
  shifted <- do.call(design_mfa, c(list(Y = shifted_Y, design = sim$design), args))

  expect_equal(shifted$full_model$sse, baseline$full_model$sse, tolerance = 2e-6)
  for (term in c("saliency", "probe")) {
    expect_equal(
      effects(shifted, term)$partial_ss,
      effects(baseline, term)$partial_ss,
      tolerance = 2e-6
    )
    expect_equal(
      effects(shifted, term)$partial_r2,
      effects(baseline, term)$partial_r2,
      tolerance = 2e-6
    )
  }
  expect_lt(effects(shifted, "probe")$partial_r2, 1)
})

test_that("thin effect SVD agrees with base SVD and reconstructs the effect", {
  set.seed(21)
  H <- matrix(rnorm(7 * 30), 7, 30)
  H <- sweep(H, 2, colMeans(H), "-")
  w <- seq_len(7)
  w <- w / sum(w)
  out <- muscal:::.dmfa_thin_svd(H, w, ncomp = 6, rank_tol = 1e-12)
  ref <- svd(H * sqrt(w), nu = 6, nv = 6)

  expect_equal(out$d, ref$d[1:6], tolerance = 1e-10)
  expect_equal(
    .dmfa_projector(out$v),
    .dmfa_projector(ref$v[, 1:6, drop = FALSE]),
    tolerance = 1e-9
  )

  sim <- .sim_design_mfa(n_subjects = 4, trials = 45, p = 10, noise_sd = 0.02, seed = 22)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    nuisance = ~ run,
    contrasts = list(saliency = "poly", probe = "sum", run = "sum"),
    store = "fitted"
  )
  effect <- effects(fit, "saliency")
  reconstruction <- effect$scores %*% t(effect$loadings)
  expect_equal(reconstruction, effect$effect_matrix, tolerance = 1e-8)
  expect_identical(effect$svd_path, "factored_small_gram")
})

test_that("factored thin SVD matches the dense path and stays on an r-by-r gram", {
  set.seed(2101)
  m <- 40L
  p <- 90L
  r <- 5L
  A <- matrix(rnorm(m * r), m, r)
  B <- matrix(rnorm(r * p), r, p)
  colnames(B) <- paste0("f", seq_len(p))
  w <- seq_len(m)
  w <- w / sum(w)
  A <- sweep(A, 2L, as.numeric(crossprod(w, A)), "-")
  H <- A %*% B

  dense <- muscal:::.dmfa_thin_svd(H, w, ncomp = 4, rank_tol = 1e-12)
  factored <- muscal:::.dmfa_thin_svd_factored(A, B, w, ncomp = 4, rank_tol = 1e-12)
  full <- muscal:::.dmfa_thin_svd_factored(A, B, w, rank_tol = 1e-12)

  expect_identical(factored$path, "factored_small_gram")
  expect_identical(factored$gram_dim, c(r, r))
  expect_equal(factored$d, dense$d, tolerance = 1e-10)
  expect_equal(factored$total_energy, dense$total_energy, tolerance = 1e-10)
  expect_equal(
    .dmfa_projector(factored$v),
    .dmfa_projector(dense$v),
    tolerance = 1e-9
  )
  expect_equal(full$scores %*% t(full$v), H, tolerance = 1e-8)
  expect_lt(sum(factored$d^2) / factored$total_energy, 1)

  rank_def <- cbind(A, A[, 1L])
  B_def <- rbind(B, 0)
  factored_def <- muscal:::.dmfa_thin_svd_factored(rank_def, B_def, w, rank_tol = 1e-12)
  expect_equal(factored_def$d, muscal:::.dmfa_thin_svd(H, w, rank_tol = 1e-12)$d, tolerance = 1e-9)
  expect_identical(factored_def$gram_dim, c(r + 1L, r + 1L))

  zero <- muscal:::.dmfa_thin_svd_factored(
    matrix(0, 6, 3),
    matrix(0, 3, 8),
    rep(1 / 6, 6),
    rank_tol = 1e-12
  )
  expect_length(zero$d, 0L)
  expect_true(is.finite(zero$total_energy))
  expect_equal(zero$total_energy, 0)
})

test_that("noiseless planted saliency field is recovered", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 90, p = 14, noise_sd = 0, seed = 23)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    nuisance = ~ run,
    contrasts = list(saliency = "poly", probe = "sum", run = "sum")
  )
  loading <- effect_map(fit, "saliency", 1)
  alignment <- abs(cor(loading, sim$truth$saliency))

  expect_gt(alignment, 0.999999)
  expect_gt(components(fit, "saliency")$explained[[1]], 0.999999)
  expect_true(all(is.finite(estimability(fit)$partial_r2)))
})

test_that("rank-truncated explained fractions use total effect energy", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 45, p = 10, noise_sd = 0.08, seed = 221)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    rank = c(saliency = 1, probe = 1),
    store = "fitted"
  )
  effect <- effects(fit, "saliency")
  total_energy <- sum(rowSums(effect$effect_matrix^2) * effect$reference_weights)

  expect_equal(effect$explained, effect$singular_values^2 / total_energy, tolerance = 1e-12)
  expect_lt(sum(effect$explained), 1 - 1e-6)
})

test_that("full-rank explained energy is complete in Euclidean and spatial metrics", {
  set.seed(204)
  H <- matrix(rnorm(7 * 9), 7, 9)
  weights <- seq_len(7)
  weights <- weights / sum(weights)
  H <- sweep(H, 2L, colSums(H * weights), "-")

  euclidean <- muscal:::.dmfa_thin_svd(
    H, weights, rank_tol = 1e-12
  )
  expect_equal(sum(euclidean$d^2) / euclidean$total_energy, 1, tolerance = 1e-11)

  adjacency <- Matrix::bandSparse(
    9,
    k = c(-1, 1),
    diagonals = list(rep(1, 8), rep(1, 8))
  )
  laplacian <- Matrix::Diagonal(x = Matrix::rowSums(adjacency)) - adjacency
  metric <- Matrix::Diagonal(9) + 2 * laplacian
  spatial <- muscal:::.dmfa_thin_svd(
    H, weights, rank_tol = 1e-12, feature_metric = metric
  )
  Hw <- H * sqrt(weights)
  energy_oracle <- sum(diag(Hw %*% solve(as.matrix(metric), t(Hw))))
  expect_equal(spatial$total_energy, energy_oracle, tolerance = 1e-10)
  expect_equal(sum(spatial$d^2) / spatial$total_energy, 1, tolerance = 1e-10)

  truncated <- muscal:::.dmfa_thin_svd(
    H, weights, ncomp = 2L, rank_tol = 1e-12, feature_metric = metric
  )
  expect_equal(
    sum(truncated$d^2) / truncated$total_energy,
    sum(spatial$d[1:2]^2) / energy_oracle,
    tolerance = 1e-11
  )
  expect_lt(sum(truncated$d^2) / truncated$total_energy, 1)

  zero <- muscal:::.dmfa_thin_svd(
    matrix(0, 4, 5), rep(0.25, 4), feature_metric = diag(5)
  )
  expect_length(zero$d, 0L)
  expect_true(is.finite(zero$total_energy))
  expect_equal(zero$total_energy, 0)
})

test_that("predictions and effect projection follow stored matrix contracts", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 35, p = 8, seed = 24)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    nuisance = ~ run,
    contrasts = list(saliency = "poly", probe = "sum", run = "sum")
  )

  pred <- predict(fit)
  L <- fit$design$X_reference
  expect_equal(pred, L %*% fit$full_model$coefficient, tolerance = 1e-10)
  expect_identical(colnames(pred), fit$feature_names)

  effect <- effects(fit, "probe")
  direct <- sim$Y[[1]] %*% effect$loadings
  projected <- multivarious::project(effect, sim$Y[[1]])
  expect_equal(projected, direct, tolerance = 1e-12)

  reordered <- sim$Y[[1]][, rev(colnames(sim$Y[[1]])), drop = FALSE]
  expect_equal(multivarious::project(effect, reordered), direct, tolerance = 1e-12)
  expect_error(multivarious::project(effect, sim$Y[[1]][, -1]), "does not match")
})
