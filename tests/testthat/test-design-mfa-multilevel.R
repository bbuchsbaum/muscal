test_that("ridge subject deviations agree with a dense equality-constrained oracle", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 28, p = 5, noise_sd = 0.04, seed = 50)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    nuisance = ~ run,
    contrasts = list(saliency = "poly", probe = "sum", run = "sum"),
    subject_model = "ridge_slopes",
    subject_penalty = 0.4,
    random_terms = "probe",
    control = design_mfa_control(multilevel_max_iter = 100, multilevel_tol = 1e-11)
  )

  X <- fit$design$X
  Y <- do.call(rbind, sim$Y)
  w <- unlist(fit$weights, use.names = FALSE)
  subject_weight <- unname(fit$subject_weights)
  random_columns <- fit$design$term_columns$probe
  r <- length(random_columns)
  S <- length(sim$Y)
  offsets <- cumsum(c(0L, vapply(sim$Y, nrow, integer(1))))
  Z <- matrix(0, nrow(X), S * r)
  for (s in seq_len(S)) {
    rows <- seq.int(offsets[[s]] + 1L, offsets[[s + 1L]])
    cols <- (s - 1L) * r + seq_len(r)
    Z[rows, cols] <- X[rows, random_columns, drop = FALSE]
  }
  A <- cbind(X, Z)
  penalty <- c(rep(0, ncol(X)), rep(0.4 * subject_weight, each = r))
  G <- crossprod(A * sqrt(w)) + diag(penalty)
  C <- matrix(0, r, ncol(A))
  for (s in seq_len(S)) {
    cols <- ncol(X) + (s - 1L) * r + seq_len(r)
    C[, cols] <- subject_weight[[s]] * diag(r)
  }
  kkt <- rbind(cbind(G, t(C)), cbind(C, matrix(0, r, r)))
  rhs <- rbind(crossprod(A, Y * w), matrix(0, r, ncol(Y)))
  oracle <- solve(kkt, rhs)
  coefficient_oracle <- oracle[seq_len(ncol(X)), , drop = FALSE]
  deviation_oracle <- oracle[ncol(X) + seq_len(ncol(Z)), , drop = FALSE]

  expect_equal(fit$full_model$coefficient, coefficient_oracle, tolerance = 2e-6)
  for (s in seq_len(S)) {
    rows <- (s - 1L) * r + seq_len(r)
    expect_equal(
      unname(subject_effects(fit, "probe")[[s]]),
      unname(deviation_oracle[rows, , drop = FALSE]),
      tolerance = 2e-6
    )
  }
  expect_true(all(diff(fit$full_model$objective_trace) <= 1e-10))
})

test_that("ridge deviations are centered and improve heterogeneous-subject fit", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 50, p = 8, noise_sd = 0.03, seed = 51)
  set.seed(52)
  heterogeneity_map <- scale(sin(seq(0, pi, length.out = 8)), center = TRUE, scale = FALSE)[, 1]
  slopes <- c(-1.2, -0.4, 0.5, 1.1)
  for (s in seq_along(sim$Y)) {
    sal <- as.numeric(sim$design[[s]]$saliency) - 3
    sim$Y[[s]] <- sim$Y[[s]] + slopes[[s]] * outer(sal, heterogeneity_map)
  }

  fixed <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    subject_model = "fixed_intercept"
  )
  pooled <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    subject_model = "ridge_slopes",
    subject_penalty = 0.02,
    random_terms = "saliency"
  )

  deviations <- subject_effects(pooled, "saliency")
  weighted_mean <- Reduce(`+`, Map(function(U, sw) sw * U, deviations, pooled$subject_weights))
  expect_lt(max(abs(weighted_mean)), 1e-10)
  expect_lt(pooled$full_model$sse, fixed$full_model$sse)
  expect_identical(effects(pooled, "saliency")$comparison, "penalized_objective")
  expect_true(is.finite(effects(pooled, "saliency")$objective_delta))
  expect_true(is.na(effects(pooled, "saliency")$partial_ss))
  expect_output(print(effects(pooled, "saliency")), "Objective fraction")
})

test_that("nonconverged partial-pooling fits are refused", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 30, p = 6, seed = 55)
  expect_error(
    design_mfa(
      sim$Y,
      sim$design,
      ~ saliency + probe,
      subject_model = "ridge_slopes",
      subject_penalty = 0.1,
      control = design_mfa_control(multilevel_max_iter = 1, multilevel_tol = 1e-15)
    ),
    "final relative objective change.*tolerance"
  )
})

test_that("the default iteration budget admits convergent weak-pooling fits", {
  sim <- .sim_design_mfa(
    n_subjects = 4,
    trials = c(25, 31, 37, 43),
    p = 6,
    noise_sd = 0.4,
    seed = 10
  )
  args <- list(
    Y = sim$Y,
    design = sim$design,
    formula = ~ saliency + probe,
    subject_model = "ridge_slopes",
    subject_penalty = 0.01
  )
  expect_error(
    do.call(
      design_mfa,
      c(args, list(control = design_mfa_control(multilevel_max_iter = 50)))
    ),
    "did not converge after 50 iterations"
  )
  fit <- do.call(design_mfa, args)

  expect_true(fit$full_model$converged)
  expect_gt(fit$full_model$iterations, 50L)
  expect_lte(fit$full_model$iterations, design_mfa_control()$multilevel_max_iter)
})

test_that("large subject penalty approaches the fixed-intercept group fit", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 35, p = 7, seed = 53)
  args <- list(
    formula = ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    rank = c(saliency = 2, probe = 2)
  )
  fixed <- do.call(design_mfa, c(list(Y = sim$Y, design = sim$design), args))
  pooled <- do.call(design_mfa, c(list(
    Y = sim$Y,
    design = sim$design,
    subject_model = "ridge_slopes",
    subject_penalty = 1e8
  ), args))

  expect_equal(predict(pooled), predict(fixed), tolerance = 1e-7)
  expect_lt(max(vapply(subject_effects(pooled), function(U) max(abs(U)), numeric(1))), 1e-7)
})

test_that("subject storage order preserves the partial-pooling group solution", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = c(30, 34, 38, 42), p = 7, seed = 54)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    contrasts = list(saliency = "poly"),
    subject_model = "ridge_slopes",
    subject_penalty = 0.2,
    rank = 2
  )
  order <- c(4, 2, 1, 3)
  reordered <- design_mfa(
    sim$Y[order],
    sim$design[order],
    ~ saliency,
    contrasts = list(saliency = "poly"),
    subject_model = "ridge_slopes",
    subject_penalty = 0.2,
    rank = 2
  )

  expect_equal(predict(fit), predict(reordered), tolerance = 1e-7)
  expect_equal(
    .dmfa_projector(effect_loadings(fit, "saliency")),
    .dmfa_projector(effect_loadings(reordered, "saliency")),
    tolerance = 1e-7
  )
})
