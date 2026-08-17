.dmfa_reference_effect_truth <- function(effect, map, variable = "saliency") {
  grid_value <- effect$reference_grid[[variable]]
  score <- if (is.factor(grid_value)) as.numeric(grid_value) - mean(seq_along(levels(grid_value))) else grid_value
  H <- outer(score, map)
  sweep(H, 2L, colSums(H * effect$reference_weights), "-")
}

.dmfa_relative_effect_error <- function(fit, term, truth) {
  estimate <- effects(fit, term)$effect_matrix
  weights <- effects(fit, term)$reference_weights
  sum(rowSums((estimate - truth)^2) * weights) /
    sum(rowSums(truth^2) * weights)
}

.dmfa_effect_energy <- function(fit, term) {
  effect <- effects(fit, term)
  sum(rowSums(effect$effect_matrix^2) * effect$reference_weights)
}

test_that("fair same-index spatial baseline retains known correspondence", {
  sim <- .sim_design_mfa_misregistration(seed = 80)
  common <- list(
    design = sim$design,
    formula = ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    rank = 1,
    store = "fitted"
  )
  soft <- do.call(design_mfa, c(
    list(Y = sim$Y),
    common,
    list(
      alignment = "soft",
      feature_map = sim$maps,
      alignment_penalty = 0.1,
      control = design_mfa_control(soft_max_iter = 1000, soft_tol = 1e-9)
    )
  ))
  naive <- do.call(design_mfa, c(list(Y = sim$Y_naive), common))
  truth <- .dmfa_reference_effect_truth(
    effects(soft, "saliency"),
    sim$truth$saliency
  )
  soft_error <- .dmfa_relative_effect_error(soft, "saliency", truth)
  naive_error <- .dmfa_relative_effect_error(naive, "saliency", truth)

  expect_lt(soft_error, 0.20)
  expect_lt(naive_error, 0.03)
  expect_lt(naive_error, soft_error)
  expect_gt(abs(cor(effect_map(soft, "saliency"), sim$truth$saliency)), 0.95)
  expect_identical(soft$computation, "matrix_free_soft_alignment")
})

test_that("spatial null-energy diagnostics remain finite but are not inferential calibration", {
  seeds <- 81:84
  energy <- vapply(seeds, function(seed) {
    sim <- .sim_design_mfa_misregistration(
      n_subjects = 6,
      trials = 36,
      p = 20,
      effect_scale = 0,
      seed = seed
    )
    common <- list(
      design = sim$design,
      formula = ~ saliency,
      contrasts = list(saliency = "poly"),
      rank = 1,
      store = "fitted"
    )
    soft <- do.call(design_mfa, c(
      list(Y = sim$Y),
      common,
      list(
        alignment = "soft",
        feature_map = sim$maps,
        alignment_penalty = 0.1,
        control = design_mfa_control(soft_max_iter = 800, soft_tol = 1e-8)
      )
    ))
    naive <- do.call(design_mfa, c(list(Y = sim$Y_naive), common))
    c(soft = .dmfa_effect_energy(soft, "saliency"), naive = .dmfa_effect_energy(naive, "saliency"))
  }, numeric(2))

  expect_true(all(is.finite(energy)))
  expect_true(all(energy >= 0))
})
