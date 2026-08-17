test_that("subject and trial storage order do not change fitted effects", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = c(29, 35, 41, 32), p = 9, seed = 30)
  args <- list(
    formula = ~ saliency + probe,
    nuisance = ~ run,
    contrasts = list(saliency = "poly", probe = "sum", run = "sum")
  )
  fit <- do.call(design_mfa, c(list(Y = sim$Y, design = sim$design), args))

  block_order <- c(3, 1, 4, 2)
  fit_blocks <- do.call(design_mfa, c(list(
    Y = sim$Y[block_order],
    design = sim$design[block_order]
  ), args))

  perm <- lapply(sim$Y, function(y) sample.int(nrow(y)))
  Y_rows <- Map(function(y, i) y[i, , drop = FALSE], sim$Y, perm)
  design_rows <- Map(function(d, i) d[i, , drop = FALSE], sim$design, perm)
  names(Y_rows) <- names(design_rows) <- names(sim$Y)
  fit_rows <- do.call(design_mfa, c(list(Y = Y_rows, design = design_rows), args))

  for (term in c("saliency", "probe")) {
    expect_equal(effects(fit, term)$partial_ss, effects(fit_blocks, term)$partial_ss, tolerance = 1e-9)
    expect_equal(effects(fit, term)$partial_ss, effects(fit_rows, term)$partial_ss, tolerance = 1e-9)
    expect_equal(
      .dmfa_projector(effect_loadings(fit, term)),
      .dmfa_projector(effect_loadings(fit_blocks, term)),
      tolerance = 1e-8
    )
    expect_equal(
      .dmfa_projector(effect_loadings(fit, term)),
      .dmfa_projector(effect_loadings(fit_rows, term)),
      tolerance = 1e-8
    )
  }
})

test_that("equivalent contrast bases preserve fitted term subspaces", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 55, p = 10, seed = 31)
  fit_sum <- design_mfa(
    sim$Y,
    sim$design,
    ~ probe,
    contrasts = list(probe = "sum")
  )
  fit_helmert <- design_mfa(
    sim$Y,
    sim$design,
    ~ probe,
    contrasts = list(probe = "helmert")
  )

  expect_equal(effects(fit_sum, "probe")$partial_ss, effects(fit_helmert, "probe")$partial_ss, tolerance = 1e-9)
  expect_equal(
    .dmfa_projector(effect_loadings(fit_sum, "probe")),
    .dmfa_projector(effect_loadings(fit_helmert, "probe")),
    tolerance = 1e-8
  )
})

test_that("default partial interaction hypotheses are contrast invariant", {
  sim <- .sim_design_mfa(n_subjects = 5, trials = 70, p = 10, seed = 311)
  fit_sum <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency * probe,
    contrasts = list(saliency = "poly", probe = "sum")
  )
  fit_treatment <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency * probe,
    contrasts = list(saliency = "poly", probe = "treatment")
  )

  expect_identical(fit_sum$term_scope, "hierarchical")
  expect_equal(
    effects(fit_sum, "saliency")$partial_ss,
    effects(fit_treatment, "saliency")$partial_ss,
    tolerance = 2e-8
  )
  expect_equal(
    .dmfa_projector(effect_loadings(fit_sum, "saliency")),
    .dmfa_projector(effect_loadings(fit_treatment, "saliency")),
    tolerance = 2e-7
  )
})

test_that("continuous reference grids vary and training transforms are reused", {
  set.seed(312)
  ids <- paste0("S", 1:4)
  design <- lapply(seq_along(ids), function(s) {
    data.frame(x = seq(-2, 2, length.out = 32) + s / 100)
  })
  Y <- lapply(seq_along(ids), function(s) {
    x <- design[[s]]$x
    out <- cbind(x, x^2, -0.5 * x + 0.2 * x^2) +
      matrix(rnorm(length(x) * 3, sd = 0.01), length(x), 3)
    colnames(out) <- paste0("v", 1:3)
    out
  })
  names(Y) <- names(design) <- ids

  fit <- design_mfa(Y, design, ~ poly(x, 2))
  effect <- effects(fit, "poly(x, 2)")
  stacked <- do.call(rbind, design)
  rebuilt <- muscal:::.dmfa_build_effect_matrix(fit$design, stacked)
  effect_columns <- unlist(fit$design$term_columns, use.names = FALSE)

  expect_equal(nrow(effect$reference_grid), fit$control$max_numeric_grid)
  expect_gt(diff(range(effect$reference_grid$x)), 0)
  expect_gt(effect$rank, 0L)
  expect_equal(
    unname(rebuilt),
    unname(fit$design$X[, effect_columns, drop = FALSE]),
    tolerance = 1e-12
  )
})

test_that("duplicating subject trials preserves equal-subject fits", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 36, p = 8, seed = 32)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum")
  )

  duplicated_Y <- sim$Y
  duplicated_design <- sim$design
  duplicated_Y[[1]] <- rbind(duplicated_Y[[1]], duplicated_Y[[1]])
  duplicated_design[[1]] <- rbind(duplicated_design[[1]], duplicated_design[[1]])
  fit_duplicated <- design_mfa(
    duplicated_Y,
    duplicated_design,
    ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum")
  )

  expect_equal(predict(fit), predict(fit_duplicated), tolerance = 1e-9)
  expect_equal(
    effect_loadings(fit, "saliency"),
    effect_loadings(fit_duplicated, "saliency"),
    tolerance = 1e-8
  )
})

test_that("partial and sequential effects expose their different contracts", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 45, p = 9, seed = 33)
  partial <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    ss_type = "partial"
  )
  sequential <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    ss_type = "sequential"
  )

  expect_identical(effects(partial, "probe")$ss_type, "partial")
  expect_identical(effects(sequential, "probe")$ss_type, "sequential")
  expect_gte(effects(partial, "probe")$partial_ss, 0)
  expect_gte(effects(sequential, "probe")$partial_ss, 0)
  expect_equal(
    estimability(sequential)$prediction_nonestimable_grid_rows,
    estimability(partial)$prediction_nonestimable_grid_rows
  )
  expect_true(all(estimability(sequential)$estimability_scope == "effect_and_prediction"))
  expect_equal(predict(sequential), predict(partial), tolerance = 1e-10)
  expect_error(
    design_mfa(
      sim$Y,
      sim$design,
      ~ saliency + probe,
      ss_type = "sequential",
      term_scope = "hierarchical"
    ),
    "only defined for partial"
  )
})

test_that("near-collinear designs remain finite and report numerical rank", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 40, p = 7, seed = 34)
  set.seed(35)
  for (s in seq_along(sim$design)) {
    x <- rnorm(nrow(sim$design[[s]]))
    sim$design[[s]]$x <- x
    sim$design[[s]]$x_near <- x + rnorm(length(x), sd = 1e-12)
  }
  fit <- design_mfa(sim$Y, sim$design, ~ x + x_near)

  expect_true(is.finite(fit$full_model$sse))
  expect_true(fit$full_model$rank < ncol(fit$design$X))
  expect_true(all(is.finite(estimability(fit)$partial_ss)))
})
