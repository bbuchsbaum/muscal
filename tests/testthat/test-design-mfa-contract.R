test_that("design_mfa validates shared-grid and design contracts", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 24, p = 7, seed = 10)

  unnamed <- lapply(sim$Y, unname)
  expect_error(
    design_mfa(unnamed, sim$design, ~ saliency),
    "feature identity is ambiguous"
  )

  bad_schema <- sim$Y
  colnames(bad_schema[[2]])[1] <- "wrong_feature"
  expect_error(
    design_mfa(bad_schema, sim$design, ~ saliency),
    "same named feature schema"
  )

  bad_design <- sim$design
  bad_design[[1]]$probe <- as.character(bad_design[[1]]$probe)
  expect_error(
    design_mfa(sim$Y, bad_design, ~ probe),
    "Convert it to a factor"
  )

  bad_levels <- sim$design
  bad_levels[[2]]$probe <- factor(bad_levels[[2]]$probe, levels = rev(levels(bad_levels[[2]]$probe)))
  expect_error(
    design_mfa(sim$Y, bad_levels, ~ probe),
    "identical levels in identical order"
  )

  nonfinite <- sim$Y
  nonfinite[[1]][1, 1] <- Inf
  expect_error(design_mfa(nonfinite, sim$design, ~ saliency), "non-finite")

  expect_error(
    design_mfa(sim$Y, sim$design, ~ saliency, alignment = "soft"),
    "requires one sparse `feature_map`"
  )
  expect_error(
    design_mfa(sim$Y, sim$design, ~ saliency, subject_model = "ridge_slopes", subject_penalty = 0),
    "positive finite"
  )
  expect_error(
    design_mfa(sim$Y, sim$design, ~ saliency + run, nuisance = ~ run),
    "both `formula` and `nuisance`"
  )
})

test_that("equal-subject weighting is explicit and insensitive to trial counts", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = c(18, 27, 36, 45), p = 8, seed = 11)
  precision <- lapply(sim$Y, function(y) seq_len(nrow(y)))
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    nuisance = ~ run,
    contrasts = list(saliency = "poly", probe = "sum", run = "sum"),
    trial_weights = precision,
    block_weight = "equal_subject"
  )

  expect_equal(unname(vapply(fit$weights, sum, numeric(1))), rep(0.25, 4), tolerance = 1e-14)
  expect_equal(sum(unlist(fit$weights)), 1, tolerance = 1e-14)
  expect_equal(fit$fit_spec$task, "design_decomposition")
  expect_equal(fit$fit_spec$oos_types, "response")
  expect_true(fit$fit_spec$refit_supported)
  expect_identical(fit$fit_spec$refit$resample_unit, "subjects")
  expect_output(print(fit), "Alignment: shared")
})

test_that("design_mfa records nonestimable missing-cell interactions", {
  ids <- paste0("S", 1:3)
  design <- lapply(ids, function(id) {
    grid <- expand.grid(a = c("a1", "a2"), b = c("b1", "b2"), rep = 1:4)
    grid <- grid[!(grid$a == "a2" & grid$b == "b2"), ]
    data.frame(
      a = factor(grid$a, levels = c("a1", "a2")),
      b = factor(grid$b, levels = c("b1", "b2"))
    )
  })
  Y <- lapply(design, function(d) {
    out <- matrix(rnorm(nrow(d) * 5), nrow(d), 5)
    colnames(out) <- paste0("v", 1:5)
    out
  })
  names(Y) <- names(design) <- ids

  fit <- design_mfa(Y, design, ~ a * b, contrasts = list(a = "sum", b = "sum"))
  interaction <- effects(fit, "a:b")

  expect_false(interaction$estimable)
  expect_true(interaction$status %in% c("nonestimable", "nonestimable_reference"))
  expect_true(interaction$estimable_df == 0L || any(!interaction$estimable_grid))
  expect_match(paste(capture.output(print(interaction)), collapse = "\n"), "Status")
})

test_that("between-subject aliases do not suppress within-subject effects", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 42, p = 8, seed = 131)
  ages <- c(22, 31, 47, 58)
  for (s in seq_along(sim$design)) {
    sim$design[[s]]$age <- rep(ages[[s]], nrow(sim$design[[s]]))
  }
  base <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    contrasts = list(saliency = "poly")
  )
  adjusted <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + age,
    contrasts = list(saliency = "poly"),
    terms = "saliency"
  )

  effect <- effects(adjusted, "saliency")
  expect_true(effect$estimable)
  expect_true(any(!effect$raw_estimable_grid))
  table <- estimability(adjusted)
  expect_equal(table$effect_nonestimable_grid_rows, 0L)
  expect_gt(table$prediction_nonestimable_grid_rows, 0L)
  expect_identical(table$estimability_scope, "effect_only")
  expect_error(predict(adjusted), "full-response prediction")
  expect_equal(effect$partial_ss, effects(base, "saliency")$partial_ss, tolerance = 1e-8)
  expect_equal(
    .dmfa_projector(effect_loadings(adjusted, "saliency")),
    .dmfa_projector(effect_loadings(base, "saliency")),
    tolerance = 1e-7
  )
  expect_length(effect_map(adjusted, "saliency", 1L), adjusted$n_features)
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(adjusted, "saliency"))
})

test_that("aliased nuisance scaling leaves the exact term functional invariant", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 42, p = 8, seed = 1311)
  ages <- c(22, 31, 47, 58)
  for (s in seq_along(sim$design)) {
    sim$design[[s]]$age <- rep(ages[[s]], nrow(sim$design[[s]]))
  }
  fit_one <- function(multiplier) {
    design <- sim$design
    for (s in seq_along(design)) design[[s]]$age <- design[[s]]$age * multiplier
    design_mfa(
      sim$Y,
      design,
      ~ saliency + age,
      contrasts = list(saliency = "poly"),
      terms = "saliency",
      store = "fitted"
    )
  }
  fits <- lapply(c(0.01, 1, 1000), fit_one)
  base <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    contrasts = list(saliency = "poly"),
    store = "fitted"
  )
  canonical_rows <- function(effect) !duplicated(effect$reference_grid$saliency)
  base_effect <- effects(base, "saliency")

  for (fit in fits) {
    effect <- effects(fit, "saliency")
    rows <- canonical_rows(effect)
    expect_true(effect$estimable)
    expect_identical(effect$status, "estimable")
    expect_equal(effect$partial_ss, base_effect$partial_ss, tolerance = 1e-10)
    expect_equal(
      effect$effect_matrix[rows, , drop = FALSE],
      base_effect$effect_matrix,
      tolerance = 1e-10
    )
    for (indices in split(seq_len(nrow(effect$reference_grid)), effect$reference_grid$saliency)) {
      expect_equal(
        effect$effect_matrix[indices, , drop = FALSE],
        effect$effect_matrix[rep(indices[[1L]], length(indices)), , drop = FALSE],
        tolerance = 1e-12
      )
    }
    expect_gt(sum(!effect$prediction_estimable_grid), 0L)
    expect_error(predict(fit), "full-response prediction")
  }

  reference <- effects(fits[[2L]], "saliency")
  for (fit in fits[c(1L, 3L)]) {
    effect <- effects(fit, "saliency")
    expect_equal(effect$singular_values, reference$singular_values, tolerance = 1e-10)
    expect_equal(effect$scores, reference$scores, tolerance = 1e-10)
    expect_equal(
      .dmfa_projector(effect$loadings),
      .dmfa_projector(reference$loadings),
      tolerance = 1e-10
    )
    expect_identical(effect$rank, reference$rank)
  }
})

test_that("WLS sufficient-statistic solver never derives SSE by subtraction", {
  expect_error(
    muscal:::.dmfa_wls_from_stats(
      XtX = diag(2),
      XtY = matrix(1, 2, 1),
      response_ss = 1e20,
      coefficient_names = c("a", "b"),
      feature_names = "v1",
      rank_tol = 1e-8,
      sse = NULL
    ),
    "must be supplied directly"
  )
})

test_that("default plot selects the first retained estimable effect", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 25, p = 6, seed = 1312)
  fit <- design_mfa(sim$Y, sim$design, ~ saliency, rank = 1)
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  plotted <- plot(fit)
  expect_equal(nrow(plotted), nrow(effects(fit, "saliency")$reference_grid))

  one_row <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    reference_grid = data.frame(
      saliency = ordered(20, levels = c(20, 40, 60, 80, 100))
    )
  )
  expect_identical(effects(one_row, "saliency")$status, "degenerate")
  expect_error(plot(one_row), "No estimable fitted effect")

  no_components <- fit
  no_components$effects[["saliency"]]$rank <- 0L
  no_components$effects[["saliency"]]$scores <- NULL
  expect_error(plot(no_components), "No estimable fitted effect")
})

test_that("control and edge arguments fail at the API boundary", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 20, p = 6, seed = 12)
  expect_error(design_mfa_control(rank_tol = 0), "positive finite")
  expect_error(design_mfa_control(max_grid_rows = 0), "positive integer")
  expect_error(design_mfa_control(feature_chunk_size = 0), "positive integer")
  expect_error(
    design_mfa(sim$Y, sim$design, ~ saliency, trial_weights = list(1, 1)),
    "one vector per subject"
  )
  expect_error(
    design_mfa(sim$Y, sim$design, ~ saliency, terms = "not_a_term"),
    "Unknown requested term"
  )
  expect_error(
    design_mfa(sim$Y, sim$design, ~ saliency, reference_weight = -1),
    "one nonnegative value"
  )
})

test_that("multidesign inputs adapt to the same list engine", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 25, p = 7, seed = 13)
  stacked_x <- do.call(rbind, sim$Y)
  stacked_design <- do.call(rbind, Map(function(d, id) {
    d$subject_id <- id
    d
  }, sim$design, sim$subjects))
  md <- multidesign::multidesign(stacked_x, stacked_design)

  fit_md <- design_mfa(
    md,
    formula = ~ saliency + probe,
    subject = subject_id,
    contrasts = list(saliency = "poly", probe = "sum")
  )
  fit_list <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum")
  )

  expect_s3_class(fit_md, "design_mfa")
  expect_equal(predict(fit_md), predict(fit_list), tolerance = 1e-9)
})

test_that("generated documentation pins estimability and method terminology", {
  package_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  testthat::skip_if_not(
    file.exists(file.path(package_root, "man", "predict.design_mfa.Rd")),
    "source Rd files are unavailable in the installed-package test copy"
  )
  predict_rd <- paste(
    readLines(file.path(package_root, "man", "predict.design_mfa.Rd"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(predict_rd, "Full-response prediction is stricter")
  expect_match(predict_rd, "prediction_nonestimable_grid_rows", fixed = TRUE)

  public_docs <- vapply(
    c("README.md", "vignettes/design_mfa.Rmd", "man/design_mfa.Rd"),
    function(path) gsub(
      "[[:space:]]+", " ",
      paste(readLines(file.path(package_root, path), warn = FALSE), collapse = "\n")
    ),
    character(1)
  )
  expect_true(all(grepl("ASCA\\+/partial-RDA", public_docs)))
  expect_true(all(grepl("block normalization", public_docs, fixed = TRUE)))

  infer_rd <- paste(
    readLines(file.path(package_root, "man", "infer_design_mfa.Rd"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(infer_rd, "must not be read as a", fixed = TRUE)
  expect_match(infer_rd, "n_success", fixed = TRUE)
})

test_that("sequential prediction estimability uses the global model", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 42, p = 8, seed = 131)
  ages <- c(22, 31, 47, 58)
  for (s in seq_along(sim$design)) {
    sim$design[[s]]$age <- rep(ages[[s]], nrow(sim$design[[s]]))
  }
  args <- list(
    Y = sim$Y,
    design = sim$design,
    formula = ~ saliency + age,
    contrasts = list(saliency = "poly"),
    terms = "saliency"
  )
  partial <- do.call(design_mfa, c(args, list(ss_type = "partial")))
  sequential <- do.call(design_mfa, c(args, list(ss_type = "sequential")))

  expect_true(effects(sequential, "saliency")$estimable)
  expect_equal(estimability(sequential)$effect_nonestimable_grid_rows, 0L)
  expect_gt(estimability(sequential)$prediction_nonestimable_grid_rows, 0L)
  expect_identical(estimability(sequential)$estimability_scope, "effect_only")
  expect_equal(
    estimability(sequential)$prediction_nonestimable_grid_rows,
    estimability(partial)$prediction_nonestimable_grid_rows
  )
  expect_identical(
    effects(sequential, "saliency")$prediction_estimable_grid,
    effects(partial, "saliency")$prediction_estimable_grid
  )
  global_check <- muscal:::.dmfa_estimable_rows(
    sequential$design$X_reference,
    sequential$full_model,
    sequential$control$estimability_tol
  )
  expect_identical(
    effects(sequential, "saliency")$prediction_estimable_grid,
    global_check$estimable
  )
  expect_error(predict(sequential), "full-response prediction")
  expect_false("nonestimable_grid_rows" %in% names(estimability(sequential)))
})

test_that("stored fit drops dead sufficient statistics and isolates the refit closure", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 24, p = 6, seed = 14)
  fit <- design_mfa(sim$Y, sim$design, ~ saliency + probe)

  expect_false("XtY" %in% names(fit$full_model))
  expect_false("XtX_equilibrated" %in% names(fit$full_model))
  expect_true(all(c("coefficient", "row_projector", "column_scales", "XtX") %in% names(fit$full_model)))

  fit_fn <- fit$fit_spec$refit$fit_fn
  env <- environment(fit_fn)
  env_names <- ls(env, all.names = TRUE)
  expect_false(any(c("Y", "Y_stacked", "fit", "global_full", "effect_list", "compiled") %in% env_names))
  expect_true(all(c("formula", "control", "ss_type") %in% env_names))
  expect_identical(parent.env(env), asNamespace("muscal"))
  expect_true(exists("Y", envir = environment(fit$fit_spec$refit$fit_fn), inherits = FALSE) == FALSE)

  refit <- fit_fn(fit$fit_spec$refit$data)
  expect_s3_class(refit, "design_mfa")
  expect_equal(effects(refit, "saliency")$partial_ss, effects(fit, "saliency")$partial_ss, tolerance = 1e-10)
})

test_that("declared xlev levels appear on the default reference grid", {
  stacked <- data.frame(
    cond = factor(c("A", "B", "A"), levels = c("A", "B")),
    stringsAsFactors = FALSE
  )
  control <- design_mfa_control()
  with_xlev <- muscal:::.dmfa_make_reference_grid(
    stacked,
    ~ cond,
    xlev = list(cond = c("A", "B", "C")),
    reference_grid = NULL,
    control = control
  )
  without_xlev <- muscal:::.dmfa_make_reference_grid(
    stacked,
    ~ cond,
    xlev = NULL,
    reference_grid = NULL,
    control = control
  )
  expect_identical(levels(with_xlev$cond), c("A", "B", "C"))
  expect_equal(nrow(with_xlev), 3L)
  expect_identical(levels(without_xlev$cond), c("A", "B"))
  expect_equal(nrow(without_xlev), 2L)

  sim <- .sim_design_mfa(n_subjects = 3, trials = 20, p = 5, seed = 15)
  for (s in seq_along(sim$design)) {
    sim$design[[s]]$probe <- factor(
      as.character(sim$design[[s]]$probe),
      levels = c("old", "new", "lure", "unseen")
    )
  }
  fit <- design_mfa(sim$Y, sim$design, ~ probe, contrasts = list(probe = "sum"))
  expect_identical(levels(fit$design$reference_grid$probe), c("old", "new", "lure", "unseen"))
  expect_equal(nrow(fit$design$reference_grid), 4L)
})
