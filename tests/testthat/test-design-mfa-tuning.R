test_that("effect reconstruction honors Euclidean and spatial metrics", {
  set.seed(80)
  H <- matrix(rnorm(6 * 12), 6, 12)
  H <- sweep(H, 2, colMeans(H), "-")
  weights <- rep(1 / 6, 6)

  euclidean <- muscal:::.dmfa_thin_svd(H, weights, ncomp = 5, rank_tol = 1e-12)
  effect <- list(scores = euclidean$scores, loadings = euclidean$v)
  expect_equal(.dmfa_reconstruct_effect(effect, NULL), H, tolerance = 1e-9)

  L <- Matrix::bandSparse(
    12,
    k = c(-1, 1),
    diagonals = list(rep(1, 11), rep(1, 11))
  )
  L <- Matrix::Diagonal(x = Matrix::rowSums(L)) - L
  metric <- Matrix::Diagonal(12) + 2 * L
  spatial <- muscal:::.dmfa_thin_svd(
    H,
    weights,
    ncomp = 5,
    rank_tol = 1e-12,
    feature_metric = metric
  )
  spatial_effect <- list(scores = spatial$scores, loadings = spatial$v)
  expect_equal(.dmfa_reconstruct_effect(spatial_effect, metric), H, tolerance = 1e-8)
})

test_that("held-out tuning is reproducible and has no subject leakage", {
  sim <- .sim_design_mfa(
    n_subjects = 8,
    trials = c(46, 50, 54, 58, 48, 52, 56, 60),
    p = 12,
    noise_sd = 0.2,
    seed = 81
  )
  grid <- data.frame(
    rank = c(1L, 2L, 3L),
    spatial_penalty = 0,
    subject_penalty = 1
  )

  set.seed(8100)
  caller_seed <- .Random.seed
  tuned1 <- tune_design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    terms = "saliency",
    contrasts = list(saliency = "poly", probe = "sum"),
    grid = grid,
    n_folds = 4,
    seed = 100
  )
  expect_identical(.Random.seed, caller_seed)
  tuned2 <- tune_design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    terms = "saliency",
    contrasts = list(saliency = "poly", probe = "sum"),
    grid = grid,
    n_folds = 4,
    seed = 100
  )

  expect_s3_class(tuned1, "design_mfa_tuning")
  expect_equal(tuned1$results, tuned2$results, tolerance = 1e-12)
  expect_identical(tuned1$best_params, tuned2$best_params)
  expect_equal(tuned1$best_params$rank, 1L)
  expect_s3_class(tuned1$fit, "design_mfa")
  expect_true(all(tuned1$results$n_success == 4L))
  for (fold in tuned1$folds) {
    expect_length(intersect(fold$analysis, fold$assessment), 0L)
    expect_setequal(c(fold$analysis, fold$assessment), sim$subjects)
  }
  for (fold_id in seq_along(tuned1$folds)) {
    expect_identical(
      tuned1$fit_subject_audit[[fold_id]]$analysis,
      sort(tuned1$folds[[fold_id]]$analysis)
    )
    expect_identical(
      tuned1$fit_subject_audit[[fold_id]]$assessment,
      sort(tuned1$folds[[fold_id]]$assessment)
    )
  }
  expect_output(print(tuned1), "Held-out-subject tuning")
  expect_error(
    muscal:::.dmfa_assert_tuning_subjects(c("S1", "S2"), c("S1", "S3"), "analysis", 1L),
    "Subject leakage detected"
  )
})

test_that("subject storage order does not change named-fold tuning", {
  sim <- .sim_design_mfa(n_subjects = 6, trials = 42, p = 9, noise_sd = 0.18, seed = 82)
  grid <- data.frame(
    rank = c(1L, 2L),
    spatial_penalty = 0,
    subject_penalty = 1
  )
  args <- list(
    formula = ~ saliency,
    terms = "saliency",
    grid = grid,
    n_folds = 3,
    seed = 101
  )
  original <- do.call(tune_design_mfa, c(list(Y = sim$Y, design = sim$design), args))
  order <- c(5, 2, 6, 1, 4, 3)
  reordered <- do.call(tune_design_mfa, c(list(
    Y = sim$Y[order],
    design = sim$design[order]
  ), args))

  expect_equal(original$results$score, reordered$results$score, tolerance = 1e-10)
  expect_identical(original$best_params, reordered$best_params)
  expect_identical(original$folds, reordered$folds)
})

test_that("assessment responses cannot alter their fold's analysis fit", {
  sim <- .sim_design_mfa(n_subjects = 6, trials = 30, p = 6, noise_sd = 0.2, seed = 8201)
  grid <- data.frame(
    rank = c(1L, 2L),
    spatial_penalty = 0,
    subject_penalty = 1
  )
  args <- list(
    formula = ~ saliency,
    terms = "saliency",
    grid = grid,
    n_folds = 3,
    seed = 9301,
    refit = FALSE
  )
  baseline <- do.call(tune_design_mfa, c(list(Y = sim$Y, design = sim$design), args))
  assessment <- baseline$folds[[1L]]$assessment
  altered_Y <- sim$Y
  for (id in assessment) {
    signal <- as.numeric(sim$design[[id]]$saliency) - 3
    altered_Y[[id]] <- altered_Y[[id]] + 1000 * outer(signal, seq_len(6))
  }
  altered <- do.call(tune_design_mfa, c(list(Y = altered_Y, design = sim$design), args))

  fold_one <- baseline$fold_scores$fold == 1L
  fingerprint_columns <- c(
    "analysis_rank",
    "analysis_coefficient_sum",
    "analysis_coefficient_sumsq",
    "analysis_coefficient_indexed",
    "analysis_design_sum",
    "analysis_design_sumsq",
    "analysis_design_indexed"
  )
  expect_equal(
    baseline$fold_scores[fold_one, fingerprint_columns],
    altered$fold_scores[fold_one, fingerprint_columns],
    tolerance = 1e-12
  )
  expect_identical(
    baseline$fit_subject_audit[[1L]]$reference_grid,
    altered$fit_subject_audit[[1L]]$reference_grid
  )
  expect_equal(
    baseline$fit_subject_audit[[1L]]$reference_weights,
    altered$fit_subject_audit[[1L]]$reference_weights,
    tolerance = 1e-15
  )
  expect_false(isTRUE(all.equal(
    baseline$fit_subject_audit[[1L]]$target_coefficient_fingerprint,
    altered$fit_subject_audit[[1L]]$target_coefficient_fingerprint,
    tolerance = 1e-12
  )))
  expect_false(isTRUE(all.equal(
    baseline$fold_scores$score[fold_one],
    altered$fold_scores$score[fold_one],
    tolerance = 1e-12
  )))
})

test_that("tuning aligns named trial weights before fold subsetting", {
  sim <- .sim_design_mfa(
    n_subjects = 4,
    trials = c(28, 34, 40, 46),
    p = 7,
    noise_sd = 0.15,
    seed = 821
  )
  weights <- lapply(sim$Y, function(y) seq_len(nrow(y)))
  weights <- weights[rev(names(weights))]
  subject_ids <- paste0("participant_", seq_along(sim$Y))
  block_weight <- stats::setNames(c(1, 2, 3, 4), rev(names(sim$Y)))
  tuned <- tune_design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    subject = subject_ids,
    terms = "saliency",
    trial_weights = weights,
    block_weight = block_weight,
    grid = data.frame(rank = 1L, spatial_penalty = 0, subject_penalty = 1),
    n_folds = 2,
    seed = 1011,
    refit = FALSE
  )

  expect_true(all(tuned$fold_scores$status == "ok"))
  expect_true(is.finite(tuned$best_score))
})

test_that("tuning rejects malformed subject-keyed arguments before folds", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 24, p = 5, seed = 8211)
  grid <- data.frame(rank = 1L, spatial_penalty = 0, subject_penalty = 1)
  weights <- lapply(sim$Y, function(y) rep(1, nrow(y)))
  maps <- lapply(sim$Y, function(y) {
    P <- Matrix::Diagonal(ncol(y))
    dimnames(P) <- list(colnames(y), colnames(y))
    P
  })

  names(weights) <- c(names(sim$Y)[-1L], "extra")
  expect_error(
    tune_design_mfa(
      sim$Y, sim$design, ~ saliency, trial_weights = weights,
      grid = grid, n_folds = 2, refit = FALSE
    ),
    "missing \\{S1\\}; extra \\{extra\\}"
  )

  names(maps) <- c("", names(sim$Y)[-1L])
  expect_error(
    tune_design_mfa(
      sim$Y, sim$design, ~ saliency, alignment = "soft", feature_map = maps,
      grid = transform(grid, alignment_penalty = 0), n_folds = 2, refit = FALSE
    ),
    "missing or empty keys"
  )

  duplicated <- stats::setNames(rep(1, 4), c("S1", "S1", "S3", "S4"))
  expect_error(
    tune_design_mfa(
      sim$Y, sim$design, ~ saliency, block_weight = duplicated,
      grid = grid, n_folds = 2, refit = FALSE
    ),
    "duplicated keys \\{S1\\}"
  )

  expect_error(
    tune_design_mfa(
      sim$Y, sim$design, ~ saliency, trial_weights = weights[1:3],
      grid = grid, n_folds = 2, refit = FALSE
    ),
    "missing"
  )
})

test_that("fold subsetting keeps every subject-keyed argument on Y block keys", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = c(20, 22, 24, 26), p = 5, seed = 822)
  weights <- lapply(sim$Y, function(y) seq_len(nrow(y)))
  maps <- lapply(seq_along(sim$Y), function(s) matrix(s, 2, 2))
  names(weights) <- names(maps) <- rev(names(sim$Y))
  block_weight <- stats::setNames(seq_along(sim$Y), rev(names(sim$Y)))
  indices <- c(2L, 4L)
  canonical <- muscal:::.dmfa_canonicalize_tuning_args(
    list(trial_weights = weights, feature_map = maps, block_weight = block_weight),
    sim$Y
  )
  out <- muscal:::.dmfa_subset_tuning_args(
    canonical,
    indices,
    sim$Y,
    sim$design,
    paste0("participant_", seq_along(sim$Y)),
    length(sim$Y)
  )

  expect_identical(names(out$trial_weights), names(sim$Y)[indices])
  expect_identical(names(out$feature_map), names(sim$Y)[indices])
  expect_identical(names(out$block_weight), names(sim$Y)[indices])
  expect_identical(out$trial_weights[[1L]], weights[[names(sim$Y)[indices[[1L]]]]])
  expect_identical(out$feature_map[[2L]], maps[[names(sim$Y)[indices[[2L]]]]])
  expect_equal(
    unname(out$block_weight),
    unname(block_weight[names(sim$Y)][indices])
  )
})

test_that("tuning evaluates spatial and partial-pooling candidates", {
  sim <- .sim_design_mfa(n_subjects = 6, trials = 45, p = 10, noise_sd = 0.15, seed = 83)
  adjacency <- Matrix::bandSparse(
    10,
    k = c(-1, 1),
    diagonals = list(rep(1, 9), rep(1, 9))
  )
  laplacian <- Matrix::Diagonal(x = Matrix::rowSums(adjacency)) - adjacency
  dimnames(laplacian) <- list(sim$feature_names, sim$feature_names)
  grid <- expand.grid(
    rank = c(1L, 2L),
    spatial_penalty = c(0, 2),
    subject_penalty = c(0.1, 2),
    KEEP.OUT.ATTRS = FALSE
  )

  tuned <- tune_design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    terms = "saliency",
    subject_model = "ridge_slopes",
    random_terms = "saliency",
    spatial_laplacian = laplacian,
    grid = grid,
    n_folds = 3,
    seed = 102,
    control = design_mfa_control(multilevel_max_iter = 100)
  )

  expect_equal(nrow(tuned$results), 8L)
  expect_true(all(is.finite(tuned$results$score)))
  expect_true(all(tuned$results$n_success == 3L))
  expect_true(tuned$best_params$rank %in% c(1L, 2L))
  expect_true(tuned$best_params$spatial_penalty %in% c(0, 2))
  expect_true(tuned$best_params$subject_penalty %in% c(0.1, 2))
})

test_that("held-out tuning subsets coordinate maps with their subjects", {
  sim <- .sim_design_mfa_misregistration(
    n_subjects = 6,
    trials = 34,
    p = 15,
    noise_sd = 0.25,
    seed = 831
  )
  tuned <- tune_design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    terms = "saliency",
    alignment = "soft",
    feature_map = sim$maps,
    contrasts = list(saliency = "poly"),
    grid = data.frame(
      rank = rep(c(1L, 2L), each = 2L),
      spatial_penalty = 0,
      subject_penalty = 1,
      alignment_penalty = rep(c(0.02, 0.1), 2L)
    ),
    n_folds = 3,
    seed = 104,
    control = design_mfa_control(soft_max_iter = 800, soft_tol = 1e-8)
  )

  expect_true(all(tuned$results$n_success == 3L))
  expect_true(all(is.finite(tuned$results$score)))
  expect_identical(tuned$fit$alignment, "soft")
  expect_identical(names(tuned$fit$feature_map), names(sim$Y))
  expect_true(tuned$best_params$alignment_penalty %in% c(0.02, 0.1))
  expect_equal(tuned$fit$alignment_penalty, tuned$best_params$alignment_penalty)
})

test_that("tuning rejects invalid folds, grids, and unsupported overrides", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 30, p = 7, seed = 84)
  expect_error(
    tune_design_mfa(sim$Y, sim$design, ~ saliency),
    "at least four subjects"
  )

  sim <- .sim_design_mfa(n_subjects = 4, trials = 30, p = 7, seed = 85)
  expect_error(
    tune_design_mfa(
      sim$Y,
      sim$design,
      ~ saliency,
      grid = data.frame(rank = 0, spatial_penalty = 0, subject_penalty = 1)
    ),
    "positive integer"
  )
  expect_error(
    tune_design_mfa(
      sim$Y,
      sim$design,
      ~ saliency,
      grid = data.frame(rank = 1, spatial_penalty = 1, subject_penalty = 1)
    ),
    "requires `spatial_laplacian`"
  )
  expect_error(
    tune_design_mfa(sim$Y, sim$design, ~ saliency, rank = 2),
    "controlled by `grid`"
  )
  expect_error(tune_design_mfa(sim$Y, sim$design, ~ saliency, n_folds = 2.5), "integer")
  expect_error(tune_design_mfa(sim$Y, sim$design, ~ saliency, refit = 1), "TRUE or FALSE")
  expect_error(
    tune_design_mfa(
      sim$Y,
      sim$design,
      ~ saliency,
      grid = data.frame(rank = 1, alignment_penalty = 0.1)
    ),
    "requires `alignment = \\\"soft\\\"`"
  )
})

test_that("tuning errors retain candidate and fold diagnostics", {
  sim <- .sim_design_mfa(
    n_subjects = 4,
    trials = c(25, 31, 37, 43),
    p = 6,
    noise_sd = 0.4,
    seed = 10
  )
  expect_error(
    tune_design_mfa(
      sim$Y,
      sim$design,
      ~ saliency + probe,
      subject_model = "ridge_slopes",
      grid = data.frame(rank = 1, spatial_penalty = 0, subject_penalty = 0.01),
      n_folds = 2,
      refit = FALSE,
      control = design_mfa_control(multilevel_max_iter = 1, multilevel_tol = 1e-15)
    ),
    "candidate 1, fold 1: The model could not be fitted: full partial-pooling solver did not converge"
  )
})

test_that("successful tuning warns about excluded candidates and retains errors", {
  sim <- .sim_design_mfa(n_subjects = 6, trials = 45, p = 6, noise_sd = 0.4, seed = 10)
  expect_warning(
    tuned <- tune_design_mfa(
      sim$Y,
      sim$design,
      ~ saliency + probe,
      subject_model = "ridge_slopes",
      grid = data.frame(
        rank = 1,
        spatial_penalty = 0,
        subject_penalty = c(0.001, 100)
      ),
      n_folds = 3,
      refit = FALSE,
      control = design_mfa_control(multilevel_max_iter = 10)
    ),
    "3 of 6 candidate-fold fits failed"
  )

  expect_identical(tuned$results$n_success, c(0L, 3L))
  expect_equal(sum(tuned$failure_summary$n_failures), 3L)
  expect_true(all(grepl("did not converge", tuned$failure_summary$message)))
  expect_identical(tuned$failure_summary$stage, "analysis_fit")
  expect_identical(tuned$failure_summary$reason, "nonconvergence")
  expect_identical(tuned$failure_summary$condition_class, "simpleError")

  messages <- character(0)
  expect_warning(
    withCallingHandlers(
      tune_design_mfa(
        sim$Y,
        sim$design,
        ~ saliency + probe,
        subject_model = "ridge_slopes",
        grid = data.frame(
          rank = 1,
          spatial_penalty = 0,
          subject_penalty = c(0.001, 100)
        ),
        n_folds = 3,
        refit = FALSE,
        verbose = TRUE,
        control = design_mfa_control(multilevel_max_iter = 10)
      ),
      message = function(condition) {
        messages <<- c(messages, conditionMessage(condition))
        invokeRestart("muffleMessage")
      }
    ),
    "3 of 6 candidate-fold fits failed"
  )
  expect_true(any(grepl("candidate 1 fold", messages, fixed = TRUE)))
})

test_that("nested tuning gives reproducible outer-subject evaluation", {
  sim <- .sim_design_mfa(
    n_subjects = 6,
    trials = c(28, 30, 32, 34, 36, 38),
    p = 8,
    noise_sd = 0.2,
    seed = 86
  )
  grid <- data.frame(
    rank = c(1L, 2L),
    spatial_penalty = 0,
    subject_penalty = 1
  )
  args <- list(
    formula = ~ saliency,
    terms = "saliency",
    grid = grid,
    outer_folds = 3,
    inner_folds = 2,
    seed = 103
  )

  nested1 <- do.call(nested_tune_design_mfa, c(list(Y = sim$Y, design = sim$design), args))
  nested2 <- do.call(nested_tune_design_mfa, c(list(Y = sim$Y, design = sim$design), args))

  expect_s3_class(nested1, "nested_design_mfa_tuning")
  expect_equal(nested1$outer_results, nested2$outer_results, tolerance = 1e-12)
  expect_equal(nested1$mean_score, mean(nested1$outer_results$score), tolerance = 1e-14)
  expect_s3_class(nested1$final_tuning, "design_mfa_tuning")
  expect_s3_class(nested1$fit, "design_mfa")
  expect_true(all(is.finite(nested1$outer_results$score)))
  for (i in seq_along(nested1$outer_folds)) {
    outer <- nested1$outer_folds[[i]]
    inner <- nested1$outer_tuning[[i]]$folds
    for (fold in inner) {
      expect_true(all(fold$analysis %in% outer$analysis))
      expect_true(all(fold$assessment %in% outer$analysis))
      expect_length(intersect(fold$analysis, outer$assessment), 0L)
      expect_length(intersect(fold$assessment, outer$assessment), 0L)
    }
  }
  expect_output(print(nested1), "Nested Design MFA tuning")
})

test_that("outer assessment responses cannot enter nested analysis or selection", {
  sim <- .sim_design_mfa(n_subjects = 6, trials = 26, p = 5, noise_sd = 0.2, seed = 8601)
  grid <- data.frame(
    rank = c(1L, 2L),
    spatial_penalty = 0,
    subject_penalty = 1
  )
  args <- list(
    formula = ~ saliency,
    terms = "saliency",
    grid = grid,
    outer_folds = 3,
    inner_folds = 2,
    seed = 9601,
    refit = FALSE
  )
  baseline <- do.call(nested_tune_design_mfa, c(list(Y = sim$Y, design = sim$design), args))
  assessment <- baseline$outer_folds[[1L]]$assessment
  altered_Y <- sim$Y
  for (id in assessment) {
    signal <- as.numeric(sim$design[[id]]$saliency) - 3
    altered_Y[[id]] <- altered_Y[[id]] + 1500 * outer(signal, seq_len(5))
  }
  altered <- do.call(nested_tune_design_mfa, c(list(Y = altered_Y, design = sim$design), args))

  base_audit <- baseline$outer_fit_subject_audit[[1L]]
  altered_audit <- altered$outer_fit_subject_audit[[1L]]
  expect_identical(base_audit$analysis, altered_audit$analysis)
  expect_equal(
    base_audit$analysis_coefficient_fingerprint,
    altered_audit$analysis_coefficient_fingerprint,
    tolerance = 1e-12
  )
  expect_equal(
    base_audit$analysis_design_fingerprint,
    altered_audit$analysis_design_fingerprint,
    tolerance = 1e-12
  )
  expect_identical(
    baseline$outer_tuning[[1L]]$best_params,
    altered$outer_tuning[[1L]]$best_params
  )
  expect_false(isTRUE(all.equal(
    base_audit$target_coefficient_fingerprint,
    altered_audit$target_coefficient_fingerprint,
    tolerance = 1e-12
  )))
  expect_false(isTRUE(all.equal(
    baseline$outer_results$score[[1L]],
    altered$outer_results$score[[1L]],
    tolerance = 1e-12
  )))
})

test_that("seeded tuning restores present and absent RNG state on success and error", {
  sim <- .sim_design_mfa(n_subjects = 6, trials = 24, p = 5, seed = 861)
  grid <- data.frame(rank = 1L, spatial_penalty = 0, subject_penalty = 1)
  original_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (original_exists) original_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (original_exists) {
      assign(".Random.seed", original_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  check_state <- function(seed_present, callback, errors = FALSE) {
    if (seed_present) {
      set.seed(9861)
      before <- .Random.seed
    } else {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
      before <- NULL
    }
    if (errors) expect_error(callback()) else expect_no_error(callback())
    if (seed_present) {
      expect_identical(.Random.seed, before)
    } else {
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    }
  }

  for (present in c(FALSE, TRUE)) {
    check_state(present, function() {
      out <- tune_design_mfa(
        sim$Y, sim$design, ~ saliency, grid = grid,
        n_folds = 3, seed = 21.0, refit = FALSE
      )
      expect_identical(out$seed, 21L)
    })
    check_state(present, function() {
      tune_design_mfa(
        sim$Y,
        sim$design,
        ~ saliency,
        subject_model = "ridge_slopes",
        grid = transform(grid, subject_penalty = 0.001),
        n_folds = 3,
        seed = 22.0,
        refit = FALSE,
        control = design_mfa_control(multilevel_max_iter = 1, multilevel_tol = 1e-15)
      )
    }, errors = TRUE)
    check_state(present, function() {
      out <- nested_tune_design_mfa(
        sim$Y, sim$design, ~ saliency, grid = grid,
        outer_folds = 3, inner_folds = 2, seed = 23.0, refit = FALSE
      )
      expect_identical(out$seed, 23L)
    })
    bad_weights <- stats::setNames(rep(1, 6), c("bad", names(sim$Y)[-1L]))
    check_state(present, function() {
      nested_tune_design_mfa(
        sim$Y, sim$design, ~ saliency, grid = grid,
        block_weight = bad_weights,
        outer_folds = 3, inner_folds = 2, seed = 24.0, refit = FALSE
      )
    }, errors = TRUE)
  }
})
