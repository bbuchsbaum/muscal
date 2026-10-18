test_that("design_mfa adapter claims effect_significance and not a component ladder", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.dmfa_multifer_effect_contract_available())

  sim <- .sim_design_mfa(n_subjects = 4, trials = 24, p = 5, seed = 91)
  fit <- design_mfa(sim$Y, sim$design, ~ saliency, rank = 1)
  adapter <- as_multifer_adapter(fit)

  expect_s3_class(adapter, "multifer_adapter")
  expect_true(multifer::adapter_supports(
    adapter, "adapter", "design_contrast", "effect_significance"
  ))
  if (muscal:::.dmfa_multifer_stability_contract_available()) {
    expect_true(multifer::adapter_supports(
      adapter, "adapter", "design_contrast", "statistic_stability"
    ))
  }
  expect_false(multifer::adapter_supports(
    adapter, "adapter", "design_contrast", "component_significance"
  ))
  expect_equal(adapter$null_spec$kind, "caller_owned")
  expect_null(adapter$null_action)
  expect_null(adapter$residualize)
  expect_equal(adapter$effects(fit), "saliency")
  expect_equal(
    adapter$effect_stat(fit, fit$fit_spec$refit$data, effect = "saliency"),
    if (is.finite(effects(fit, "saliency")$partial_ss)) {
      effects(fit, "saliency")$partial_ss
    } else {
      effects(fit, "saliency")$objective_delta
    }
  )
})

test_that("serial permutation uses multifer infer and keeps the public result shape", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.dmfa_multifer_effect_contract_available())

  sim <- .sim_design_mfa(n_subjects = 4, trials = 28, p = 6, seed = 92)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    contrasts = list(saliency = "poly"),
    rank = 1
  )
  within_subject_shuffle <- function(data, term, replicate) {
    out <- data
    out$design <- lapply(out$design, function(d) {
      d$saliency <- sample(d$saliency)
      d
    })
    out
  }
  perm <- infer_design_mfa(
    fit,
    "saliency",
    method = "permutation",
    exchangeability = within_subject_shuffle,
    nrep = 5,
    seed = 93
  )

  expect_identical(perm$engine, "multifer")
  expect_identical(perm$validity, "user_supplied_exchangeability")
  expect_identical(perm$statistics, c("partial_ss", "singular_values"))
  expect_equal(perm$n_success, 5)
  expect_true(!is.null(perm$summary$singular_values))
  null <- vapply(perm$replicates, `[[`, numeric(1), "partial_ss")
  expect_equal(
    perm$summary$partial_ss$p_value,
    (1 + sum(null >= perm$observed$partial_ss)) / (length(null) + 1)
  )
})

test_that("unchanged permutation draws retain every requested statistic", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.dmfa_multifer_effect_contract_available())

  sim <- .sim_design_mfa(n_subjects = 4, trials = 24, p = 5, seed = 97)
  fit <- design_mfa(sim$Y, sim$design, ~ saliency, rank = 1)
  for (unchanged_only in c(TRUE, FALSE)) {
    null_action <- function(data, term, replicate) {
      if (unchanged_only || replicate %% 2L == 1L) return(data)
      data$design <- lapply(data$design, function(d) {
        d$saliency <- sample(d$saliency)
        d
      })
      data
    }
    run <- function() infer_design_mfa(
      fit, "saliency", method = "permutation", nrep = 4, seed = 98,
      exchangeability = null_action
    )
    perm <- run()
    native <- with_mocked_bindings(
      run(),
      .dmfa_multifer_effect_contract_available = function() FALSE,
      .package = "muscal"
    )
    expect_identical(perm$engine, "multifer")
    expect_identical(native$engine, "native")
    expect_equal(perm$n_success, 4L)
    expect_equal(perm$replicates, native$replicates)
    expect_equal(perm$summary, native$summary)
    if (unchanged_only) {
      expect_equal(perm$summary$partial_ss$p_value, 1)
      expect_true(all(perm$summary$singular_values$p_value == 1))
    }
  }
})

test_that("design_mfa schema plans subject-list bootstrap without bootstrap_action", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.dmfa_multifer_effect_contract_available())

  sim <- .sim_design_mfa(n_subjects = 3, trials = c(8, 11, 9), p = 4, seed = 94)
  fit <- design_mfa(sim$Y, sim$design, ~ saliency, rank = 1)
  adapter <- as_multifer_adapter(fit)
  expect_null(adapter$bootstrap_action)

  recipe <- multifer::infer_recipe(
    geometry = "adapter",
    relation = "design_contrast",
    design = multifer::clustered_units(unit_role = "Y"),
    targets = "effect_significance",
    effects = "saliency",
    adapter = adapter
  )
  payload <- fit$fit_spec$refit$data
  plan <- multifer:::compile_bootstrap_plan(
    recipe,
    adapter,
    payload,
    original_fit = fit,
    store_aligned_scores = FALSE
  )
  idx <- c(3L, 1L, 3L)
  rep_data <- plan$resample_data(idx)

  expect_equal(plan$data_shape, "adapter_list")
  expect_equal(rep_data$Y, payload$Y[idx])
  expect_equal(rep_data$design, payload$design[idx])
  expect_identical(rep_data$reference_grid, payload$reference_grid)
})

test_that("serial bootstrap uses the statistic_stability sidecar and keeps the public result shape", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.dmfa_multifer_stability_contract_available())

  sim <- .sim_design_mfa(n_subjects = 4, trials = 24, p = 5, seed = 95)
  fit <- design_mfa(sim$Y, sim$design, ~ saliency, rank = 1)
  boot <- infer_design_mfa(
    fit,
    "saliency",
    statistics = c("partial_ss", "singular_values"),
    nrep = 5,
    seed = 96
  )

  expect_identical(boot$engine, "multifer")
  expect_identical(boot$validity, "subject_cluster_bootstrap")
  expect_equal(boot$n_success, 5)
  expect_false("p_value" %in% names(boot$summary$partial_ss))
  expect_true(all(c("observed", "mean", "sd", "lower", "upper") %in%
                    names(boot$summary$partial_ss)))
  expect_equal(nrow(boot$summary$singular_values), 1L)
  expect_true(all(vapply(boot$sampled_subjects, length, integer(1)) == 4L))
})
