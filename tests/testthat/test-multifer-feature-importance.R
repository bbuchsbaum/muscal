library(testthat)
library(muscal)

test_that("MCCA-family multifer adapters expose squared-loading feature importance", {
  skip_if_not_installed("multifer")

  set.seed(901)
  n <- 18
  X1 <- matrix(rnorm(n * 7), n, 7)
  X2 <- matrix(rnorm(n * 6), n, 6)
  Y <- matrix(rnorm(n * 5), n, 5)
  colnames(X1) <- paste0("x1_", seq_len(ncol(X1)))
  colnames(X2) <- paste0("x2_", seq_len(ncol(X2)))
  colnames(Y) <- paste0("y_", seq_len(ncol(Y)))

  mcca_fit <- mcca(list(X1 = X1, X2 = X2), ncomp = 3, ridge = 1e-6)
  aligned_fit <- aligned_mcca(
    list(X1 = X1, X2 = X2),
    row_index = list(X1 = seq_len(n), X2 = seq_len(n)),
    N = n,
    ncomp = 3,
    ridge = 1e-6
  )
  anchored_fit <- anchored_mcca(
    Y,
    list(X1 = X1, X2 = X2),
    row_index = list(X1 = seq_len(n), X2 = seq_len(n)),
    ncomp = 3,
    ridge = 1e-6
  )

  for (fit in list(mcca_fit, aligned_fit, anchored_fit)) {
    adapter <- as_multifer_adapter(fit)
    units <- multifer::form_units(adapter$roots(fit), selected = c(TRUE, TRUE, FALSE))
    report <- multifer::check_feature_importance_adapter(
      adapter = adapter,
      fit = fit,
      data = fit$fit_spec$refit$data,
      units = units,
      statistic = "both",
      normalize = "none"
    )

    expect_s3_class(report, "multifer_feature_importance_adapter_check")
    expect_true(all(report$passed))

    from_loadings <- multifer::feature_importance_from_fit(
      fit = fit,
      adapter = adapter,
      units = units,
      data = fit$fit_spec$refit$data,
      k = "selected",
      scope = "both",
      statistic = "loadings",
      normalize = "none"
    )
    from_adapter <- multifer::feature_importance_from_fit(
      fit = fit,
      adapter = adapter,
      units = units,
      data = fit$fit_spec$refit$data,
      k = "selected",
      scope = "both",
      statistic = "adapter",
      normalize = "none"
    )

    expect_equal(from_adapter$estimate, from_loadings$estimate, tolerance = 1e-10)
    expect_true(all(from_adapter$estimate >= 0))
  }
})

test_that("subspace feature importance is invariant to rotations inside MCCA units", {
  skip_if_not_installed("multifer")

  set.seed(902)
  n <- 16
  X1 <- matrix(rnorm(n * 5), n, 5)
  X2 <- matrix(rnorm(n * 4), n, 4)
  fit <- mcca(list(X1 = X1, X2 = X2), ncomp = 2, ridge = 1e-6)
  adapter <- as_multifer_adapter(fit)
  units <- multifer::infer_units(
    unit_id = "sub",
    unit_type = "subspace",
    members = list(1:2),
    identifiable = FALSE,
    selected = TRUE
  )

  theta <- pi / 5
  Q <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  fit_rot <- fit
  fit_rot$scaled_loadings <- fit$scaled_loadings %*% Q
  fit_rot$cor_loadings <- fit$cor_loadings %*% Q
  fit_rot$scaled_loadings_by_block <- lapply(fit$scaled_loadings_by_block, function(L) L %*% Q)

  base <- multifer::feature_importance_from_fit(
    fit, adapter, units = units, k = "selected", statistic = "loadings",
    scope = "both", normalize = "none"
  )
  rotated <- multifer::feature_importance_from_fit(
    fit_rot, adapter, units = units, k = "selected", statistic = "loadings",
    scope = "both", normalize = "none"
  )

  expect_equal(rotated$estimate, base$estimate, tolerance = 1e-10)
})

test_that("MCCA-family multifer adapters support feature-importance p-values", {
  skip_if_not_installed("multifer")

  set.seed(904)
  n <- 12
  X1 <- matrix(rnorm(n * 4), n, 4)
  X2 <- matrix(rnorm(n * 5), n, 5)
  Y <- matrix(rnorm(n * 3), n, 3)
  fits <- list(
    mcca = mcca(list(X1 = X1, X2 = X2), ncomp = 2, ridge = 1e-6),
    aligned_mcca = aligned_mcca(
      list(X1 = X1, X2 = X2),
      row_index = list(X1 = seq_len(n), X2 = seq_len(n)),
      N = n,
      ncomp = 2,
      ridge = 1e-6
    ),
    anchored_mcca = anchored_mcca(
      Y,
      list(X1 = X1, X2 = X2),
      row_index = list(X1 = seq_len(n), X2 = seq_len(n)),
      ncomp = 2,
      ridge = 1e-6
    )
  )

  for (fit in fits) {
    adapter <- as_multifer_adapter(fit)
    geometry <- if (inherits(fit, "mcca") && !inherits(fit, "aligned_mcca")) {
      "multiblock"
    } else {
      "adapter"
    }
    recipe <- multifer::infer_recipe(
      geometry = geometry,
      relation = "correlation",
      targets = "variable_stability",
      adapter = adapter
    )
    units <- multifer::form_units(adapter$roots(fit), selected = c(TRUE, FALSE))
    out <- multifer::feature_importance_pvalues(
      recipe = recipe,
      adapter = adapter,
      data = fit$fit_spec$refit$data,
      units = units,
      original_fit = fit,
      B = 2L,
      statistic = "adapter"
    )

    expect_s3_class(out, "multifer_feature_importance_pvalues")
    expect_true(all(out$p_value >= 0 & out$p_value <= 1))
  }
})

test_that("existing muscal multifer adapters expose adapter-owned feature statistics", {
  skip_if_not_installed("multifer")

  set.seed(903)
  subjects <- factor(rep(paste0("S", 1:4), each = 12))
  y <- factor(rep(rep(c("A", "B", "C"), each = 4), times = 4))
  class_signal <- model.matrix(~ y - 1)
  subject_shift <- model.matrix(~ subjects - 1) %*% matrix(rnorm(4 * 6, sd = 0.2), 4, 6)
  load <- matrix(c(
    1.2, 0.1, -0.4, 0.0, 0.2, 0.0,
   -0.8, 0.5,  0.3, 0.1, 0.0, 0.2,
    0.1, 1.1,  0.2, 0.4, 0.3, 0.0
  ), 3, 6, byrow = TRUE)
  X <- class_signal %*% load + subject_shift +
    matrix(rnorm(length(y) * 6, sd = 0.25), length(y), 6)
  md <- multidesign::multidesign(X, data.frame(y = y, subj_id = subjects))
  fit <- bada(md, y = y, subject = subj_id, ncomp = 2, preproc = multivarious::pass())
  blocks <- muscal:::.bada_subject_barycenter_blocks(fit, md)
  adapter <- as_multifer_adapter(fit)
  adapter_fit <- muscal:::.bada_fit_from_barycenter_blocks(blocks, reference = fit)
  units <- multifer::form_units(adapter$roots(adapter_fit), selected = c(TRUE, FALSE))

  report <- multifer::check_feature_importance_adapter(
    adapter = adapter,
    fit = adapter_fit,
    data = blocks,
    units = units,
    statistic = "both",
    normalize = "none"
  )

  expect_true(all(report$passed))
})
