library(testthat)
library(muscal)
library(multidesign)

skip_if_not_installed("multidesign")

# --- bada.multidesign basic tests --------------------------------------------

test_that("bada.multidesign works with basic splitting", {
  set.seed(123)
  x <- matrix(rnorm(20), nrow = 10)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = 5),
    y = rep(c("A", "B"), 5)
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 1)

  expect_s3_class(res, "bada")
  expect_equal(levels(res$subjects), c("S1", "S2"))
  expect_equal(levels(res$labels), c("A", "B"))
  expect_length(res$block_indices, 2)
  expect_true(multivarious::ncomp(res) >= 1)
})

test_that("bada.multidesign produces valid projector", {
  set.seed(42)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 5), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 2)

  expect_s3_class(res, "bada")
  expect_true(inherits(res, "projector"))
  expect_true(!is.null(res$v))
  expect_true(!is.null(res$s))
  expect_true(!is.null(res$sdev))
  expect_equal(nrow(res$s), n_obs)
})

test_that("bada.multidesign with multiple classes", {
  set.seed(1)
  n_obs <- 30
  x <- matrix(rnorm(n_obs * 6), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2", "S3"), each = n_obs / 3),
    y = factor(rep(c("A", "B", "C"), n_obs / 3))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 2)

  expect_s3_class(res, "bada")
  expect_equal(length(levels(res$labels)), 3)
  expect_equal(length(levels(res$subjects)), 3)
  # Should have barycenters for each class
  expect_equal(nrow(res$barycenters), 3)
})

# --- bada.multidesign with residual analysis ---------------------------------

test_that("bada.multidesign with residual components", {
  set.seed(99)
  n_obs <- 40
  x <- matrix(rnorm(n_obs * 10), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  # Using rescomp > 0 enables residual analysis
  res <- bada(md, y = y, subject = subj_id, ncomp = 1, resdim = 5, rescomp = 2)

  expect_s3_class(res, "bada")
  expect_equal(res$resdim, 5)
  expect_equal(res$rescomp, 2)
  # Total components = ncomp (class) + rescomp (residual)
  expect_equal(ncol(res$v), 3)  # 1 + 2
})

# --- reprocess.bada tests ----------------------------------------------------

test_that("reprocess.bada with no block applies averaged preprocessing", {
  set.seed(123)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 8), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 1)

  # Create new data with same dimensions
  new_data <- matrix(rnorm(5 * 8), nrow = 5)
  reprocessed <- reprocess(res, new_data)

  expect_true(is.matrix(reprocessed))
  expect_equal(nrow(reprocessed), 5)
  expect_equal(ncol(reprocessed), ncol(x))
})

test_that("reprocess.bada with block uses block-specific preprocessing", {
  set.seed(456)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 8), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 1)

  new_data <- matrix(rnorm(5 * 8), nrow = 5)
  reprocessed <- reprocess(res, new_data, block = "S1")

  expect_true(is.matrix(reprocessed))
  expect_equal(nrow(reprocessed), 5)
})

test_that("reprocess.bada validates block name", {
  set.seed(789)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 8), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 1)
  new_data <- matrix(rnorm(5 * 8), nrow = 5)

  expect_error(
    reprocess(res, new_data, block = "INVALID"),
    class = "error"
  )
})

test_that("reprocess.bada validates data dimensions", {
  set.seed(111)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 8), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 1)

  # Wrong number of columns should error
  new_data_wrong <- matrix(rnorm(5 * 3), nrow = 5)  # 3 cols instead of 8
  expect_error(reprocess(res, new_data_wrong))
})

# --- project.bada tests ------------------------------------------------------

test_that("project.bada projects new data", {
  set.seed(222)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 8), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 2)

  new_data <- matrix(rnorm(5 * 8), nrow = 5)
  projected <- project(res, new_data)

  expect_true(is.matrix(projected))
  expect_equal(nrow(projected), 5)
  # ncomp may be capped by number of classes (2) minus 1
  expect_true(ncol(projected) >= 1)
})

test_that("project.bada with block parameter", {
  set.seed(333)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 8), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 2)

  new_data <- matrix(rnorm(5 * 8), nrow = 5)
  projected <- project(res, new_data, block = "S2")

  expect_true(is.matrix(projected))
  expect_equal(nrow(projected), 5)
  # ncomp may be capped by number of classes
  expect_true(ncol(projected) >= 1)
})

# --- summarize_boot internal function ----------------------------------------

test_that("summarize_boot computes correct statistics", {
  set.seed(100)

  # Create mock bootstrap result structure
  boot_ret <- data.frame(i = 1:5)
  boot_ret$boot_i <- lapply(1:5, function(x) matrix(rnorm(6), 2, 3))
  boot_ret$boot_j <- lapply(1:5, function(x) matrix(rnorm(6), 2, 3))

  result <- muscal:::summarize_boot(boot_ret, alpha = 0.05)

  expect_type(result, "list")
  expect_named(result, c(
    "boot_scores_mean", "boot_scores_sd", "boot_scores_upper", "boot_scores_lower",
    "boot_lds_mean", "boot_lds_sd", "boot_lds_upper", "boot_lds_lower"
  ))

  # Check dimensions
  expect_equal(dim(result$boot_scores_mean), c(2, 3))
  expect_equal(dim(result$boot_lds_mean), c(2, 3))

  # Check that sd is non-negative
  expect_true(all(result$boot_scores_sd >= 0))
  expect_true(all(result$boot_lds_sd >= 0))
})

# --- summarize_matrix internal function --------------------------------------

test_that("summarize_matrix applies function across matrices", {
  matrices <- list(
    matrix(1:4, 2, 2),
    matrix(5:8, 2, 2),
    matrix(9:12, 2, 2)
  )

  result <- muscal:::summarize_matrix(matrices, mean)

  expected <- matrix(c(5, 6, 7, 8), 2, 2)
  expect_equal(result, expected)
})

test_that("summarize_matrix with sd function", {
  set.seed(42)
  matrices <- lapply(1:10, function(x) matrix(rnorm(4), 2, 2))

  result <- muscal:::summarize_matrix(matrices, sd)

  expect_equal(dim(result), c(2, 2))
  expect_true(all(result >= 0))
})

# --- bootstrap.bada test (simplified) ----------------------------------------

test_that("bootstrap.bada requires data argument", {
  set.seed(444)
  n_obs <- 20
  x <- matrix(rnorm(n_obs * 6), nrow = n_obs)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = n_obs / 2),
    y = factor(rep(c("A", "B"), n_obs / 2))
  )
  md <- multidesign::multidesign(x, design)

  res <- bada(md, y = y, subject = subj_id, ncomp = 1)

  expect_error(
    multivarious::bootstrap(res, nboot = 10),
    "data.*required"
  )
})

# --- multifer inference adapter ----------------------------------------------

make_bada_multifer_fixture <- function() {
  set.seed(5501)
  subjects <- factor(rep(paste0("S", 1:4), each = 12))
  y <- factor(rep(rep(c("A", "B", "C"), each = 4), times = 4))
  class_signal <- model.matrix(~ y - 1)
  subject_shift <- model.matrix(~ subjects - 1) %*% matrix(rnorm(4 * 6, sd = 0.2), 4, 6)
  load <- matrix(c(
    1.2, 0.1, -0.4, 0.0, 0.2, 0.0,
   -0.8, 0.5,  0.3, 0.1, 0.0, 0.2,
    0.1, 1.1,  0.2, 0.4, 0.3, 0.0
  ), 3, 6, byrow = TRUE)
  x <- class_signal %*% load + subject_shift + matrix(rnorm(length(y) * 6, sd = 0.25), length(y), 6)
  design <- data.frame(subj_id = subjects, y = y)
  md <- multidesign::multidesign(x, design)
  list(data = md, fit = bada(md, y = y, subject = subj_id, ncomp = 2))
}

test_that("BaDA reconstructs subject barycenter blocks for multifer", {
  fixture <- make_bada_multifer_fixture()
  blocks <- muscal:::.bada_subject_barycenter_blocks(fixture$fit, fixture$data)

  expect_length(blocks, 4)
  expect_true(all(vapply(blocks, is.matrix, logical(1))))
  expect_equal(unique(vapply(blocks, nrow, integer(1))), 3L)
  expect_equal(unique(vapply(blocks, ncol, integer(1))), 6L)
  expect_equal(rownames(blocks[[1]]), fixture$fit$label_set)
})

test_that("infer_bada delegates BaDA inference to multifer multiblock hooks", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.bada_multifer_contract_available())

  fixture <- make_bada_multifer_fixture()
  res <- infer_bada(
    fixture$fit,
    data = fixture$data,
    B = 5L,
    R = 3L,
    seed = 5502,
    mc_batch_size = 5L
  )

  expect_s3_class(res, "infer_result")
  expect_equal(res$provenance$adapter_id, "muscal_bada")
  expect_gt(nrow(res$component_tests), 0L)
  expect_true(all(unique(res$variable_stability$domain) == "global"))
  expect_true(all(unique(res$score_stability$domain) == "global"))
})

test_that("BaDA multifer adapter declares and enforces its hook contract", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.bada_multifer_contract_available())

  fixture <- make_bada_multifer_fixture()
  adapter <- muscal:::.bada_multifer_adapter()
  blocks <- muscal:::.bada_subject_barycenter_blocks(fixture$fit, fixture$data)

  expect_true(multifer::adapter_supports(
    adapter, "multiblock", "variance", "component_significance"
  ))
  expect_true(multifer::adapter_supports(
    adapter, "multiblock", "variance", "score_stability"
  ))
  expect_equal(adapter$domains(fixture$fit, blocks), "global")
  expect_error(adapter$scores(fixture$fit, domain = "subject"), "global")

  check <- adapter$checked_assumptions[[1]]$check
  expect_true(check(blocks))
  bad_blocks <- blocks
  bad_blocks[[1]][1, 1] <- NA_real_
  expect_false(check(bad_blocks))
})

test_that("BaDA multifer statistic matches the barycenter SVD oracle", {
  fixture <- make_bada_multifer_fixture()
  blocks <- muscal:::.bada_subject_barycenter_blocks(fixture$fit, fixture$data)
  fit <- muscal:::.bada_fit_from_barycenter_blocks(blocks, reference = fixture$fit)

  Xc <- Reduce("+", blocks) / length(blocks)
  roots <- svd(Xc, nu = 0L, nv = 0L)$d^2
  oracle <- roots[[1L]] / sum(roots)

  expect_equal(muscal:::.bada_leading_root_ratio(blocks), oracle, tolerance = 1e-10)
  expect_equal(fit$roots, roots[seq_len(ncol(fit$v))], tolerance = 1e-10)
  expect_true(all(is.finite(fit$v)))
  expect_true(all(is.finite(fit$fscores)))
})

test_that("BaDA multifer residual and null actions preserve payload validity", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.bada_multifer_contract_available())

  fixture <- make_bada_multifer_fixture()
  adapter <- muscal:::.bada_multifer_adapter()
  blocks <- muscal:::.bada_subject_barycenter_blocks(fixture$fit, fixture$data)

  stat_before <- adapter$component_stat(fixture$fit, blocks, 1L)
  residual <- adapter$residualize(fixture$fit, 1L, blocks)
  null_blocks <- adapter$null_action(fixture$fit, blocks)
  residual_average <- Reduce("+", residual) / length(residual)
  removed_direction <- fixture$fit$v[, 1L, drop = FALSE]

  expect_true(muscal:::.bada_valid_barycenter_blocks(residual))
  expect_true(muscal:::.bada_valid_barycenter_blocks(null_blocks))
  expect_named(residual, names(blocks))
  expect_named(null_blocks, names(blocks))
  expect_true(is.finite(stat_before))
  expect_equal(unname(residual_average %*% removed_direction),
               matrix(0, nrow(residual_average), 1L),
               tolerance = 1e-8)
})

test_that("BaDA multifer leading statistic is invariant to subject-block order", {
  fixture <- make_bada_multifer_fixture()
  blocks <- muscal:::.bada_subject_barycenter_blocks(fixture$fit, fixture$data)
  reversed <- rev(blocks)

  expect_equal(
    muscal:::.bada_leading_root_ratio(reversed),
    muscal:::.bada_leading_root_ratio(blocks),
    tolerance = 1e-12
  )
})

test_that("BaDA multifer bootstrap uses subject-block resampling and projected scores", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.bada_multifer_contract_available())

  fixture <- make_bada_multifer_fixture()
  adapter <- muscal:::.bada_multifer_adapter()
  blocks <- muscal:::.bada_subject_barycenter_blocks(fixture$fit, fixture$data)
  recipe <- multifer::infer_recipe(
    geometry = "multiblock",
    relation = "variance",
    adapter = adapter,
    targets = c("variable_stability", "score_stability", "subspace_stability")
  )
  units <- multifer::form_units(adapter$roots(fixture$fit))

  artifact <- multifer::bootstrap_fits(
    recipe = recipe,
    adapter = adapter,
    data = blocks,
    original_fit = fixture$fit,
    units = units,
    R = 2L,
    seed = 5503
  )

  expect_s3_class(artifact, "multifer_bootstrap_artifact")
  expect_true(artifact$used_bootstrap_action)
  expect_equal(artifact$score_source, "project_scores")
  expect_equal(artifact$domains, "global")
  expect_length(artifact$reps[[1]]$resample_indices, length(blocks))
})

test_that("infer_bada refuses residual BaDA components for now", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.bada_multifer_contract_available())

  set.seed(5504)
  n_obs <- 48
  x <- matrix(rnorm(n_obs * 8), nrow = n_obs)
  design <- data.frame(
    subj_id = factor(rep(c("S1", "S2", "S3", "S4"), each = 12)),
    y = factor(rep(rep(c("A", "B", "C"), each = 4), times = 4))
  )
  md <- multidesign::multidesign(x, design)
  fit <- bada(md, y = y, subject = subj_id, ncomp = 1, resdim = 3, rescomp = 1)

  expect_error(
    infer_bada(fit, data = md, B = 3L, R = 1L),
    "rescomp = 0"
  )
})
