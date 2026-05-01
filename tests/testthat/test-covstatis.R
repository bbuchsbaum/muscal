library(testthat)
library(muscal)
library(multivarious)

# Helper to create a random symmetric matrix
make_sym_matrix <- function(p, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- matrix(rnorm(p * p), p, p)
  (m + t(m)) / 2
}

# Helper for Frobenius norm
frobenius_norm <- function(m) {
  sqrt(sum(m^2))
}

context("covstatis: helper functions")

test_that("double_center works correctly", {
  M <- matrix(1:9, 3, 3)
  M_dc <- muscal:::double_center(M)
  
  # Row and column sums should be zero
  expect_true(all(abs(rowSums(M_dc)) < 1e-12))
  expect_true(all(abs(colSums(M_dc)) < 1e-12))
  
  # It should be idempotent
  expect_equal(M_dc, muscal:::double_center(M_dc), tolerance = 1e-12)
})

test_that("norm_crossprod normalizes to unit Frobenius norm", {
  S <- make_sym_matrix(10)
  S_norm <- muscal:::norm_crossprod(S)
  expect_equal(frobenius_norm(S_norm), 1)
})

test_that("compute_prodmat is correct", {
  S1 <- make_sym_matrix(5, seed = 1)
  S2 <- make_sym_matrix(5, seed = 2)
  S_list <- list(S1 = S1, S2 = S2)
  
  # Test fast vs slow paths
  C_fast <- muscal:::compute_prodmat(S_list, fast = TRUE, normalize = FALSE)
  C_slow <- muscal:::compute_prodmat(S_list, fast = FALSE, normalize = FALSE)
  expect_equal(C_fast, C_slow)
  
  # Test with normalization
  S_list_norm <- lapply(S_list, muscal:::norm_crossprod)
  C_norm <- muscal:::compute_prodmat(S_list_norm, normalize = TRUE)
  expect_equal(diag(C_norm), c(1, 1))
  
  # Test without normalization
  C_no_norm <- muscal:::compute_prodmat(S_list, normalize = FALSE)
  expect_equal(diag(C_no_norm), c(sum(S1^2), sum(S2^2)))
})

context("covstatis: main function and properties")

# --- Setup for main validation tests ---
p <- 30
n_subjects <- 10
n_comp_truth <- 5
n_comp_fit <- 3

# 1. Create a ground-truth low-rank matrix
set.seed(123)
V_truth <- qr.Q(qr(matrix(rnorm(p * n_comp_truth), p, n_comp_truth)))
lambda_truth <- (p - 1:n_comp_truth)^2
S_truth <- V_truth %*% diag(lambda_truth) %*% t(V_truth)

# 2. Create subject matrices as noisy versions of the truth
noise_level <- 0.01
subject_data <- lapply(1:n_subjects, function(i) {
  # start from a common orthogonal basis ...
  B  <- V_truth
  # ... then apply a SMALL random orthogonal rotation
  set.seed(i)
  Q_rand <- matrix(rnorm(n_comp_truth * n_comp_truth, sd = 0.1), n_comp_truth, n_comp_truth)
  Q <- qr.Q(qr(diag(n_comp_truth) + Q_rand))
  V_i <- B %*% Q
  
  noise <- make_sym_matrix(p, seed=i) * noise_level * frobenius_norm(S_truth)
  (V_i %*% diag(lambda_truth) %*% t(V_i)) + noise
})
names(subject_data) <- paste0("Subj_", 1:n_subjects)
# --- End Setup ---

test_that("covstatis model output has correct structure and properties", {
  res <- covstatis(subject_data, ncomp = n_comp_fit, norm_method = "none", dcenter = FALSE)
  
  expect_s3_class(res, "covstatis")
  expect_s3_class(res, "projector")
  
  # Check dimensions
  expect_equal(ncomp(res), n_comp_fit)
  expect_equal(nrow(coefficients(res)), p)
  expect_equal(ncol(coefficients(res)), n_comp_fit)
  
  # Check orthogonality of basis vectors
  expect_equal(crossprod(coefficients(res)), diag(n_comp_fit), tolerance = 1e-9)
  
  # Check alpha weights
  expect_true(all(res$alpha >= 0))
  expect_equal(sum(res$alpha), 1)
  
  # Partial scores
  expect_equal(length(res$partial_scores), n_subjects)
  expect_equal(dim(res$partial_scores[[1]]), c(p, n_comp_fit))
})

test_that("covstatis recovers the ground truth structure with low noise", {
  res <- covstatis(subject_data, ncomp = n_comp_fit, norm_method = "none", dcenter = FALSE)
  
  # The principal angles between the true subspace and fitted subspace should be near zero.
  # This is a stronger test than comparing the compromise matrices.
  proj_V <- t(coefficients(res)) %*% V_truth[, 1:n_comp_fit]
  cos_angles <- svd(proj_V)$d # cosines of principal angles
  
  # All cosines should be very close to 1 (i.e., angles close to 0)
  expect_true(all(cos_angles > 0.99))
})

test_that("dcenter=TRUE makes results invariant to constant shifts", {
  # Add a large constant to each matrix
  subject_data_shifted <- lapply(subject_data, function(m) m + 100)
  
  res_orig <- covstatis(subject_data, ncomp = n_comp_fit, dcenter = TRUE, norm_method = "frobenius")
  res_shifted <- covstatis(subject_data_shifted, ncomp = n_comp_fit, dcenter = TRUE, norm_method = "frobenius")
  
  # The key outputs should be identical up to sign flips in eigenvectors
  expect_equal(res_orig$C, res_shifted$C, tolerance = 1e-9)
  expect_equal(res_orig$alpha, res_shifted$alpha, tolerance = 1e-9)
  expect_equal(abs(coefficients(res_orig)), abs(coefficients(res_shifted)), tolerance = 1e-9)
  expect_equal(res_orig$sdev, res_shifted$sdev, tolerance = 1e-9)
})

test_that("default subject scores use nondegenerate RV-space coordinates", {
  res <- covstatis(subject_data, ncomp = n_comp_fit)

  G <- muscal:::.get_table_scores(res)
  T_scores <- muscal:::rv_subject_scores(res)

  expect_equal(G, T_scores, tolerance = 1e-12)
  expect_true(max(abs(G)) > 1e-6)
})


context("covstatis: projection and summary methods")

test_that("project_subjects works correctly for known cases", {
  res <- covstatis(subject_data, ncomp = n_comp_fit, dcenter = FALSE, norm_method = "frobenius")
  
  # Test Case 1: Round-trip projection of a training matrix
  # project_cov on a training matrix should recover its partial_scores
  proj_subj1 <- project_subjects(res, subject_data[[1]], subject_ids = "Subj_1")
  expect_equal(proj_subj1$roi_scores$Subj_1, res$partial_scores[[1]], tolerance = 1e-9)
  expect_equal(
    unname(drop(proj_subj1$subject_scores[1, ])),
    unname(drop(muscal:::rv_subject_scores(res)[1, ])),
    tolerance = 1e-9
  )
  
  # Test Case 2: Projecting a matrix that lies perfectly in the compromise space
  S_compromise_processed <- reconstruct(res)
  # Create a "raw" matrix that will become the processed compromise after normalization
  S_compromise_unnorm <- S_compromise_processed * 100
  
  proj_compromise <- project_subjects(res, S_compromise_unnorm)
  
  # Distance should be zero and RV should be 1
  expect_equal(proj_compromise$scalar_summaries$distance_to_compromise, 0, tolerance = 1e-9)
  expect_equal(proj_compromise$scalar_summaries$rv_coefficient, 1, tolerance = 1e-9)
  
  # Test Case 3: Projecting a matrix orthogonal to the compromise space
  U <- coefficients(res)
  # Find an orthonormal basis for the space orthogonal to U's columns
  Q_all <- qr.Q(qr(cbind(U, matrix(rnorm(p * (p-ncomp(res))), p, p-ncomp(res)))))
  U_ortho_basis <- Q_all[, (ncomp(res) + 1):p]
  
  # Create a matrix in this orthogonal space
  S_ortho_raw <- U_ortho_basis %*% diag(rnorm(p - ncomp(res))) %*% t(U_ortho_basis)
  
  proj_ortho <- project_subjects(res, S_ortho_raw)
  S_ortho_processed <- muscal:::.pre_process_new_cov(res, S_ortho_raw)
  
  # RV should be 0, distance should be the matrix's full norm, and ROI scores
  # should vanish because the matrix is orthogonal to the retained compromise basis.
  expect_equal(proj_ortho$scalar_summaries$rv_coefficient, 0, tolerance = 1e-9)
  expect_equal(proj_ortho$scalar_summaries$distance_to_compromise, 
               frobenius_norm(S_ortho_processed), 
               tolerance = 1e-9)
  expect_true(all(abs(proj_ortho$roi_scores[[1]]) < 1e-9))
})

test_that("project_table_covariate uses table scores that remain informative with double centering", {
  res <- covstatis(subject_data, ncomp = n_comp_fit)
  y <- seq_along(subject_data)

  dim_cos <- project_table_covariate(res, y, what = "dimension", scale = "cosine")
  T_scores <- muscal:::rv_subject_scores(res)
  manual_cos <- as.numeric(t(T_scores) %*% y) /
    (sqrt(colSums(T_scores^2)) * sqrt(sum(y^2)))
  names(manual_cos) <- paste0("Dim", seq_along(manual_cos))

  expect_equal(dim_cos, manual_cos, tolerance = 1e-10)
  expect_true(any(abs(dim_cos) > 1e-6))
})

test_that("project_feature_covariate aligns with feature-space scores", {
  res <- covstatis(subject_data, ncomp = n_comp_fit)
  z <- seq_len(nrow(res$roi_scores))

  dim_cos <- project_feature_covariate(res, z, scale = "cosine")
  F_scores <- res$roi_scores
  manual_cos <- as.numeric(t(F_scores) %*% z) /
    (sqrt(colSums(F_scores^2)) * sqrt(sum(z^2)))
  names(manual_cos) <- paste0("Dim", seq_along(manual_cos))

  expect_equal(dim_cos, manual_cos, tolerance = 1e-10)
  expect_true(any(abs(dim_cos) > 1e-6))
})

test_that("project_covariate remains backward compatible for table-side input", {
  res <- covstatis(subject_data, ncomp = n_comp_fit)
  y <- seq_along(subject_data)

  expect_equal(
    project_covariate(res, y, what = "dimension", scale = "cosine"),
    project_table_covariate(res, y, what = "dimension", scale = "cosine")
  )
  expect_equal(
    project_covariate(res, y, what = "observation", scale = "beta"),
    project_table_covariate(res, y, what = "feature", scale = "beta")
  )
})

test_that("project_subjects round-trips subject scores under default preprocessing", {
  res <- covstatis(subject_data, ncomp = n_comp_fit)
  proj_subj1 <- project_subjects(res, subject_data[[1]], subject_ids = "Subj_1")

  expect_equal(
    unname(drop(proj_subj1$subject_scores[1, ])),
    unname(drop(muscal:::rv_subject_scores(res)[1, ])),
    tolerance = 1e-9
  )
  expect_true(any(abs(proj_subj1$subject_scores[1, ]) > 1e-6))
})

context("covstatis: multifer component tests")

test_that("COVSTATIS multifer payload preserves processed table shape", {
  res <- covstatis(subject_data, ncomp = n_comp_fit)
  payload <- muscal:::.covstatis_multifer_payload(res, subject_data)

  expect_true(muscal:::.covstatis_valid_payload(payload))
  expect_length(payload$tables, length(subject_data))
  expect_equal(unique(vapply(payload$tables, nrow, integer(1))), p)
  expect_equal(unique(vapply(payload$tables, ncol, integer(1))), p)
})

test_that("COVSTATIS multifer adapters declare the expected component contract", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.covstatis_multifer_adapter_geometry_available())

  for (axes in c("compromise", "interstructure")) {
    adapter <- muscal:::.covstatis_multifer_adapter(axes)
    expect_true(multifer::adapter_supports(
      adapter, "adapter", "variance", "component_significance"
    ))
    expect_equal(adapter$validity_level, "conditional")
    expect_true("roi_permutation_null" %in% adapter$declared_assumptions)
  }
})

test_that("COVSTATIS multifer statistics match adapter-fit root ratios", {
  res <- covstatis(subject_data, ncomp = n_comp_fit)
  payload <- muscal:::.covstatis_multifer_payload(res, subject_data)

  for (axes in c("compromise", "interstructure")) {
    fit <- muscal:::.covstatis_adapter_fit(payload, axes = axes)
    oracle <- fit$roots[[1L]] / sum(fit$roots)

    expect_equal(
      muscal:::.covstatis_leading_root_ratio(payload, axes = axes),
      oracle,
      tolerance = 1e-12
    )
    expect_true(all(fit$roots >= 0))
    expect_true(all(is.finite(fit$v)))
  }
})

test_that("COVSTATIS multifer null and residual payloads remain valid", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.covstatis_multifer_adapter_geometry_available())

  res <- covstatis(subject_data, ncomp = n_comp_fit)
  payload <- muscal:::.covstatis_multifer_payload(res, subject_data)

  for (axes in c("compromise", "interstructure")) {
    adapter <- muscal:::.covstatis_multifer_adapter(axes)
    fit <- muscal:::.covstatis_adapter_fit(payload, axes = axes)
    residual <- adapter$residualize(fit, 1L, payload)
    null_payload <- adapter$null_action(fit, payload)

    expect_true(muscal:::.covstatis_valid_payload(residual))
    expect_true(muscal:::.covstatis_valid_payload(null_payload))
    expect_named(residual$tables, names(payload$tables))
    expect_named(null_payload$tables, names(payload$tables))
  }
})

test_that("COVSTATIS multifer statistics are invariant to table order", {
  res <- covstatis(subject_data, ncomp = n_comp_fit)
  payload <- muscal:::.covstatis_multifer_payload(res, subject_data)
  reversed <- payload
  reversed$tables <- rev(payload$tables)

  for (axes in c("compromise", "interstructure")) {
    expect_equal(
      muscal:::.covstatis_leading_root_ratio(reversed, axes = axes),
      muscal:::.covstatis_leading_root_ratio(payload, axes = axes),
      tolerance = 1e-10
    )
  }
})

test_that("COVSTATIS multifer payload validation rejects contaminated tables", {
  res <- covstatis(subject_data, ncomp = n_comp_fit)
  payload <- muscal:::.covstatis_multifer_payload(res, subject_data)

  bad <- payload
  bad$tables[[1]][1, 2] <- bad$tables[[1]][1, 2] + 1
  expect_false(muscal:::.covstatis_valid_payload(bad))

  bad <- payload
  bad$tables[[1]][1, 1] <- NA_real_
  expect_false(muscal:::.covstatis_valid_payload(bad))
})

test_that("infer_covstatis tests compromise ROI eigencomponents with multifer", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.covstatis_multifer_adapter_geometry_available())

  res <- covstatis(subject_data, ncomp = n_comp_fit)
  out <- infer_covstatis(
    res,
    data = subject_data,
    axes = "compromise",
    B = 5L,
    seed = 6101,
    mc_batch_size = 5L
  )

  expect_s3_class(out, "infer_result")
  expect_equal(out$provenance$adapter_id, "muscal_covstatis_compromise")
  expect_gt(nrow(out$component_tests), 0L)
  expect_match(out$mc$statistic_label, "adapter component_stat")
})

test_that("infer_covstatis tests interstructure RV axes with multifer", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.covstatis_multifer_adapter_geometry_available())

  res <- covstatis(subject_data, ncomp = n_comp_fit)
  out <- infer_covstatis(
    res,
    data = subject_data,
    axes = "interstructure",
    B = 5L,
    seed = 6102,
    mc_batch_size = 5L
  )

  expect_s3_class(out, "infer_result")
  expect_equal(out$provenance$adapter_id, "muscal_covstatis_interstructure")
  expect_gt(nrow(out$component_tests), 0L)
})

test_that("infer_covstatis can return both COVSTATIS significance families", {
  skip_if_not_installed("multifer")
  skip_if_not(muscal:::.covstatis_multifer_adapter_geometry_available())

  res <- covstatis(subject_data, ncomp = n_comp_fit)
  out <- infer_covstatis(
    res,
    data = subject_data,
    axes = "both",
    B = 5L,
    seed = 6103,
    mc_batch_size = 5L
  )

  expect_s3_class(out, "covstatis_multifer_result")
  expect_named(out, c("interstructure", "compromise"))
  expect_s3_class(out$interstructure, "infer_result")
  expect_s3_class(out$compromise, "infer_result")
})
