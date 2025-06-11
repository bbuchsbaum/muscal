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


context("covstatis: projection and summary methods")

test_that("project_subjects works correctly for known cases", {
  res <- covstatis(subject_data, ncomp = n_comp_fit, dcenter = FALSE, norm_method = "frobenius")
  
  # Test Case 1: Round-trip projection of a training matrix
  # project_cov on a training matrix should recover its partial_scores
  proj_subj1 <- project_subjects(res, subject_data[[1]], subject_ids = "Subj_1")
  expect_equal(proj_subj1$roi_scores$Subj_1, res$partial_scores[[1]], tolerance = 1e-9)
  
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
  
  # RV should be 0, distance should be the matrix's full norm, and scores should be 0
  expect_equal(proj_ortho$scalar_summaries$rv_coefficient, 0, tolerance = 1e-9)
  expect_equal(proj_ortho$scalar_summaries$distance_to_compromise, 
               frobenius_norm(S_ortho_processed), 
               tolerance = 1e-9)
  expect_true(all(abs(proj_ortho$subject_scores) < 1e-9))
}) 