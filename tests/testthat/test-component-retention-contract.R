.common_loading_fixture <- function() {
  V <- matrix(c(1, 0, 0, 0, 0), 5, 1,
               dimnames = list(paste0("f", 1:5), "component1"))
  list(V_list = list(V, V), sdev = 3)
}

test_that("heterogeneous loadings cannot produce cross-space coherence evidence", {
  fit <- list(V_list = list(matrix(c(1, 0, 0), 3, 1),
                            matrix(c(1, 0, 0, 0, 0, 0), 6, 1)), sdev = 1)
  result <- significant_components(fit, n = 20, noise_variance = 1,
    spectrum = "covariance", noise_model = "isotropic")
  expect_true(all(is.na(result$icc)))
  expect_true(all(is.na(result$icc_pvalue)))
  expect_identical(result$icc_status, "disabled")
  expect_error(significant_components(fit, n = 20, check_icc = TRUE,
    common_features = paste0("f", 1:3)), "common feature domain")
})

test_that("ICC requires an explicit common feature identity contract", {
  fit <- .common_loading_fixture()
  ids <- rownames(fit$V_list[[1]])
  expect_error(significant_components(fit, n = 20, check_icc = TRUE),
               "common_features")
  fit$V_list[[2]] <- `rownames<-`(fit$V_list[[2]], paste0("other", 1:5))
  expect_error(significant_components(fit, n = 20, check_icc = TRUE,
    common_features = ids), "common feature domain")
})

test_that("matched ICC has analytic cosines and is feature-order invariant", {
  fit <- .common_loading_fixture()
  ids <- rownames(fit$V_list[[1]])
  same <- suppressMessages(significant_components(fit, n = 20, check_rmt = FALSE,
    check_icc = TRUE, common_features = ids))
  expect_equal(same$icc, 1)
  expect_identical(same$icc_status, "evaluated_experimental")
  fit$V_list[[2]] <- fit$V_list[[2]][c(5, 3, 1, 4, 2), , drop = FALSE]
  permuted <- suppressMessages(significant_components(fit, n = 20, check_rmt = FALSE,
    check_icc = TRUE, common_features = ids))
  expect_equal(permuted$icc, same$icc)
  expect_equal(permuted$icc_pvalue, same$icc_pvalue)
  fit <- .common_loading_fixture()
  fit$V_list[[2]][, 1] <- c(0, 1, 0, 0, 0)
  orthogonal <- suppressMessages(significant_components(fit, n = 20,
    check_rmt = FALSE, check_icc = TRUE, common_features = ids))
  expect_equal(orthogonal$icc, 0)
})

test_that("invalid common spaces cannot be coerced into coherence results", {
  fit <- .common_loading_fixture()
  ids <- rownames(fit$V_list[[1]])
  for (bad in list(rep("same", 5), c(ids[-1], ""), ids[-1])) {
    expect_error(significant_components(fit, n = 20, check_icc = TRUE,
      common_features = bad), "common")
  }
  fit$V_list[[2]] <- unname(fit$V_list[[2]])
  expect_error(significant_components(fit, n = 20, check_icc = TRUE,
    common_features = ids), "common feature domain")
  fit <- .common_loading_fixture()
  fit$V_list[[2]][, 1] <- 0
  expect_error(significant_components(fit, n = 20, check_icc = TRUE,
    common_features = ids), "zero.norm")
  fit$V_list[[2]][1, 1] <- Inf
  expect_error(significant_components(fit, n = 20, check_icc = TRUE,
    common_features = ids), "finite")
})

test_that("component retention uses feature rows independently of retained component count", {
  V <- diag(20)[, 1:3]
  fit <- list(V_list = list(V, V), sdev = sqrt(c(9, 1, 0.5)), s = matrix(0, 10, 3))
  result <- estimate_components(fit, noise_variance = 0.1,
    spectrum = "covariance", noise_model = "isotropic")
  expect_equal(result$diagnostics$k_vec, c(20, 20))
  expect_equal(result$diagnostics$gamma, 4)
  expect_equal(result$diagnostics$mp_edge, 0.9, tolerance = 1e-12)
  expect_equal(result$criterion, c(TRUE, TRUE, FALSE))
  expect_equal(result$keep, 1:2)
  fit$V_list <- lapply(fit$V_list, function(V) V[, 1:2])
  fit$sdev <- fit$sdev[1:2]
  fewer <- estimate_components(fit, noise_variance = 0.1,
    spectrum = "covariance", noise_model = "isotropic")
  expect_equal(fewer$diagnostics$gamma, result$diagnostics$gamma)
})

test_that("heterogeneous feature dimensions are retained without implicit ICC", {
  fit <- list(V_list = lapply(c(3, 6, 9), function(p) matrix(seq_len(p), p, 1)),
               sdev = 3, s = matrix(0, 100, 1))
  result <- estimate_components(fit, noise_variance = 0.1,
    spectrum = "covariance", noise_model = "isotropic")
  expect_equal(result$diagnostics$k_vec, c(3, 6, 9))
  expect_identical(result$diagnostics$icc_status, "disabled")
  fit$block_indices <- list(1:4, 5:10, 11:19)
  expect_error(estimate_components(fit), "block.*indices.*feature counts")
})

test_that("missing RMT information never silently selects the variance method", {
  fit <- .common_loading_fixture()
  result <- estimate_components(fit)
  expect_identical(result$method_used, "rmt_unavailable")
  expect_true(all(is.na(result$criterion)))
  expect_length(result$keep, 0)
})

test_that("retention wrapper passes feature rows across the criterion boundary", {
  captured <- NULL
  testthat::local_mocked_bindings(significant_components = function(fit, n, k_vec, ...) {
    captured <<- k_vec
    list(keep = integer(0), rmt_pass = rep(FALSE, 3), rmt_status = "evaluated")
  }, .package = "muscal")
  V <- matrix(1, 20, 3)
  estimate_components(list(V_list = list(V, V), sdev = c(3, 2, 1), s = matrix(0, 10, 3)))
  expect_equal(unname(captured), c(20, 20))
})

test_that("retention alpha controls explicitly requested common-feature ICC", {
  fit <- .common_loading_fixture()
  fit$s <- matrix(0, 20, 1)
  ids <- rownames(fit$V_list[[1]])
  run <- function(alpha) estimate_components(fit, alpha = alpha,
    check_icc = TRUE, common_features = ids, noise_variance = 0.1,
    spectrum = "covariance", noise_model = "isotropic")
  strict <- run(0.001)
  loose <- run(0.01)
  expect_equal(strict$diagnostics$alpha, 0.001)
  expect_equal(loose$diagnostics$alpha, 0.01)
  expect_identical(strict$criterion, FALSE)
  expect_identical(loose$criterion, TRUE)
  expect_length(strict$keep, 0)
  expect_equal(loose$keep, 1L)
  direct <- significant_components(fit, n = 20, alpha = 0.001,
    check_icc = TRUE, common_features = ids, noise_variance = 0.1,
    spectrum = "covariance", noise_model = "isotropic")
  expect_equal(strict$keep, direct$keep)
  expect_equal(strict$criterion, direct$rmt_pass & direct$icc_pass)
})

test_that("invalid or irrelevant retention alpha cannot silently govern an edge", {
  fit <- .common_loading_fixture()
  fit$s <- matrix(0, 20, 1)
  for (alpha in list(0, 1, -0.1, NA_real_, Inf, c(0.1, 0.2))) {
    expect_error(estimate_components(fit, alpha = alpha), "alpha")
    expect_error(significant_components(fit, n = 20, alpha = alpha), "alpha")
  }
  expect_warning(estimate_components(fit, alpha = 0.1), "alpha.*only.*ICC")
  expect_warning(estimate_components(fit, method = "variance", alpha = 0.1),
                  "alpha.*only.*ICC")
  expect_warning(a <- significant_components(fit, n = 20, alpha = 0.1,
    noise_variance = 0.1, spectrum = "covariance", noise_model = "isotropic"),
    "alpha.*only.*ICC")
  b <- significant_components(fit, n = 20, noise_variance = 0.1,
    spectrum = "covariance", noise_model = "isotropic")
  expect_equal(a$mp_edge, b$mp_edge)
})
