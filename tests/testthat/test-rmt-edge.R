.rmt_fixture <- function(roots = c(4, 4.5, 5)) {
  V <- diag(8)[, seq_along(roots), drop = FALSE]
  list(V_list = list(V, V), sdev = sqrt(roots))
}

test_that("automatic RMT is unavailable rather than a disguised median threshold", {
  fit <- .rmt_fixture()
  result <- significant_components(fit, n = 64)
  expect_identical(result$rmt_status, "unavailable")
  expect_true(all(is.na(result$rmt_pass)))
  expect_true(is.na(result$mp_edge))
  expect_length(result$keep, 0)
  expect_match(result$rmt_reason, "noise variance.*spectral units")
})

test_that("declared isotropic covariance edges have independent analytic values", {
  fit <- .rmt_fixture()
  small <- significant_components(fit, n = 64, noise_variance = 2,
    spectrum = "covariance", noise_model = "isotropic", n_eff = 64)
  square <- significant_components(fit, n = 64, noise_variance = 2,
    spectrum = "covariance", noise_model = "isotropic", n_eff = 16)
  wide <- significant_components(fit, n = 64, noise_variance = 2,
    spectrum = "covariance", noise_model = "isotropic", n_eff = 4)
  expect_equal(small$mp_edge, 4.5, tolerance = 1e-12)
  expect_equal(square$mp_edge, 8, tolerance = 1e-12)
  expect_equal(wide$mp_edge, 18, tolerance = 1e-12)
  expect_identical(small$rmt_pass, c(FALSE, FALSE, TRUE))
  expect_identical(small$rmt_status, "evaluated")
  expect_equal(small$gamma, 0.25)
})

test_that("RMT unit conversion and scale transformations preserve decisions", {
  fit <- .rmt_fixture(c(4, 4.6, 5))
  cov <- significant_components(fit, n = 64, noise_variance = 2,
    spectrum = "covariance", noise_model = "isotropic")
  raw <- fit
  raw$sdev <- raw$sdev * 8
  sing <- significant_components(raw, n = 64, noise_variance = 2,
    spectrum = "singular", noise_model = "isotropic")
  expect_equal(sing$lambda, cov$lambda)
  expect_equal(sing$rmt_pass, cov$rmt_pass)
  fit$sdev <- fit$sdev * 3
  scaled <- significant_components(fit, n = 64, noise_variance = 18,
    spectrum = "covariance", noise_model = "isotropic")
  expect_equal(scaled$mp_edge, 9 * cov$mp_edge)
  expect_equal(scaled$rmt_pass, cov$rmt_pass)
})

test_that("RMT does not certify unsupported representations or invalid inputs", {
  fit <- .rmt_fixture()
  unknown <- significant_components(fit, n = 64, noise_variance = 2)
  expect_identical(unknown$rmt_status, "unavailable")
  expect_true(all(is.na(unknown$rmt_pass)))
  for (variance in list(0, -1, NA_real_, Inf, c(1, 2))) {
    expect_error(significant_components(fit, n = 64, noise_variance = variance,
      spectrum = "covariance", noise_model = "isotropic"), "noise_variance")
  }
  expect_error(significant_components(fit, n = 64, noise_variance = 2,
    spectrum = "covariance", noise_model = "isotropic", n_eff = 0), "n_eff")
  fit$sdev[1] <- Inf
  expect_error(significant_components(fit, n = 64, noise_variance = 2,
    spectrum = "covariance", noise_model = "isotropic"), "finite")
})
