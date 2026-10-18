test_that("one-axis COVSTATIS preserves rank-one geometry across units and normalization", {
  x <- c(-2, -1, 1, 2)
  for (scale in c(0.04, 4.5)) {
    B <- scale * tcrossprod(x)
    for (method in c("none", "frobenius", "mfa")) {
      expected <- if (method == "none") B else B / (scale * sum(x^2))
      fit <- covstatis(list(B, B), ncomp = 1, dcenter = FALSE,
                       norm_method = method)
      expect_equal(dim(fit$roi_scores), c(4L, 1L))
      expect_equal(dim(fit$projmat), c(4L, 1L))
      expect_true(all(is.finite(fit$roi_scores)))
      expect_true(all(is.finite(fit$projmat)))
      expect_equal(unname(tcrossprod(fit$roi_scores)), expected, tolerance = 1e-12)
      expect_equal(unname(multivarious::reconstruct(fit)), expected, tolerance = 1e-12)
      barycenter <- Reduce(`+`, Map(`*`, fit$partial_scores, fit$alpha))
      expect_equal(unname(barycenter), unname(fit$roi_scores), tolerance = 1e-12)
      expect_equal(as.numeric(crossprod(fit$roi_scores)), sum(diag(expected)),
                   tolerance = 1e-12)
    }
  }
})

test_that("COVSTATIS scalar-axis repair preserves multi-axis and zero-rank contracts", {
  x <- c(-2, -1, 1, 2)
  z <- c(1, -1, -1, 1)
  B <- tcrossprod(x) + 0.3 * tcrossprod(z)
  fit <- covstatis(list(B, B), ncomp = 2, dcenter = FALSE, norm_method = "none")
  expect_equal(dim(fit$roi_scores), c(4L, 2L))
  expect_equal(unname(multivarious::reconstruct(fit)), B, tolerance = 1e-12)
  expect_error(covstatis(list(B * 0, B * 0), ncomp = 1, dcenter = FALSE,
                         norm_method = "none"), "no positive components")
})
