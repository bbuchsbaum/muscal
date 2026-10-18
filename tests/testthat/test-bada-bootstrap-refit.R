.bada_boot_fixture <- function(angles = c(0, 0.2, 0.4), amplitudes = rep(1, 3)) {
  U <- cbind(c(1, -1, 0) / sqrt(2), c(1, 1, -2) / sqrt(6))
  blocks <- lapply(seq_along(angles), function(s) {
    angle <- angles[[s]]
    R <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2)
    X <- amplitudes[[s]] * U %*% diag(c(3, 1)) %*% R
    X <- sweep(X, 2, c(10 * s, -5 * s), "+")
    colnames(X) <- c("f1", "f2")
    X
  })
  design <- data.frame(group = factor(rep(c("a", "b", "c"), length(blocks))),
                       subject = factor(rep(paste0("s", seq_along(blocks)), each = 3)))
  md <- multidesign::multidesign(do.call(rbind, blocks), design)
  fit <- bada(md, y = group, subject = subject, ncomp = 2)
  list(data = md, fit = fit, blocks = blocks)
}

test_that("identical-subject bootstrap refits preserve loading units", {
  f <- .bada_boot_fixture(angles = c(0, 0, 0))
  f$fit$ncomp <- 2L # Isolate the old scale error from its missing-count metadata bug.
  result <- multivarious::bootstrap(f$fit, data = f$data, nboot = 2)
  expect_equal(result$boot_lds_mean, unname(f$fit$v), tolerance = 1e-10)
  expect_equal(result$boot_lds_sd, matrix(0, 2, 2), tolerance = 1e-10)
})

test_that("BaDA bootstrap refits match independent class-mean SVD oracles", {
  f <- .bada_boot_fixture()
  schedule <- list(c(1, 1, 2), c(2, 2, 3), c(3, 3, 1), c(1, 2, 3))
  result <- multivarious::bootstrap(f$fit, data = f$data, nboot = 4,
    resample_indices = schedule, keep_replicates = TRUE, seed = 117)
  expect_equal(result$n_success, 4L)
  expect_length(result$failures, 0)
  expect_equal(result$resample_indices, schedule)
  expected_bases <- expected_scores <- vector("list", 4)
  for (i in seq_along(schedule)) {
    # Independent means/SVD and a two-permutation oracle; no production refit/alignment helper.
    centered <- lapply(f$blocks[schedule[[i]]], function(X) sweep(X, 2, colMeans(X), "-"))
    B <- Reduce(`+`, centered) / length(centered)
    V <- svd(B, nu = 0, nv = 2)$v
    choices <- list(1:2, 2:1)
    quality <- vapply(choices, function(p) sum(abs(colSums(f$fit$v * V[, p]))), numeric(1))
    V <- V[, choices[[which.max(quality)]], drop = FALSE]
    signs <- sign(colSums(f$fit$v * V))
    V <- sweep(V, 2, signs, "*")
    expected_bases[[i]] <- V
    expected_scores[[i]] <- B %*% V
    expect_equal(unname(result$replicates[[i]]$loadings), V, tolerance = 1e-10)
    expect_equal(unname(result$replicates[[i]]$scores), unname(B %*% V), tolerance = 1e-10)
    expect_equal(result$replicates[[i]]$n_subjects, 3L)
  }
  expect_equal(result$boot_lds_mean, Reduce(`+`, expected_bases) / 4, tolerance = 1e-10)
  expect_equal(result$boot_scores_mean, unname(Reduce(`+`, expected_scores) / 4), tolerance = 1e-10)
  expect_gt(max(abs(expected_bases[[1]] - unname(f$fit$v))), 1e-3)
})

test_that("BaDA bootstrap validates targets and records degenerate refits", {
  f <- .bada_boot_fixture(angles = c(0, 0, 0), amplitudes = c(1, 2, -2))
  schedule <- list(c(1, 1, 3), c(1, 1, 1), c(2, 2, 2))
  result <- multivarious::bootstrap(f$fit, data = f$data, nboot = 3,
    resample_indices = schedule)
  expect_equal(result$n_success, 2L)
  expect_length(result$failures, 1)
  expect_equal(result$failures[[1]]$index, 1L)
  expect_error(multivarious::bootstrap(f$fit, data = f$data, nboot = 2,
    resample_indices = rep(schedule[1], 2)), "valid|failed")
  unsupported <- f$fit
  unsupported$rescomp <- 1L
  expect_error(multivarious::bootstrap(unsupported, data = f$data, nboot = 2), "rescomp = 0")
  expect_error(multivarious::bootstrap(f$fit, data = f$data, nboot = 0), "nboot")
})

test_that("BaDA seeded bootstrap and subject identities are reproducible", {
  f <- .bada_boot_fixture()
  first <- multivarious::bootstrap(f$fit, data = f$data, nboot = 4, seed = 118)
  second <- multivarious::bootstrap(f$fit, data = f$data, nboot = 4, seed = 118)
  expect_equal(first$resample_indices, second$resample_indices)
  expect_equal(first$boot_lds_mean, second$boot_lds_mean)
  bad <- f$data
  bad$design$group[1] <- "b"
  expect_error(multivarious::bootstrap(f$fit, data = bad, nboot = 2), "fitted")
})

test_that("BaDA scalar and tied-root refits report appropriate coordinate summaries", {
  f <- .bada_boot_fixture(angles = c(0, 0, 0))
  scalar <- bada(f$data, y = group, subject = subject, ncomp = 1)
  out <- multivarious::bootstrap(scalar, data = f$data, nboot = 2, seed = 7)
  expect_equal(dim(out$boot_lds_mean), c(2L, 1L))
  expect_equal(out$boot_lds_mean, unname(scalar$v), tolerance = 1e-10)
  f$data$x[, 2] <- f$data$x[, 2] * 3
  tied <- bada(f$data, y = group, subject = subject, ncomp = 2)
  out <- multivarious::bootstrap(tied, data = f$data, nboot = 2, seed = 7)
  expect_true(all(is.na(out$boot_lds_mean)))
  expect_true(all(is.na(out$boot_scores_mean)))
  expect_equal(out$alignment$identifiability, "ambiguous_reference_roots")
  expect_equal(as.vector(out$alignment$subspace_distance), c(0, 0), tolerance = 1e-10)
})
