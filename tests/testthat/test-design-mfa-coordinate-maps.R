test_that("exact coordinate matches produce an identity sparse map", {
  reference <- cbind(x = 1:6, y = 0, z = 0)
  rownames(reference) <- paste0("v", 1:6)
  out <- coordinate_feature_map(reference, reference, k = 1, kernel = "nearest")

  expect_s4_class(out, "sparseMatrix")
  expect_equal(unname(as.matrix(out)), diag(6), tolerance = 0)
  expect_identical(rownames(out), rownames(reference))
  expect_identical(colnames(out), rownames(reference))
  expect_lte(length(out@x), 6L)
})

test_that("Gaussian coordinate maps agree with a dense distance oracle", {
  reference <- cbind(
    x = c(-2, -1, 0, 1, 2),
    y = c(0, 1, 0, -1, 0),
    z = c(0, 0, 1, 0, -1)
  )
  subject <- rbind(c(-1.4, 0.4, 0.1), c(0.2, -0.3, 0.7), c(1.6, -0.2, -0.4))
  rownames(reference) <- paste0("r", 1:5)
  rownames(subject) <- paste0("s", 1:3)
  bandwidth <- 1.7

  out <- coordinate_feature_map(
    reference,
    subject,
    k = 3,
    bandwidth = bandwidth,
    kernel = "gaussian"
  )

  distance <- as.matrix(stats::dist(rbind(subject, reference)))[1:3, 4:8]
  oracle <- matrix(0, 3, 5)
  for (i in 1:3) {
    keep <- order(distance[i, ], seq_len(5))[1:3]
    weight <- exp(-0.5 * (distance[i, keep] / bandwidth)^2)
    oracle[i, keep] <- weight / sum(weight)
  }
  dimnames(oracle) <- list(rownames(subject), rownames(reference))

  expect_equal(as.matrix(out), oracle, tolerance = 1e-12)
  expect_equal(unname(Matrix::rowSums(out)), rep(1, 3), tolerance = 1e-14)
})

test_that("coordinate maps are equivariant to rigid transforms and storage order", {
  set.seed(70)
  reference <- matrix(rnorm(24), 8, 3, dimnames = list(paste0("r", 1:8), c("x", "y", "z")))
  subject <- matrix(rnorm(15), 5, 3, dimnames = list(paste0("s", 1:5), c("x", "y", "z")))
  rotation <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  translation <- matrix(c(12, -4, 9), 1, 3)

  original <- coordinate_feature_map(reference, subject, k = 4, bandwidth = 0.8)
  transformed <- coordinate_feature_map(
    sweep(reference %*% rotation, 2, translation, "+"),
    sweep(subject %*% rotation, 2, translation, "+"),
    k = 4,
    bandwidth = 0.8
  )
  reference_order <- c(5, 1, 8, 3, 2, 7, 4, 6)
  subject_order <- c(4, 2, 5, 1, 3)
  permuted <- coordinate_feature_map(
    reference[reference_order, , drop = FALSE],
    subject[subject_order, , drop = FALSE],
    k = 4,
    bandwidth = 0.8
  )

  expect_equal(as.matrix(transformed), as.matrix(original), tolerance = 1e-12)
  expect_equal(
    as.matrix(permuted)[rownames(subject), rownames(reference)],
    as.matrix(original),
    tolerance = 1e-12
  )
})

test_that("coordinate-derived maps restore permuted subject fields", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 28, p = 9, seed = 71)
  reference <- cbind(
    x = seq(-20, 20, length.out = 9),
    y = sin(seq(0, pi, length.out = 9)),
    z = cos(seq(0, pi, length.out = 9))
  )
  rownames(reference) <- sim$feature_names
  subject_coords <- Y_permuted <- vector("list", 3)
  set.seed(72)
  for (s in seq_along(sim$Y)) {
    permutation <- sample.int(9)
    observed_names <- paste0("S", s, "_voxel_", seq_len(9))
    subject_coords[[s]] <- reference[permutation, , drop = FALSE]
    rownames(subject_coords[[s]]) <- observed_names
    Y_permuted[[s]] <- sim$Y[[s]][, permutation, drop = FALSE]
    colnames(Y_permuted[[s]]) <- observed_names
  }
  names(subject_coords) <- names(Y_permuted) <- names(sim$Y)
  maps <- coordinate_feature_maps(
    reference,
    subject_coords,
    k = 1,
    kernel = "nearest"
  )

  shared <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    rank = 2
  )
  soft <- design_mfa(
    Y_permuted,
    sim$design,
    ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    alignment = "soft",
    feature_map = maps,
    rank = 2,
    control = design_mfa_control(soft_tol = 1e-10)
  )

  expect_equal(predict(soft), predict(shared), tolerance = 1e-8)
  expect_equal(
    .dmfa_projector(effect_loadings(soft, "saliency")),
    .dmfa_projector(effect_loadings(shared, "saliency")),
    tolerance = 1e-7
  )
  expect_equal(alignment_coverage(soft)$coverage_fraction, rep(1, 9))
})

test_that("coordinate map validation exposes unsupported geometry", {
  reference <- cbind(x = 1:5, y = 0, z = 0)
  subject <- cbind(x = c(1, 2), y = 0, z = 0)
  rownames(reference) <- paste0("r", 1:5)
  rownames(subject) <- paste0("s", 1:2)

  duplicate <- reference
  duplicate[2, ] <- duplicate[1, ]
  expect_error(coordinate_feature_map(duplicate, subject), "duplicate coordinates")
  contaminated <- subject
  contaminated[1, 1] <- Inf
  expect_error(coordinate_feature_map(reference, contaminated), "finite numeric")
  expect_error(
    coordinate_feature_map(reference, cbind(subject, extra = 1)),
    "same coordinate dimensions"
  )
  expect_error(
    coordinate_feature_map(reference, subject + 100, max_distance = 0.1),
    "no reference neighbor"
  )
  expect_error(
    coordinate_feature_maps(reference, list(S1 = subject), k = 1),
    "not covered"
  )
  expect_error(coordinate_feature_map(reference, subject, k = 1.5), "positive integer")
  expect_error(coordinate_feature_map(reference, subject, normalize_rows = 1), "TRUE or FALSE")
  expect_error(coordinate_feature_map(reference, subject, bandwidth = "near"), "positive finite")
  expect_error(
    coordinate_feature_maps(reference, list(S1 = subject), require_full_reference = 1, k = 1),
    "TRUE or FALSE"
  )
})

test_that("coordinate mapping remains sparse at representative scale", {
  set.seed(73)
  reference <- matrix(rnorm(5000 * 3), 5000, 3)
  subject <- matrix(rnorm(4000 * 3), 4000, 3)
  rownames(reference) <- paste0("r", seq_len(nrow(reference)))
  rownames(subject) <- paste0("s", seq_len(nrow(subject)))
  out <- coordinate_feature_map(reference, subject, k = 4, bandwidth = 0.3)

  expect_lte(length(out@x), 4L * nrow(subject))
  expect_lt(as.numeric(object.size(out)), nrow(subject) * nrow(reference))
  expect_true(all(Matrix::rowSums(out) > 0))
})
