# Published data/targets from Abdi et al. (2007), supplied manuscript pp. 8-10:
# Table 1, eqs. 21-22, and Table 2. Tolerances reflect printed rounding.
distatis_beer_fixture <- function() {
  sorting <- matrix(c(
    1,4,3,4,1,1,2,2,1,3,
    4,5,2,5,2,3,1,1,4,3,
    3,1,2,3,2,4,3,1,1,2,
    4,2,3,3,1,1,1,2,1,4,
    1,5,3,5,2,1,1,2,1,3,
    2,3,1,1,3,5,4,4,3,1,
    1,4,3,4,1,1,2,2,2,4,
    5,2,4,2,4,2,5,3,4,5
  ), 8, 10, byrow = TRUE, dimnames = list(
    c("Affligen", "Budweiser", "Buckler Blonde", "Killian", "St. Landelin",
      "Buckler Highland", "Fruit Defendu", "EKU28"), paste0("A", 1:10)))
  RV <- matrix(c(
    1,.57,.60,.57,.35,.68,.67,.60,.37,.45,
    .57,1,.42,.83,.65,.49,.76,.42,.28,.57,
    .60,.42,1,.40,.63,.93,.46,1,.44,.46,
    .57,.83,.40,1,.60,.49,.76,.40,.42,.57,
    .35,.65,.63,.60,1,.59,.53,.63,.33,.53,
    .68,.49,.93,.49,.59,1,.56,.93,.56,.56,
    .67,.76,.46,.76,.53,.56,1,.46,.35,.45,
    .60,.42,1,.40,.63,.93,.46,1,.44,.46,
    .37,.28,.44,.42,.33,.56,.35,.44,1,.35,
    .45,.57,.46,.57,.53,.56,.45,.46,.35,1
  ), 10, 10, byrow = TRUE)
  scores <- t(matrix(c(
    -.39,.23,.28,-.16,-.15,.30,-.40,.29,
    .07,-.37,-.22,-.03,-.24,.38,.15,.26,
    -.03,-.02,-.07,.09,-.04,-.40,-.01,.48,
    .03,-.25,.45,.08,-.22,-.10,.07,-.06,
    .16,.06,.07,-.42,-.01,-.05,.11,.06,
    .02,-.20,.05,-.05,.26,.01,-.13,.04,
    .16,.02,-.01,.03,-.06,.00,-.14,-.00
  ), 7, 8, byrow = TRUE))
  list(sorting = sorting, RV = RV, scores = scores,
       alpha = c(.100,.101,.109,.101,.099,.116,.101,.109,.074,.090),
       roots = c(.66,.49,.40,.34,.23,.13,.05))
}

test_that("DISTATIS reproduces the published beer example", {
  beer <- distatis_beer_fixture()
  fit <- distatis(sorting_distances(beer$sorting), ncomp = 10)
  expect_s3_class(fit, c("distatis", "covstatis", "projector"))
  expect_lt(max(abs(fit$C - beer$RV)), .0051)
  expect_lt(max(abs(fit$alpha - beer$alpha)), .00051)
  expect_lt(max(abs(fit$sdev^2 - beer$roots)), .0051)
  aligned <- sweep(fit$object_scores, 2,
                   sign(colSums(fit$object_scores * beer$scores)), `*`)
  expect_lt(max(abs(aligned - beer$scores)), .0051)
  expect_equal(fit$compromise_quality, .5982768, tolerance = 1e-7)
  expect_equal(sum(fit$inertia[1:4]), .8228512, tolerance = 1e-7)
  expect_equal(unname(diag(fit$C)), rep(1, 10), tolerance = 1e-12)
  expect_true(all(fit$alpha >= 0))
  expect_equal(sum(fit$alpha), 1, tolerance = 1e-12)
  expect_equal(multivarious::scores(fit), fit$object_scores)
  expect_equal(multivarious::sdev(fit), fit$sdev)
  expect_equal(multivarious::ncomp(fit), 7L)
  expect_equal(unname(crossprod(coef(fit))), diag(7), tolerance = 1e-12)
  expect_equal(Reduce(`+`, Map(`*`, fit$partial_scores, fit$alpha)),
               fit$object_scores, tolerance = 1e-12)
  expect_equal(tcrossprod(fit$subject_scores), fit$C, tolerance = 1e-12)
  expect_equal(unname(colSums(fit$subject_scores^2)),
               fit$interstructure_eigenvalues[fit$subject_axes], tolerance = 1e-12)
  expect_equal(multivarious::reconstruct(fit), fit$compromise_full, tolerance = 1e-12)
  expect_output(print(fit), "8 objects, 10 assessors")
})

test_that("nonbinary distances reduce to classical MDS and preserve squared semantics", {
  X <- cbind(c(0, 1, 3, 5), c(0, 2, 1, 4), c(0, 0, 2, 1))
  rownames(X) <- letters[1:4]
  D <- as.matrix(dist(X))
  B <- tcrossprod(scale(X, scale = FALSE))
  fit <- distatis(dist(X), ncomp = 4, norm_method = "none")
  expect_equal(fit$compromise_full, B, tolerance = 1e-12)
  expect_equal(multivarious::reconstruct(fit, what = "distance"), D, tolerance = 1e-12)
  expect_equal(multivarious::reconstruct(fit, what = "distance", squared = TRUE), D^2, tolerance = 1e-12)
  squared <- distatis(D^2, ncomp = 4, norm_method = "none", squared = TRUE)
  expect_equal(squared$compromise_full, B, tolerance = 1e-12)
  expect_equal(project_distance(fit, D^2, squared = TRUE), fit$partial_scores[[1]], tolerance = 1e-12)
  expect_equal(project_distance(squared, D^2), squared$partial_scores[[1]], tolerance = 1e-12)
  expect_equal(dim(fit$subject_scores), c(1L, 1L))
  expect_equal(dim(fit$object_scores), c(4L, 3L))
  # A coordinate-derived oracle checks the leading-axis, single-component path.
  sv <- svd(scale(X, scale = FALSE))
  one <- distatis(D, ncomp = 1, norm_method = "none")
  expect_equal(unname(multivarious::reconstruct(one)),
               tcrossprod(sv$u[, 1] * sv$d[1]), tolerance = 1e-12)
  expect_equal(dim(one$object_scores), c(4L, 1L))
  expect_equal(one$normalization_constants, c(Assessor_1 = 1))
})

test_that("lists, arrays, and scale or order transformations agree", {
  Ds <- sorting_distances(distatis_beer_fixture()$sorting)
  cube <- array(unlist(Ds, use.names = FALSE), c(8, 8, 10),
                dimnames = c(dimnames(Ds[[1]]), list(names(Ds))))
  fit <- distatis(Ds)
  expect_equal(distatis(cube)$compromise_full, fit$compromise_full, tolerance = 1e-12)
  expect_equal(distatis(lapply(Ds, as.dist))$compromise_full, fit$compromise_full, tolerance = 1e-12)
  for (method in c("mfa", "frobenius")) {
    base <- distatis(Ds, norm_method = method)
    scaled <- distatis(Map(`*`, Ds, seq_along(Ds) / 5), norm_method = method)
    expect_equal(scaled$C, base$C, tolerance = 1e-12)
    expect_equal(scaled$compromise_full, base$compromise_full, tolerance = 1e-12)
    tiny <- distatis(lapply(Ds, `*`, 1e-60), norm_method = method)
    expect_equal(tiny$compromise_full, base$compromise_full, tolerance = 1e-12)
  }
  scaled_none <- distatis(lapply(Ds, `*`, 3), norm_method = "none")
  none <- distatis(Ds, norm_method = "none")
  expect_equal(scaled_none$compromise_full, none$compromise_full * 9, tolerance = 1e-12)
  perm <- c(8, 3, 1, 7, 2, 6, 4, 5)
  objects <- distatis(lapply(Ds, function(D) D[perm, perm]))
  expect_equal(objects$compromise_full, fit$compromise_full[perm, perm], tolerance = 1e-12)
  tables <- distatis(rev(Ds))
  expect_equal(unname(tables$alpha), rev(unname(fit$alpha)), tolerance = 1e-12)
  expect_equal(tables$compromise_full, fit$compromise_full, tolerance = 1e-12)
})

test_that("supplementary training tables reproduce both maps with independent ranks", {
  Ds <- sorting_distances(distatis_beer_fixture()$sorting)
  fit <- distatis(Ds, ncomp = 10)
  projected <- project_subjects(fit, Ds)
  expect_equal(projected$subject_scores, fit$subject_scores, tolerance = 1e-12)
  expect_equal(projected$object_scores, fit$partial_scores, tolerance = 1e-12)
  expect_equal(project_distance(fit, Ds[[1]]), fit$partial_scores[[1]], tolerance = 1e-12)
  expect_equal(project(fit, Ds[[1]]), project_distance(fit, Ds[[1]]))
  expect_equal(unname(project_subjects(fit, Ds[[1]], subject_ids = "new")$subject_scores),
               unname(projected$subject_scores[1, , drop = FALSE]), tolerance = 1e-12)
  # Exact equality to a full-compromise table is easier to test with K=1.
  single <- distatis(Ds[[1]], ncomp = 8)
  result <- project_subjects(single, Ds[[1]])
  expect_equal(result$scalar_summaries$rv_coefficient, 1, tolerance = 1e-12)
  expect_lt(result$scalar_summaries$distance_to_compromise, 1e-12)
  expect_lt(result$scalar_summaries$distance_to_subspace, 1e-12)
  expect_equal(result$subject_scores, single$subject_scores, tolerance = 1e-12)
  expect_equal(dim(result$object_scores[[1]]), dim(single$object_scores))
  thin <- distatis(Ds, ncomp = 2)
  H <- diag(8) - matrix(1/8, 8, 8)
  B <- -0.5 * H %*% Ds[[1]]^2 %*% H
  B <- B / eigen(B, symmetric = TRUE, only.values = TRUE)$values[1]
  diagnostic <- project_subjects(thin, Ds[[1]])$scalar_summaries
  expect_equal(diagnostic$rv_coefficient,
               sum(B * thin$compromise_full) / sqrt(sum(B^2) * sum(thin$compromise_full^2)),
               tolerance = 1e-12)
  expect_equal(diagnostic$distance_to_compromise,
               sqrt(sum((B - thin$compromise_full)^2)), tolerance = 1e-12)
  expect_gt(abs(diagnostic$distance_to_compromise - diagnostic$distance_to_subspace), 1e-4)
  squared <- distatis(lapply(Ds, `^`, 2), squared = TRUE)
  expect_equal(project_distance(squared, Ds[[1]]^2), squared$partial_scores[[1]], tolerance = 1e-12)
  expect_equal(project_distance(squared, Ds[[1]], squared = FALSE), squared$partial_scores[[1]], tolerance = 1e-12)
  expect_error(project_cov(fit, Ds[[1]]), "project_distance")
  expect_error(infer_covstatis(fit, Ds), "DISTATIS inference")
  expect_error(as_multifer_adapter(fit), "DISTATIS inference")
  expect_error(check_duality(fit), "separate spectra")
})

test_that("sorting group names are categorical and degenerate sorts are explicit", {
  p <- data.frame(A = c("x", "x", "y", "z"), B = factor(c(1, 1, 2, 2)))
  Ds <- sorting_distances(p)
  renamed <- data.frame(A = c(90, 90, 7, -3), B = c("a", "a", "b", "b"))
  expect_equal(Ds, sorting_distances(renamed))
  expect_equal(Ds[[1]][1, 2], 0)
  expect_equal(Ds[[1]][1, 3], 1)
  expect_equal(unname(diag(Ds[[1]])), rep(0, 4))
  close_ids <- sorting_distances(data.frame(A = c(1, 1 + 1e-15, 2)))
  expect_equal(close_ids[[1]][1, 2], 1)
  singleton <- distatis(sorting_distances(data.frame(A = 1:4)))
  expect_true(all(is.finite(singleton$object_scores)))
  expect_error(distatis(sorting_distances(data.frame(A = rep(1, 4)))), "zero Gram")
  expect_error(sorting_distances(data.frame(A = c(1, NA, 2))), "nonmissing")
  expect_error(sorting_distances(matrix(1, 1, 2)), "at least two")
  expect_error(sorting_distances(data.frame(A = I(list(1, 2, 3)))), "categorical")
  expect_equal(dim(distatis(list(Ds[[1]], Ds[[1]]), ncomp = 1)$object_scores), c(4L, 1L))
  # Orthogonal centered coordinate vectors give C = I, with no unique leading root.
  expect_error(distatis(list(as.matrix(dist(c(-1, -1, 1, 1))),
                            as.matrix(dist(c(-1, 1, -1, 1))))), "leading RV eigenvalue is tied")
})

test_that("DISTATIS rejects invalid geometry, input metadata, and options", {
  D <- as.matrix(dist(c(0, 1, 3)))
  for (bad in list(list(), data.frame(a = 1:3), matrix(1, 2, 3), matrix(0, 1, 1),
                   D * NA_real_, D * Inf, list(D, diag(4)))) {
    expect_error(distatis(bad))
  }
  asymmetric <- D; asymmetric[1, 2] <- 4
  negative <- D; negative[1, 2] <- negative[2, 1] <- -1
  diagonal <- D; diag(diagonal) <- 1
  noneuclidean <- matrix(c(0, 1, 1, 1, 0, 3, 1, 3, 0), 3)
  expect_error(distatis(asymmetric), "symmetric")
  expect_error(distatis(negative), "nonnegative")
  expect_error(distatis(diagonal), "zero diagonal")
  expect_error(distatis(noneuclidean), "non-Euclidean.*negative inertia")
  for (k in list(0, -1, 1.2, NA_real_, Inf, TRUE, numeric())) expect_error(distatis(D, ncomp = k), "positive integer")
  expect_error(distatis(D, squared = 1), "TRUE or FALSE")
  expect_error(distatis(D, tol = 0), "tol")
  expect_error(distatis(D, tol = NA_real_), "tol")
  expect_error(distatis(D, labels = c("x", "x", "y")), "unique")
  expect_error(distatis(D, typo = TRUE), "Unused arguments")
  expect_error(distatis(list(A = D, A = D)), "assessor names")
  named <- D; dimnames(named) <- list(letters[1:3], letters[1:3])
  reversed <- named[3:1, 3:1]
  expect_error(distatis(list(named, reversed)), "names/order")
  mismatch <- named; colnames(mismatch) <- rev(colnames(named))
  expect_error(distatis(mismatch), "row and column")
  fit <- distatis(named, labels = LETTERS[1:3], ncomp = 1)
  expect_equal(rownames(fit$object_scores), LETTERS[1:3])
  expect_equal(rownames(project_distance(fit, named)), LETTERS[1:3])
  expect_error(project_distance(fit, reversed), "names/order")
  expect_error(project_distance(fit, diag(4)), "dimensions")
  expect_error(project_subjects(fit, list(named), subject_ids = c("a", "b")), "subject_ids")
  expect_error(multivarious::reconstruct(fit, comp = c(1, 1)), "distinct")
  expect_error(multivarious::reconstruct(fit, comp = 0), "indices")
  expect_error(multivarious::reconstruct(fit, comp = numeric()), "indices")
})

test_that("covariates and plots respect separate spaces and full inertia", {
  Ds <- sorting_distances(distatis_beer_fixture()$sorting)
  fit <- distatis(Ds)
  y <- seq_along(Ds)
  expected <- drop(crossprod(y, fit$subject_scores)) /
    sqrt(colSums(fit$subject_scores^2) * sum(y^2))
  expect_equal(project_table_covariate(fit, y), expected)
  expect_equal(length(project_feature_covariate(fit, 1:8)), 2L)
  expect_equal(dim(project_table_covariate(fit, y, what = "feature")), c(8L, 2L))
  expect_equal(unname(project_table_covariate(fit, y * 0)), c(0, 0))
  expect_error(project_table_covariate(fit, c(NA, 1:9)), "finite")
  expect_error(project_feature_covariate(fit, 1:8, what = "feature"), "only support")
  object_plot <- ggplot2::autoplot(fit, labels = TRUE)
  subject_plot <- ggplot2::autoplot(fit, type = "subjects")
  expect_equal(object_plot$data$label, fit$labels)
  expect_equal(subject_plot$data$label, fit$block_labels)
  expect_equal(object_plot$coordinates$ratio, 1)
  expect_equal(subject_plot$coordinates$ratio, 1)
  expect_match(object_plot$labels$x, "28.8%", fixed = TRUE)
  expect_match(subject_plot$labels$x, "59.8%", fixed = TRUE)
  for (p in list(object_plot, subject_plot, plot_partial_scores(fit),
                 plot_block_weights(fit), plot_block_similarity(fit),
                 plot_variance(fit), plot_variance(fit, axes = "subjects", type = "line"))) {
    expect_s3_class(p, "ggplot")
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  }
  variance <- plot_variance(fit)
  expect_equal(variance$data$pct[1], 100 * fit$inertia[1], tolerance = 1e-12)
  single <- distatis(Ds[[1]], ncomp = 5)
  expect_equal(length(project_table_covariate(single, 1)), 1L)
  expect_equal(dim(project_table_covariate(single, 1, what = "feature", scale = "beta")),
               dim(single$object_scores))
  expect_error(ggplot2::autoplot(single, type = "subjects"), "map has 1")
  expect_error(ggplot2::autoplot(fit, dims = c(1, 1)), "distinct")
  expect_error(plot_partial_scores(fit, blocks = 99), "assessor indices")
  expect_error(plot_block_similarity(fit, method = "RV2"), "not RV2")
  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output)
  tryCatch({
    expect_equal(plot(fit, labels = TRUE), fit$object_scores)
    expect_equal(plot_subjects(fit), fit$subject_scores)
  }, finally = { grDevices::dev.off(); unlink(output) })
})

test_that("COVSTATIS one-component Gram fits retain the correct geometry", {
  X <- cbind(0:3, c(0, 2, 1, 5))
  B <- tcrossprod(scale(X, scale = FALSE))
  fit <- covstatis(list(B, B), ncomp = 1, dcenter = FALSE, norm_method = "none")
  ev <- eigen(B, symmetric = TRUE)
  expect_equal(unname(multivarious::reconstruct(fit)),
               tcrossprod(ev$vectors[, 1] * sqrt(ev$values[1])), tolerance = 1e-12)
})
