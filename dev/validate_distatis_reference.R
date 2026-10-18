# Optional development check; DistatisR is not a runtime dependency.
# From repository root, with DistatisR 1.1.2 installed:
# LC_ALL=C Rscript --vanilla dev/validate_distatis_reference.R
stopifnot(requireNamespace("DistatisR", quietly = TRUE),
          as.character(utils::packageVersion("DistatisR")) == "1.1.2")
pkgload::load_all(quiet = TRUE)

beer <- matrix(c(
  1,4,3,4,1,1,2,2,1,3,
  4,5,2,5,2,3,1,1,4,3,
  3,1,2,3,2,4,3,1,1,2,
  4,2,3,3,1,1,1,2,1,4,
  1,5,3,5,2,1,1,2,1,3,
  2,3,1,1,3,5,4,4,3,1,
  1,4,3,4,1,1,2,2,2,4,
  5,2,4,2,4,2,5,3,4,5
), 8, 10, byrow = TRUE)
beer_distances <- muscal::sorting_distances(beer)

# Nonbinary Euclidean distances test the squared-input convention independently
# of sorting, whose binary distances equal their own squares.
X <- cbind(c(0, 1, 3, 5, 2, -1), c(0, 2, 1, 4, -1, 1), c(0, 1, 2, 0, 4, 3))
general <- list(A = as.matrix(dist(X)),
                B = as.matrix(dist(sweep(X, 2, c(2, .5, 1.5), `*`))),
                C = as.matrix(dist(X[, c(2, 3, 1)] + diag(6)[, 1:3])))

compare_reference <- function(Ds, method, reference_method, name) {
  n <- nrow(Ds[[1]])
  cube <- array(unlist(lapply(Ds, `^`, 2), use.names = FALSE), c(n, n, length(Ds)),
                dimnames = list(rownames(Ds[[1]]), colnames(Ds[[1]]), names(Ds)))
  reference <- DistatisR::distatis(cube, Norm = reference_method, Distance = TRUE,
                                  RV = TRUE, nfact2keep = 3)
  fit <- muscal::distatis(Ds, ncomp = 3, norm_method = method)
  signs <- sign(colSums(fit$object_scores * reference$res4Splus$F))
  partial_error <- max(vapply(seq_along(Ds), function(k) {
    max(abs(sweep(fit$partial_scores[[k]], 2, signs, `*`) - reference$res4Splus$PartialF[, , k]))
  }, numeric(1)))
  errors <- c(
    RV = max(abs(fit$C - reference$res4Cmat$C)),
    alpha = max(abs(fit$alpha - reference$res4Cmat$alpha)),
    compromise = max(abs(fit$compromise_full - reference$res4Splus$Splus)),
    roots = max(abs(fit$eigenvalues - reference$res4Splus$eigValues)),
    object_scores = max(abs(sweep(fit$object_scores, 2, signs, `*`) - reference$res4Splus$F)),
    partial_scores = partial_error,
    assessor_gram = max(abs(tcrossprod(fit$subject_scores) -
                             tcrossprod(reference$res4Cmat$G[, fit$subject_axes, drop = FALSE])))
  )
  if (any(!is.finite(errors)) || max(errors) >= 1e-10) {
    print(data.frame(case = name, normalization = method, quantity = names(errors),
                     maximum_absolute_error = unname(errors)), row.names = FALSE)
    stop("Installed-reference comparison failed.")
  }
  data.frame(case = name, normalization = method, quantity = names(errors),
             maximum_absolute_error = unname(errors))
}

results <- do.call(rbind, lapply(names(list(beer = beer_distances, general = general)), function(name) {
  Ds <- if (name == "beer") beer_distances else general
  do.call(rbind, list(compare_reference(Ds, "mfa", "MFA", name),
                      compare_reference(Ds, "none", "NONE", name)))
}))
cat("DISTATIS installed-reference parity\n")
cat("R:", as.character(getRversion()), "DistatisR:", as.character(utils::packageVersion("DistatisR")), "\n")
cat("Reference input: squared distances; RV = TRUE; three retained axes.\n")
print(results, row.names = FALSE, digits = 5)
cat("All reference comparisons passed at absolute tolerance 1e-10.\n")
cat("Direct parity covers mfa and none. Reference SUMPCA is trace normalization, not Frobenius.\n")
