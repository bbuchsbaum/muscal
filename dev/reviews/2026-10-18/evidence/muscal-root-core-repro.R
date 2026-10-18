pkgload::load_all('.', compile = FALSE, quiet = TRUE)

expect_error_text <- function(expr, pattern, id) {
  msg <- tryCatch({ force(expr); "<no error>" }, error = conditionMessage)
  stopifnot(grepl(pattern, msg, fixed = TRUE))
  cat(id, ": ", msg, "\n", sep = "")
}

# C1: default k_nn = 6 exceeds the available pooled neighbours.
expect_error_text(
  penalized_mfa_clusterwise(
    list(matrix(rnorm(4), 2, 2), matrix(rnorm(4), 2, 2)),
    list(matrix(c(0,0,0, 1,0,0), 2, 3, byrow = TRUE),
         matrix(c(0,1,0, 1,1,0), 2, 3, byrow = TRUE)),
    ncomp = 1, max_iter = 1),
  "NA' indices", "C1")

# C2: unequal row counts are accepted but consensus-score addition fails.
expect_error_text(
  penalized_mfa(list(matrix(rnorm(12), 3, 4), matrix(rnorm(20), 5, 4)), ncomp = 1),
  "non-conformable arrays", "C2")

# C3: bootstrap code fixes the fitted loading basis for every replicate.
bada_lines <- readLines("R/bada.R")[84:138]
bada_src <- paste(bada_lines, collapse = "\n")
stopifnot(grepl("v_boot   <- x\\$v", bada_src),
          grepl("fs_boot  <- x\\$fscores", bada_src),
          !grepl("[^#]pca_group <-", bada_src, perl = TRUE))
cat("C3: bootstrap replicate uses fixed x$v and x$fscores; no active PCA refit\n")

# C7: normalizing an already formed Laplacian with rowSums(L + I) is a no-op.
A <- Matrix::Matrix(c(0,1,1, 1,0,0, 1,0,0), 3, 3, sparse = TRUE)
L <- Matrix::Diagonal(x = Matrix::rowSums(A)) - A
Dinv <- Matrix::Diagonal(x = 1 / sqrt(pmax(Matrix::rowSums(L + Matrix::Diagonal(3)), 1e-10)))
stopifnot(isTRUE(all.equal(as.matrix(Dinv %*% L %*% Dinv), as.matrix(L))))
cat("C7: documented normalized branch equals unnormalized Laplacian\n")

# C8: forceSymmetric with uplo='U' deletes lower-only directed edges.
Aone <- Matrix::Matrix(0, 3, 3, sparse = TRUE); Aone[3, 1] <- 1
stopifnot(Matrix::forceSymmetric(Aone, uplo = "U")[3, 1] == 0)
cat("C8: lower-only directed edge deleted by forceSymmetric(uplo='U')\n")

# C9: fit centers inputs but stores a fitted pass() preprocessor.
set.seed(2)
xs <- list(matrix(rnorm(20, mean = 10), 5, 4), matrix(rnorm(20, mean = -5), 5, 4))
cs <- list(cbind(1:4, 0, 0), cbind(1:4, 1, 0))
f <- penalized_mfa_clusterwise(xs, cs, ncomp = 1, lambda = 0, max_iter = 1)
z <- multivarious::transform(f$preproc, xs[[1]], colind = f$block_indices[[1]])
stopifnot(identical(unname(z), unname(xs[[1]])))
cat("C9: returned preprocessor leaves nonzero-mean training block unchanged\n")

# C10: variable local rank and padded spatial rank cannot be added.
expect_error_text(
  penalized_mfa_clusterwise(
    list(matrix(rnorm(20), 10, 2), matrix(rnorm(50), 10, 5)),
    list(cbind(1:2,0,0), cbind(1:5,1,0)), ncomp = 3, lambda = 1, max_iter = 1),
  "number of rows are not compatible for +", "C10")

cat("All 7 confirmed core-method repro checks passed.\n")
