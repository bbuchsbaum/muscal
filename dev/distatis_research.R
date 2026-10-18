# Research reproducer, not a public DISTATIS implementation.
# Run from the repository root: Rscript --vanilla dev/distatis_research.R
# Table 1 of Abdi, Valentin, Chollet & Chrea (2007), transcribed from
# Analyzing_assessors_and_products_in_sort.pdf, PDF page 8.
sorting <- matrix(c(
  1,4,3,4,1,1,2,2,1,3,
  4,5,2,5,2,3,1,1,4,3,
  3,1,2,3,2,4,3,1,1,2,
  4,2,3,3,1,1,1,2,1,4,
  1,5,3,5,2,1,1,2,1,3,
  2,3,1,1,3,5,4,4,3,1,
  1,4,3,4,1,1,2,2,2,4,
  5,2,4,2,4,2,5,3,4,5
), nrow = 8L, byrow = TRUE,
dimnames = list(c("Affligen", "Budweiser", "Buckler Blonde", "Killian",
                   "St. Landelin", "Buckler Highland", "Fruit Defendu", "EKU28"),
                 paste0("Assessor_", 1:10)))

# Deliberately use explicit centering matrices, independently of double_center().
H <- diag(8) - matrix(1/8, 8, 8)
distances <- lapply(seq_len(ncol(sorting)), function(k) {
  1 * outer(sorting[, k], sorting[, k], `!=`)
})
grams <- lapply(distances, function(D) -0.5 * H %*% (D^2) %*% H)
tables <- lapply(grams, function(B) {
  B / eigen(B, symmetric = TRUE, only.values = TRUE)$values[1L]
})
M <- vapply(tables, as.vector, numeric(64))
raw_C <- crossprod(M)
RV <- raw_C / outer(sqrt(diag(raw_C)), sqrt(diag(raw_C)))
rv_eig <- eigen(RV, symmetric = TRUE)
w <- rv_eig$vectors[, 1L]
w <- w * sign(sum(w))
alpha <- w / sum(w)
compromise <- Reduce(`+`, Map(`*`, tables, alpha))
eig <- eigen(compromise, symmetric = TRUE)
keep <- which(eig$values > 1e-10)
F <- sweep(eig$vectors[, keep, drop = FALSE], 2, sqrt(eig$values[keep]), `*`)
projection <- sweep(eig$vectors[, keep, drop = FALSE], 2,
                    sqrt(eig$values[keep]), `/`)
partial <- lapply(tables, function(B) B %*% projection)

# Independent, rounded published targets: eq. 21, eq. 22, and Table 2.
published_RV <- matrix(c(
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
published_alpha <- c(.100,.101,.109,.101,.099,.116,.101,.109,.074,.090)
published_roots <- c(.66,.49,.40,.34,.23,.13,.05)
published_F <- t(matrix(c(
  -.39,.23,.28,-.16,-.15,.30,-.40,.29,
  .07,-.37,-.22,-.03,-.24,.38,.15,.26,
  -.03,-.02,-.07,.09,-.04,-.40,-.01,.48,
  .03,-.25,.45,.08,-.22,-.10,.07,-.06,
  .16,.06,.07,-.42,-.01,-.05,.11,.06,
  .02,-.20,.05,-.05,.26,.01,-.13,.04,
  .16,.02,-.01,.03,-.06,.00,-.14,-.00
), 7, 8, byrow = TRUE))
F_aligned <- sweep(F, 2, sign(colSums(F * published_F)), `*`)
stopifnot(max(abs(RV - published_RV)) < .0051,
          max(abs(alpha - published_alpha)) < .00051,
          max(abs(eig$values[keep] - published_roots)) < .0051,
          max(abs(F_aligned - published_F)) < .0051,
          max(abs(Reduce(`+`, Map(`*`, partial, alpha)) - F)) < 1e-12,
          max(abs(rowSums(compromise))) < 1e-12)

# Euclidean coordinates with non-binary distances reveal accidental re-squaring.
X <- cbind(c(0, 1, 3, 5), c(0, 2, 1, 4))
D <- as.matrix(stats::dist(X))
H4 <- diag(4) - matrix(1/4, 4, 4)
B <- -0.5 * H4 %*% D^2 %*% H4
stopifnot(max(abs(B - tcrossprod(scale(X, scale = FALSE)))) < 1e-12)

# Compare the current checkout, rather than a potentially stale installation.
env <- new.env(parent = globalenv())
env$`%||%` <- function(x, y) if (is.null(x)) y else x
sys.source("R/covstatis.R", envir = env)
current <- env$covstatis.list(grams, ncomp = 7, norm_method = "mfa", dcenter = FALSE)
stopifnot(max(abs(current$C - raw_C)) < 1e-12)

# Current helpers use compromise roots to scale the assessor eigenvectors.
current_G <- env$rv_subject_scores(current)
current_G_norms <- colSums(current_G^2)
standard_G_norms <- current$eigC$values[seq_len(ncol(current_G))]

# Rank-one inputs are valid sorting data. Capture the current scalar-diag error
# without normalization (unit-root normalization can hide the scalar-diag bug).
H3 <- diag(3) - matrix(1/3, 3, 3)
D3 <- 1 * outer(c(1, 1, 2), c(1, 1, 2), `!=`)
B3 <- -0.5 * H3 %*% D3 %*% H3
rank_one <- tryCatch({
  env$covstatis.list(list(B3, B3), ncomp = 1, norm_method = "none", dcenter = FALSE)
  "NO ERROR"
}, error = function(e) conditionMessage(e))
one_component <- tryCatch({
  env$covstatis.list(grams, ncomp = 1, norm_method = "mfa", dcenter = FALSE)
  "NO ERROR"
}, error = function(e) conditionMessage(e))

cat("DISTATIS research reproducer\n")
cat("R version:", as.character(getRversion()), "\n")
cat("All published-target and geometry assertions passed.\n\n")
cat("Assessor weights (paper method vs current covstatis MFA):\n")
print(data.frame(assessor = 1:10, paper_method = alpha,
                 covstatis_mfa = current$alpha), row.names = FALSE, digits = 7)
cat("\nMaximum weight difference:", max(abs(alpha - current$alpha)), "\n")
cat("RV diagonal:", diag(RV), "\n")
cat("Current MFA interstructure diagonal:", diag(current$C), "\n")
cat("RV leading roots:", head(rv_eig$values, 3), "\n")
cat("Compromise positive roots:", eig$values[keep], "\n")
cat("First four compromise axes, percent inertia:",
    100 * sum(eig$values[1:4]) / sum(eig$values[keep]), "\n")
cat("Compromise quality (leading RV root / trace):",
    rv_eig$values[1] / sum(rv_eig$values), "\n")
cat("Maximum published RV rounding error:", max(abs(RV - published_RV)), "\n")
cat("Maximum published score rounding error:", max(abs(F_aligned - published_F)), "\n")
cat("Partial-score barycenter error:",
    max(abs(Reduce(`+`, Map(`*`, partial, alpha)) - F)), "\n")
cat("\nCurrent assessor coordinate squared norms:", current_G_norms, "\n")
cat("Standard interstructure coordinate squared norms:", standard_G_norms, "\n")
cat("Current rank-one fit outcome:", rank_one, "\n")
cat("Current one-component MFA beer fit outcome:", one_component, "\n")
cat("\nDistatisR installed:", requireNamespace("DistatisR", quietly = TRUE), "\n")
cat("No DistatisR numerical parity claim is made by this script.\n")
