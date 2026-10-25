# Multiple Factor Analysis (generic)

A generic front–end for \*\*Multiple Factor Analysis\*\* (MFA). Concrete
methods should perform the estimation and return an object inheriting
from class \``mfa`\`.

## Usage

``` r
mfa(data, preproc, ncomp = 2, normalization = "MFA", A = NULL, M = NULL, ...)

# S3 method for class 'muscal_block_source'
mfa(
  data,
  preproc = center(),
  ncomp = 2,
  normalization = c("MFA", "RV", "None", "Frob", "custom"),
  A = NULL,
  M = NULL,
  gram_rank = NULL,
  gram_tol = 1e-10,
  ...
)

# S3 method for class 'list'
mfa(
  data,
  preproc = center(),
  ncomp = 2,
  normalization = c("MFA", "RV", "None", "Frob", "custom"),
  A = NULL,
  M = NULL,
  ...
)

# S3 method for class 'multiblock'
mfa(
  data,
  preproc = center(),
  ncomp = 2,
  normalization = c("MFA", "RV", "None", "Frob", "custom"),
  A = NULL,
  M = NULL,
  missing = c("error", "regularized", "em", "nipals"),
  ncp_impute = ncomp,
  missing_tol = 1e-06,
  missing_maxiter = 100,
  return_imputed = FALSE,
  ...
)
```

## Arguments

- data:

  A data object for which a BaDA method is defined (see individual
  methods for details).

- preproc:

  A preprocessing pipeline from the multivarious package. Defaults to
  \`multivarious::center()\` in the concrete methods.

- ncomp:

  Integer; the number of discriminant components to compute.

- normalization:

  Character string specifying the block‑weighting scheme (see
  \[mfa.multiblock\]).

- A, M:

  Optional user‑supplied column/row weight matrices used when
  \`normalization = "custom"\`.

- ...:

  Additional arguments passed to the method.

- gram_rank:

  Block sources only: optional maximum number of Gram eigenpairs kept
  per block (one value or one per block); \`NULL\` keeps all, which is
  exact. The fraction of each block's variance left out is returned as
  \`gram_discarded\`.

- gram_tol:

  Block sources only: relative eigenvalue tolerance for dropping
  numerically null directions of each Gram matrix.

- missing:

  Missing-data handling mode. `"error"` (default) rejects `NA` values.
  `"regularized"` and `"em"` run iterative MFA imputation before the
  final complete-data MFA fit (see the *Missing data* section).
  `"nipals"` is reserved for a future direct available-data backend.

- ncp_impute:

  Integer; number of MFA components used inside the iterative imputation
  loop. Values above `min(n - 2, p - 1)` are capped.

- missing_tol:

  Numeric convergence tolerance for iterative imputation.

- missing_maxiter:

  Integer maximum number of imputation iterations.

- return_imputed:

  Logical; if `TRUE`, attach completed blocks as `imputed_data`.

## Value

An object of class \`mfa\`.

## Details

The \`mfa.muscal_block_source\` method fits MFA out of core from a
\[block_source()\]: each block is read once, chunk by chunk, compressed
to its row Gram matrix, fitted in memory, and read once more to map
loadings and correlation loadings back to the original features. Results
match the in-memory fit (up to the sign of each component). Only
column-local preprocessing is supported, missing values are not, and a
custom per-feature column metric \`A\` is rejected because it acts on
individual features.

The `mfa.list` method applies the MFA algorithm to a list of data
matrices or data frames. This method first converts the list to a
multiblock object and then calls `mfa.multiblock`.

The `mfa.multiblock` method implements Multiple Factor Analysis for a
collection of data blocks. This method handles data preprocessing, block
normalization, and integration of multiple data tables that share the
same observations.

Normalization options include:

- `MFA`: Scales each block by its first singular value (default)

- `RV`: Normalizes blocks based on RV matrix correlation

- `None`: No scaling applied

- `Frob`: Uses Frobenius norm for scaling

- `custom`: Uses custom weight matrices provided via A and M parameters

## Missing data

With `missing = "regularized"` or `"em"`, `NA`/`NaN` cells are imputed
by iterative MFA (Josse & Husson, 2012; the algorithm of
`missMDA::imputeMFA()`), and the returned model is an ordinary MFA fit
to the completed blocks. Missing cells start at their column means; each
iteration then refits the preprocessing, block weights, and an
`ncp_impute`-component MFA to the current completed data and replaces
only the missing cells with the low-rank reconstruction. `"em"` uses the
plain rank-`ncp_impute` reconstruction; `"regularized"` shrinks each
retained singular value `d_k` to `(d_k^2 - sigma2) / d_k`, where
`sigma2` is the residual inertia per residual degree of freedom, which
guards against overfitting when data are noisy or many cells are
missing.

Rows may be missing entire blocks, but every column needs at least one
observed value, blocks must share the same rows (see
[`linked_mfa()`](https://bbuchsbaum.github.io/muscal/reference/linked_mfa.md)
for blocks observed on different row sets), and infinite values are
rejected. `ncp_impute` is capped at `min(n - 2, p - 1)` with a warning,
and a warning is issued if the imputation does not converge within
`missing_maxiter` iterations. Convergence is the root-mean-square change
of the imputed cells, in units of each variable's observed standard
deviation, falling below `missing_tol`. Imputed cells influence the
final scores and loadings, so choose `ncp_impute` with care (it defaults
to `ncomp`). The result records
`missing = list(method, ncp_impute, iterations, converged, delta, mask)`.

## References

Abdi, H., Williams, L. J., & Valentin, D. (2013). \*Multiple factor
analysis: principal component analysis for multi‑table and multi‑block
data sets\*. \*\*Wiley Interdisciplinary Reviews: Computational
Statistics, 5\*\*(2), 149–179.

Josse, J., & Husson, F. (2012). Handling missing values in exploratory
multivariate data analysis methods. \*Journal de la Société Française de
Statistique, 153\*(2), 79–99.

## Examples

``` r
# \donttest{
# Apply MFA to a list of matrices
X <- replicate(5, { matrix(rnorm(10*10), 10, 10) }, simplify=FALSE)
res <- mfa(X, ncomp=3, normalization="MFA")
#> Applying the same preprocessor definition independently to each block.
#> normalization type:MFA
#> Registered S3 method overwritten by 'genpca':
#>   method                   from        
#>   transfer.cross_projector multivarious
# }
# \donttest{
# Create 5 random matrices of the same size
X <- replicate(5, { matrix(rnorm(10*10), 10, 10) }, simplify=FALSE)

# Apply MFA with MFA normalization
res <- mfa(X, ncomp=3, normalization="MFA")
#> Applying the same preprocessor definition independently to each block.
#> normalization type:MFA

# Project a block onto the model
p <- multivarious::project_block(res, X[[1]], 1)

# Verify number of components
stopifnot(ncol(multivarious::scores(res)) == 3)

# Create a classifier
labs <- letters[1:10]
cfier <- multivarious::classifier(res, new_data=do.call(cbind, X), labels=labs)
pred <- predict(cfier, do.call(cbind, X)[1:2,])

# Create a classifier using a specific block
cfier2 <- multivarious::classifier(res, new_data=X[[2]], labels=labs,
                                  colind=res$block_indices[[2]])
pred2 <- predict(cfier2, X[[2]][1:2,])
# }
```
