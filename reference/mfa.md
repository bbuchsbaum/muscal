# Multiple Factor Analysis (generic)

A generic front–end for \*\*Multiple Factor Analysis\*\* (MFA). Concrete
methods should perform the estimation and return an object inheriting
from class \``mfa`\`.

## Usage

``` r
mfa(data, preproc, ncomp = 2, normalization = "MFA", A = NULL, M = NULL, ...)

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

- missing:

  Missing-data handling mode. `"error"` preserves the historical
  complete-data contract. `"regularized"` and `"em"` run iterative MFA
  imputation before the final complete-data MFA fit. `"nipals"` is
  reserved for a future direct available-data backend.

- ncp_impute:

  Integer; number of MFA components used inside the iterative imputation
  loop.

- missing_tol:

  Numeric convergence tolerance for iterative imputation.

- missing_maxiter:

  Integer maximum number of imputation iterations.

- return_imputed:

  Logical; if `TRUE`, attach completed blocks as `imputed_data`.

## Value

An object of class \`mfa\`.

## Details

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

## References

Abdi, H., Williams, L. J., & Valentin, D. (2013). \*Multiple factor
analysis: principal component analysis for multi‑table and multi‑block
data sets\*. \*\*Wiley Interdisciplinary Reviews: Computational
Statistics, 5\*\*(2), 149–179.

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
