# Multiblock Canonical Correlation Analysis (generic)

A generic front-end for \*\*multiblock canonical correlation
analysis\*\* (MCCA), also known as generalized CCA (GCCA). Concrete
methods should estimate a compromise score space shared across multiple
data blocks that have the same observations (rows) but potentially
different variables (columns).

## Usage

``` r
mcca(
  data,
  preproc,
  ncomp = 2,
  ridge = 1e-06,
  block_weights = NULL,
  use_future = FALSE,
  ...
)

# S3 method for class 'list'
mcca(
  data,
  preproc = multivarious::center(),
  ncomp = 2,
  ridge = 1e-06,
  block_weights = NULL,
  use_future = FALSE,
  ...
)

# S3 method for class 'multiblock'
mcca(
  data,
  preproc = multivarious::center(),
  ncomp = 2,
  ridge = 1e-06,
  block_weights = NULL,
  use_future = FALSE,
  ...
)
```

## Arguments

- data:

  A data object for which an MCCA method is defined. Typically a
  \`list\` of matrices/data frames or a \`multiblock\` object.

- preproc:

  A preprocessing pipeline from the multivarious package. Each block is
  preprocessed independently.

- ncomp:

  Integer; the number of canonical components to compute.

- ridge:

  Non-negative numeric scalar (or vector of length equal to the number
  of blocks) controlling ridge stabilization. The effective ridge used
  per block is scaled by the block's overall energy so the default works
  across variable scalings.

- block_weights:

  Optional numeric vector of non-negative weights (length = number of
  blocks) controlling each block's influence on the compromise.

- use_future:

  Logical; if \`TRUE\`, block-wise computations are performed via
  \`furrr::future_map()\` when available.

- ...:

  Additional arguments passed to the underlying method.

## Value

An object inheriting from class \`mcca\`.

## Details

The reference implementation in this package uses a ridge-stabilized
MAXVAR GCCA formulation that remains well-posed when blocks have more
variables than rows (\\p \gg n\\) or are rank-deficient after
preprocessing (e.g., centering).

The `mcca.list` method applies MCCA to a list of data matrices or data
frames. This method first converts the list to a `multiblock` object and
then calls `mcca.multiblock()`.

The `mcca.multiblock` method implements a ridge-stabilized MAXVAR
generalized CCA. Internally it computes a compromise score space by
eigendecomposition of a weighted sum of block projection operators:
\$\$H = \sum_s w_s \\ X_s (X_s^T X_s + \kappa_s I)^{-1} X_s^T.\$\$

The ridge parameter is applied per block as \\\kappa_s = \text{ridge}\_s
\cdot \bar{d}\_s\\, where \\\bar{d}\_s\\ is the mean diagonal of the
block Gram matrix \\X_s X_s^T\\. This scaling makes the default ridge
behave reasonably across different variable scalings.

## References

Kettenring, J. R. (1971). Canonical analysis of several sets of
variables. \*Biometrika\*, 58(3), 433–451.

## Examples

``` r
# \donttest{
set.seed(1)
X <- replicate(3, matrix(rnorm(50 * 20), 50, 20), simplify = FALSE)
fit <- mcca(X, ncomp = 3)
#> Applying the same preprocessor definition independently to each block.
# }
# \donttest{
set.seed(1)
blocks <- replicate(3, matrix(rnorm(20 * 50), 20, 50), simplify = FALSE)
fit <- mcca(blocks, ncomp = 2)
#> Applying the same preprocessor definition independently to each block.
stopifnot(ncol(multivarious::scores(fit)) == 2)
# }
```
