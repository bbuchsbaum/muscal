# Integrated Principal Components Analysis (generic)

A generic front-end for \*\*integrated principal components analysis\*\*
(iPCA) for multiple aligned data blocks. iPCA models each block with a
matrix-normal covariance structure with shared row covariance and
block-specific column covariance.

## Usage

``` r
ipca(
  data,
  preproc,
  ncomp = 2,
  lambda = 1,
  method = c("auto", "gram", "dense"),
  max_iter = 100,
  tol = 1e-06,
  normalize_trace = TRUE,
  use_future = FALSE,
  eig_solver = c("auto", "full", "truncated"),
  eig_rank = NULL,
  eig_trunc_min_n = 400,
  ...
)

# S3 method for class 'list'
ipca(
  data,
  preproc = multivarious::center(),
  ncomp = 2,
  lambda = 1,
  method = c("auto", "gram", "dense"),
  max_iter = 100,
  tol = 1e-06,
  normalize_trace = TRUE,
  use_future = FALSE,
  eig_solver = c("auto", "full", "truncated"),
  eig_rank = NULL,
  eig_trunc_min_n = 400,
  ...
)

# S3 method for class 'multiblock'
ipca(
  data,
  preproc = multivarious::center(),
  ncomp = 2,
  lambda = 1,
  method = c("auto", "gram", "dense"),
  max_iter = 100,
  tol = 1e-06,
  normalize_trace = TRUE,
  use_future = FALSE,
  eig_solver = c("auto", "full", "truncated"),
  eig_rank = NULL,
  eig_trunc_min_n = 400,
  .init_state = NULL,
  .return_state = FALSE,
  ...
)
```

## Arguments

- data:

  A data object for which an iPCA method is defined. Typically a
  \`list\` of matrices/data frames or a \`multiblock\` object.

- preproc:

  A preprocessing pipeline from the multivarious package. Each block is
  preprocessed independently.

- ncomp:

  Integer; number of joint components to return.

- lambda:

  Positive numeric scalar (recycled) or vector of length equal to the
  number of blocks. These are the multiplicative Frobenius penalties.

- method:

  One of \`"auto"\`, \`"gram"\`, or \`"dense"\`. - \`"auto"\` is
  recommended: it switches to sample-space updates when any block has
  \`p_k \> n\`. - \`"gram"\` is typically preferred for high-dimensional
  blocks. - \`"dense"\` can be faster when all blocks are moderate in
  size.

- max_iter:

  Integer; maximum Flip-Flop iterations.

- tol:

  Positive numeric convergence tolerance.

- normalize_trace:

  Logical; if \`TRUE\`, enforce \`mean(diag(Sigma)) = 1\` after each
  iteration for numerical stability.

- use_future:

  Logical; if \`TRUE\`, block-wise updates are parallelized via
  \`furrr::future_map()\` when available.

- eig_solver:

  Character; eigensolver policy for dense-mode Sigma updates. One of
  \`"auto"\`, \`"full"\`, or \`"truncated"\`.

- eig_rank:

  Optional integer; truncated eigendecomposition rank.

- eig_trunc_min_n:

  Integer; minimum \`n\` at which \`eig_solver = "auto"\` switches to
  truncated eigendecomposition.

- ...:

  Additional arguments passed to the underlying method.

- .init_state:

  Optional list of warm-start state from a previous fit.

- .return_state:

  Logical; if `TRUE`, include warm-start state in the result.

## Value

An object inheriting from class \`ipca\`.

## Details

In practice, iPCA is most useful when you want: - a shared sample-space
representation across all blocks (joint scores), - block-specific
feature loadings, and - stable integration in high-dimensional regimes
(\`p_k \>\> n\`).

Relative to methods such as MCCA/MFA, iPCA is often more robust when
block covariance/noise structures differ, but can be slower and is not
always best for purely correlation-maximizing objectives.

`ipca.list()` converts a list of blocks to a `multiblock` object and
then dispatches to `ipca.multiblock()`.

`ipca.multiblock()` implements multiplicative Frobenius iPCA from Tang
and Allen (2021) using a Flip-Flop algorithm. In `method = "gram"` mode,
fitting is performed in sample space (`n x n`) and avoids `p_k x p_k`
eigendecompositions. This is exact for the multiplicative Frobenius
updates.

## References

Tang, T. M. and Allen, G. I. (2021). Integrated Principal Components
Analysis. \*Journal of Machine Learning Research\*, 22(198), 1-81.

## Examples

``` r
# \donttest{
set.seed(1)
X <- list(
  X1 = matrix(rnorm(60 * 40), 60, 40),
  X2 = matrix(rnorm(60 * 80), 60, 80),
  X3 = matrix(rnorm(60 * 50), 60, 50)
)
fit <- ipca(X, ncomp = 3, lambda = 1)
#> Applying the same preprocessor definition independently to each block.
stopifnot(ncol(multivarious::scores(fit)) == 3)
# }
```
