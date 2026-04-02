# STATIS for Covariance Matrices (generic)

A generic function implementing the STATIS approach for a collection of
covariance matrices. Method implementations should return an object of
class \``covstatis`\`.

## Usage

``` r
covstatis(data, ncomp = 2, normalize = TRUE, dcenter = TRUE, ...)

# S3 method for class 'list'
covstatis(data, ncomp = 2, normalize = TRUE, dcenter = TRUE, ...)
```

## Arguments

- data:

  A data object for which a BaDA method is defined (see individual
  methods for details).

- ncomp:

  Integer; number of components to compute.

- normalize:

  Logical; if \`TRUE\` (default) each covariance matrix is scaled to
  have unit Frobenius norm before analysis.

- dcenter:

  Logical; if \`TRUE\` (default) each covariance matrix is double
  centred (Gower transformation) prior to analysis.

- ...:

  Additional arguments passed to the method.

## Value

An object inheriting from class \`covstatis\`.

## Details

The `covstatis.list` method implements STATIS analysis for a list of
covariance matrices.

The method operates as follows:

1.  Optionally double-centers each matrix (Gower transformation)

2.  Optionally normalizes each matrix to unit Frobenius norm

3.  Computes an RV matrix of inner products between matrices

4.  Determines optimal weights (`alpha`) via first eigenvector of the RV
    matrix. These weights are accessible in the output and can be useful
    for outlier diagnostics.

5.  Creates a compromise matrix as weighted sum of normalized matrices

6.  Performs eigendecomposition on the compromise matrix

The `partial_scores` in the returned object are a list of matrices (one
for each subject), with each matrix containing the ROI-level factor
scores (ROI x D).

Additional arguments can be passed via `...`:

- labels:

  Optional character vector of labels for the rows/columns of the
  covariance matrices. If NULL, tries to use row names from the first
  matrix, or generates sequential labels.

- norm_method:

  The normalization method to apply to each matrix. One of `"frobenius"`
  (default when `normalize=TRUE`), `"mfa"`, or `"none"`. Frobenius norm
  scales each matrix to have a total sum-of-squares of 1. MFA (multiple
  factor analysis) scales each matrix so its first eigenvalue is 1.
  `"none"` applies no normalization.

## See also

\[covstatis.list\] for the reference implementation.

## Examples

``` r
# Create a list of correlation matrices
Xs <- lapply(1:5, function(i) matrix(rnorm(10*10), 10, 10))
Xs <- lapply(Xs, cor)

# Apply COVSTATIS
res <- covstatis(Xs, ncomp=3)

# Project a new correlation matrix
new_mat <- cor(matrix(rnorm(10*10), 10, 10))
proj <- project_cov(res, new_mat)
```
