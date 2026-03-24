# Verify duality between ROI- and RV-space coordinates

This function checks whether the squared norms of the ROI-space
coordinates (F-scores) are equal to the squared norms of the RV-space
coordinates (T-scores) for each dimension, which is a fundamental
property of COVSTATIS.

## Usage

``` r
check_duality(x, tol = 1e-10)
```

## Arguments

- x:

  A fitted \`covstatis\` object.

- tol:

  The tolerance for the comparison.

## Value

\`TRUE\` if the duality property holds within the given tolerance.

## Examples

``` r
# Create a list of correlation matrices
Xs <- lapply(1:5, function(i) matrix(rnorm(10*10), 10, 10))
Xs <- lapply(Xs, cor)
# Apply COVSTATIS
res <- covstatis(Xs, ncomp=3)
# Check duality
check_duality(res)
#> [1] TRUE
```
