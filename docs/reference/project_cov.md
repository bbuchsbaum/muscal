# Project a Covariance Matrix (generic)

Generic for projecting a new covariance matrix onto a previously fitted
\`covstatis\` model.

## Usage

``` r
project_cov(x, new_data, ...)

# S3 method for class 'covstatis'
project_cov(x, new_data, ...)
```

## Arguments

- x:

  A model object for which a \`project_cov\` method exists (e.g. an
  object of class \`covstatis\`).

- new_data:

  A symmetric numeric matrix to be projected.

- ...:

  Additional arguments passed to methods.

## Value

A numeric matrix of projected scores.

A numeric matrix of projected scores with dimensions `nrow(new_data)` ×
`ncomp`.

## Details

For `covstatis` objects, a new covariance/correlation matrix is
transformed using the same preprocessing steps (centering and
normalization) as were applied during model fitting, then projected onto
the model space. This returns the "partial factor scores" or ROI-level
coordinates for the new subject.
