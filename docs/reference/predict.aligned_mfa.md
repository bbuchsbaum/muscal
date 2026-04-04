# Predict from an Aligned MFA Fit

Predict from an Aligned MFA Fit

## Usage

``` r
# S3 method for class 'aligned_mfa'
predict(
  object,
  new_data,
  block,
  type = c("scores", "reconstruction"),
  preprocess = TRUE,
  ...
)
```

## Arguments

- object:

  A fitted \`aligned_mfa\` object.

- new_data:

  Numeric matrix/data.frame of new rows from a known block.

- block:

  Character or integer identifying the block.

- type:

  One of \`"scores"\` or \`"reconstruction"\`.

- preprocess:

  Logical; if \`TRUE\` (default), applies the fitted block preprocessing
  pipeline before projection.

- ...:

  Unused.

## Value

A numeric matrix of latent scores or reconstructed block rows.
