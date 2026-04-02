# Predict from an Anchored MFA Fit

Predict from an Anchored MFA Fit

## Usage

``` r
# S3 method for class 'anchored_mfa'
predict(
  object,
  new_data,
  block,
  type = c("response", "scores", "reconstruction"),
  preprocess = TRUE,
  ...
)
```

## Arguments

- object:

  A fitted \`anchored_mfa\` object.

- new_data:

  Numeric matrix/data.frame of new rows from a known auxiliary block.

- block:

  Character or integer identifying the auxiliary block.

- type:

  One of \`"response"\`, \`"scores"\`, or \`"reconstruction"\`.

- preprocess:

  Logical; if \`TRUE\` (default), applies the fitted block preprocessing
  pipeline before projection.

- ...:

  Unused.

## Value

A numeric matrix of predicted anchor responses, latent scores, or
reconstructed block rows.
