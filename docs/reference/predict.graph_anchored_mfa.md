# Predict from a Graph-Anchored MFA Fit

Predict from a Graph-Anchored MFA Fit

## Usage

``` r
# S3 method for class 'graph_anchored_mfa'
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

  A fitted \`graph_anchored_mfa\` object.

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

A numeric matrix of predicted \`Y\` rows when \`type = "response"\`,
reconstructed block rows when \`type = "reconstruction"\`, or latent
scores when \`type = "scores"\`.
