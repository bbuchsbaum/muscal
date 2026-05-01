# Predict from a Response-Aligned MFA Fit

Predict from a Response-Aligned MFA Fit

## Usage

``` r
# S3 method for class 'response_aligned_mfa'
predict(
  object,
  new_data,
  block,
  type = c("response", "scores", "reconstruction"),
  new_response = NULL,
  response_weight = NULL,
  new_anchor_map = NULL,
  anchor_weight = NULL,
  conditional = FALSE,
  preprocess = TRUE,
  ...
)
```

## Arguments

- object:

  A fitted \`response_aligned_mfa\` object.

- new_data:

  Numeric matrix/data.frame of new rows from a known predictor block.

- block:

  Character or integer identifying the block.

- type:

  One of \`"response"\`, \`"scores"\`, or \`"reconstruction"\`.

- new_response:

  Optional observed multivariate response rows used to refine the
  latent-score solve when \`conditional = TRUE\`.

- response_weight:

  Optional scalar or vector of non-negative row weights paired with
  \`new_response\`.

- new_anchor_map:

  Optional anchor assignments for the new rows used to refine the
  latent-score solve.

- anchor_weight:

  Optional scalar or vector of non-negative row weights paired with
  \`new_anchor_map\`.

- conditional:

  Logical; if \`TRUE\`, allow conditional completion with
  \`new_response\`. The default \`FALSE\` keeps prediction target-pure.

- preprocess:

  Logical; if \`TRUE\` (default), applies the fitted predictor and
  response preprocessing pipelines before projection.

- ...:

  Unused.

## Value

A numeric matrix of predicted responses, latent scores, or reconstructed
block rows.
