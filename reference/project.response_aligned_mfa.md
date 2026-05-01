# Project New Rows into a Response-Aligned MFA Space

Project New Rows into a Response-Aligned MFA Space

## Usage

``` r
# S3 method for class 'response_aligned_mfa'
project(
  x,
  new_data,
  block,
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

- x:

  A fitted \`response_aligned_mfa\` object.

- new_data:

  Numeric matrix/data.frame of new rows for a known predictor block.

- block:

  Character or integer identifying the block used for projection.

- new_response:

  Optional observed multivariate response rows used to refine the
  latent-score solve when \`conditional = TRUE\`.

- response_weight:

  Optional scalar or vector of non-negative row weights paired with
  \`new_response\`.

- new_anchor_map:

  Optional anchor assignments for the new rows. May be an integer vector
  of anchor-state ids (with \`NA\` for unanchored rows) or a
  non-negative row-stochastic matrix.

- anchor_weight:

  Optional scalar or vector of non-negative row weights paired with
  \`new_anchor_map\`.

- conditional:

  Logical; if \`TRUE\`, allow conditional completion with
  \`new_response\`. The default \`FALSE\` keeps projection pure with
  respect to the response target.

- preprocess:

  Logical; if \`TRUE\` (default), applies the fitted predictor and
  response preprocessing pipelines before projection.

- ...:

  Unused.

## Value

A numeric matrix of latent scores.
