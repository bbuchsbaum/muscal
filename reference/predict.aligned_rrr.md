# Predict from an Aligned RRR Fit

Predict from an Aligned RRR Fit

## Usage

``` r
# S3 method for class 'aligned_rrr'
predict(
  object,
  new_data,
  block,
  type = c("response", "scores"),
  preprocess = TRUE,
  ...
)
```

## Arguments

- object:

  A fitted \`aligned_rrr\` object.

- new_data:

  Numeric matrix/data.frame of new rows from a known predictor block.

- block:

  Character or integer identifying the block.

- type:

  One of \`"response"\` or \`"scores"\`.

- preprocess:

  Logical; if \`TRUE\` (default), applies the fitted predictor
  preprocessing pipeline before projection.

- ...:

  Unused.

## Value

A numeric matrix of predicted responses or latent scores.
