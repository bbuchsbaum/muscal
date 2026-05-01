# Project New Rows into an Aligned RRR Space

Project New Rows into an Aligned RRR Space

## Usage

``` r
# S3 method for class 'aligned_rrr'
project(x, new_data, block, preprocess = TRUE, ...)
```

## Arguments

- x:

  A fitted \`aligned_rrr\` object.

- new_data:

  Numeric matrix/data.frame of new rows for a known predictor block.

- block:

  Character or integer identifying the block used for projection.

- preprocess:

  Logical; if \`TRUE\` (default), applies the fitted predictor
  preprocessing pipeline before projection.

- ...:

  Unused.

## Value

A numeric matrix of latent scores.
