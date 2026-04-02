# Project New Rows into an Aligned MFA Space

Project New Rows into an Aligned MFA Space

## Usage

``` r
# S3 method for class 'aligned_mfa'
project(x, new_data, block, preprocess = TRUE, ...)
```

## Arguments

- x:

  A fitted \`aligned_mfa\` object.

- new_data:

  Numeric matrix/data.frame of new rows for a known block.

- block:

  Character or integer identifying the block used for projection.

- preprocess:

  Logical; if \`TRUE\` (default), applies the fitted block preprocessing
  pipeline before projection.

- ...:

  Unused.

## Value

A numeric matrix of latent scores.
