# Project New Rows into a Graph-Anchored MFA Space

Project New Rows into a Graph-Anchored MFA Space

## Usage

``` r
# S3 method for class 'graph_anchored_mfa'
project(x, new_data, block, preprocess = TRUE, ...)
```

## Arguments

- x:

  A fitted \`graph_anchored_mfa\` object.

- new_data:

  Numeric matrix/data.frame of new rows for a known auxiliary block.

- block:

  Character or integer identifying the auxiliary block whose loading
  matrix should be used.

- preprocess:

  Logical; if \`TRUE\` (default), applies the fitted block preprocessing
  pipeline before projection.

- ...:

  Unused.

## Value

A numeric matrix of latent scores.
