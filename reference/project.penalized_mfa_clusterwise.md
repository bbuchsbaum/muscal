# Project Clusterwise Penalized MFA Observations

Project Clusterwise Penalized MFA Observations

## Usage

``` r
# S3 method for class 'penalized_mfa_clusterwise'
project(x, new_data, block = NULL, preprocess = TRUE, ...)
```

## Arguments

- x:

  A current \`penalized_mfa_clusterwise\` fit with fitted preprocessing.

- new_data:

  One raw matrix/data.frame with \`block\`, or a complete block list.
  Named lists are matched by name; unnamed lists use training order.

- block:

  Block name or index for a single matrix.

- preprocess:

  Apply the fitted training transformation exactly once (default TRUE);
  FALSE declares that the supplied rows are already transformed.

- ...:

  Unused.

## Value

A score matrix for a single block, or a named list of score matrices. No
shared-row consensus is inferred across independent blocks.
