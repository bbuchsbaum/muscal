# Predict Clusterwise Penalized MFA Scores or Reconstructions

Predict Clusterwise Penalized MFA Scores or Reconstructions

## Usage

``` r
# S3 method for class 'penalized_mfa_clusterwise'
predict(
  object,
  new_data,
  block = NULL,
  type = c("scores", "reconstruction"),
  preprocess = TRUE,
  ...
)
```

## Arguments

- object:

  A current \`penalized_mfa_clusterwise\` fit.

- new_data:

  New observations; see \[project.penalized_mfa_clusterwise()\].

- block:

  Block name/index for a single input matrix.

- type:

  \`"scores"\` or \`"reconstruction"\`. Reconstruction is returned in
  raw training feature units via the fitted inverse transformation.

- preprocess:

  Whether new observations need the fitted forward transform.

- ...:

  Unused.

## Value

A matrix or complete named list of per-block matrices.
