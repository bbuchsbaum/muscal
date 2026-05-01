# Orient Component Signs for Fitted muscal Models

Post-process a fitted model by applying a shared per-component sign
choice across scores, loadings, and cached derived quantities. This does
not refit the model. Reconstructions remain unchanged under the sign
flip; projected scores follow the new orientation.

## Usage

``` r
orient_components(
  x,
  basis = NULL,
  reference = NULL,
  signs = NULL,
  rule = c("max_abs", "sum"),
  tol = 1e-08,
  ...
)

# S3 method for class 'mfa'
orient_components(
  x,
  basis = NULL,
  reference = NULL,
  signs = NULL,
  rule = c("max_abs", "sum"),
  tol = 1e-08,
  ...
)

# S3 method for class 'aligned_mfa'
orient_components(
  x,
  basis = NULL,
  reference = NULL,
  signs = NULL,
  rule = c("max_abs", "sum"),
  tol = 1e-08,
  ...
)

# S3 method for class 'linked_mfa'
orient_components(
  x,
  basis = NULL,
  reference = NULL,
  signs = NULL,
  rule = c("max_abs", "sum"),
  tol = 1e-08,
  ...
)

# S3 method for class 'aligned_interbattery'
orient_components(
  x,
  basis = NULL,
  reference = NULL,
  signs = NULL,
  rule = c("max_abs", "sum"),
  tol = 1e-08,
  ...
)

# S3 method for class 'mcca'
orient_components(
  x,
  basis = NULL,
  reference = NULL,
  signs = NULL,
  rule = c("max_abs", "sum"),
  tol = 1e-08,
  ...
)

# S3 method for class 'aligned_mcca'
orient_components(
  x,
  basis = NULL,
  reference = NULL,
  signs = NULL,
  rule = c("max_abs", "sum"),
  tol = 1e-08,
  ...
)

# S3 method for class 'response_aligned_mfa'
orient_components(
  x,
  basis = NULL,
  reference = NULL,
  signs = NULL,
  rule = c("max_abs", "sum"),
  tol = 1e-08,
  ...
)

# S3 method for class 'aligned_rrr'
orient_components(
  x,
  basis = NULL,
  reference = NULL,
  signs = NULL,
  rule = c("max_abs", "sum"),
  tol = 1e-08,
  ...
)

# S3 method for class 'ipca'
orient_components(
  x,
  basis = NULL,
  reference = NULL,
  signs = NULL,
  rule = c("max_abs", "sum"),
  tol = 1e-08,
  ...
)

# S3 method for class 'penalized_mfa'
orient_components(
  x,
  basis = NULL,
  reference = NULL,
  signs = NULL,
  rule = c("max_abs", "sum"),
  tol = 1e-08,
  ...
)

# Default S3 method
orient_components(
  x,
  basis = NULL,
  reference = NULL,
  signs = NULL,
  rule = c("max_abs", "sum"),
  tol = 1e-08,
  ...
)
```

## Arguments

- x:

  A fitted model object.

- basis:

  Basis used to choose the component signs. May be \`NULL\` (use a
  class-specific default), a character scalar naming a stored field in
  \`x\`, a numeric matrix/data frame, or a list of matrices that will be
  row-bound before computing the sign rule.

- reference:

  Optional reference used to align signs. May be another fitted model, a
  numeric matrix/data frame, or a list of matrices. When supplied, signs
  are chosen to maximize per-component agreement between the resolved
  basis for \`x\` and the resolved basis for \`reference\`.

- signs:

  Optional numeric vector of length equal to the number of components.
  Non-zero values are coerced to \`-1\` or \`+1\`. When supplied,
  \`basis\` and \`reference\` are ignored.

- rule:

  Rule used when \`reference = NULL\` and \`signs = NULL\`.
  \`"max_abs"\` makes the largest-magnitude entry in each basis column
  positive. \`"sum"\` makes the column sum positive.

- tol:

  Non-negative tolerance used when a sign decision is numerically
  ambiguous.

- ...:

  Unused; reserved for future extensions.

## Value

An oriented copy of \`x\`.

## Details

Supported first-wave classes are \`mfa\`, \`aligned_mfa\`,
\`anchored_mfa\` / \`linked_mfa\` (including graph-anchored variants),
\`mcca\`, \`aligned_mcca\` / \`anchored_mcca\`,
\`response_aligned_mfa\`, \`aligned_rrr\`, \`aligned_interbattery\`,
\`ipca\`, and \`penalized_mfa\`.

The sign convention can be supplied explicitly via \`signs\`, inferred
from a basis matrix or list of matrices via \`basis\`, or aligned to a
reference fit or matrix via \`reference\`.
