# Autoplot for Penalized MFA

Creates a ggplot2 visualization of Penalized MFA scores.

## Usage

``` r
# S3 method for class 'penalized_mfa'
autoplot(
  object,
  dims = c(1, 2),
  block = "consensus",
  labels = FALSE,
  color = NULL,
  alpha = 0.8,
  size = 2,
  ...
)
```

## Arguments

- object:

  A fitted \`penalized_mfa\` object.

- dims:

  Integer vector of length 2 specifying components to plot.

- block:

  Which block's scores to show, or \`"consensus"\` for averaged.

- labels:

  Logical; if \`TRUE\`, add text labels.

- color:

  Optional variable for point coloring.

- alpha:

  Point transparency.

- size:

  Point size.

- ...:

  Additional arguments (unused).

## Value

A ggplot object.
