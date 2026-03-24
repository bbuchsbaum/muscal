# Autoplot for Linked MFA

Creates a ggplot2 score plot for Linked MFA results.

## Usage

``` r
# S3 method for class 'linked_mfa'
autoplot(
  object,
  dims = c(1, 2),
  labels = FALSE,
  color = "coverage",
  alpha = 0.8,
  size = 2,
  ...
)
```

## Arguments

- object:

  An object of class \`linked_mfa\`.

- dims:

  Integer vector of length 2 specifying components to plot.

- labels:

  Logical; if \`TRUE\`, add text labels.

- color:

  Optional variable for point coloring. If \`"coverage"\`, colors by
  number of blocks contributing to each observation.

- alpha:

  Point transparency.

- size:

  Point size.

- ...:

  Additional arguments (unused).

## Value

A ggplot object.
