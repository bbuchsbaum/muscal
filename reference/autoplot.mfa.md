# Autoplot for MFA

Creates a ggplot2 score plot for MFA results.

## Usage

``` r
# S3 method for class 'mfa'
autoplot(
  object,
  dims = c(1, 2),
  labels = FALSE,
  color = NULL,
  alpha = 0.8,
  size = 2,
  ...
)
```

## Arguments

- object:

  An object of class \`mfa\`.

- dims:

  Integer vector of length 2 specifying components to plot.

- labels:

  Logical; if \`TRUE\`, add text labels.

- color:

  Optional variable for point coloring (vector of length n).

- alpha:

  Point transparency.

- size:

  Point size.

- ...:

  Additional arguments (unused).

## Value

A ggplot object.
