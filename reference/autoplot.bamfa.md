# Autoplot for BaMFA

Creates a ggplot2 visualization of BaMFA global scores.

## Usage

``` r
# S3 method for class 'bamfa'
autoplot(
  object,
  dims = c(1, 2),
  block = "mean",
  labels = FALSE,
  color = NULL,
  alpha = 0.8,
  size = 2,
  ...
)
```

## Arguments

- object:

  A fitted \`bamfa\` object.

- dims:

  Integer vector of length 2 specifying components to plot.

- block:

  Which block's scores to show, or \`"mean"\` for averaged scores.

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
