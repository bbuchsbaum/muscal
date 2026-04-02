# Autoplot for BaDA

Creates a ggplot2 visualization of BaDA scores, colored by class.

## Usage

``` r
# S3 method for class 'bada'
autoplot(
  object,
  dims = c(1, 2),
  labels = FALSE,
  show_barycenters = TRUE,
  alpha = 0.7,
  size = 2,
  ...
)
```

## Arguments

- object:

  A fitted \`bada\` object.

- dims:

  Integer vector of length 2 specifying components to plot.

- labels:

  Logical; if \`TRUE\`, add text labels.

- show_barycenters:

  Logical; if \`TRUE\`, show class centroids.

- alpha:

  Point transparency.

- size:

  Point size.

- ...:

  Additional arguments (unused).

## Value

A ggplot object.
