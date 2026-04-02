# Autoplot for COVSTATIS

Creates a ggplot2 visualization of COVSTATIS results.

## Usage

``` r
# S3 method for class 'covstatis'
autoplot(
  object,
  dims = c(1, 2),
  type = c("roi", "subjects"),
  labels = FALSE,
  color = NULL,
  alpha = 0.8,
  size = 2,
  ...
)
```

## Arguments

- object:

  A fitted \`covstatis\` object.

- dims:

  Integer vector of length 2 specifying components to plot.

- type:

  \`"roi"\` for ROI-level scores, \`"subjects"\` for subject G-scores.

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
