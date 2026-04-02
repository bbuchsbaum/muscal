# Plot BaDA Results

Displays scores from a Barycentric Discriminant Analysis, colored by
class.

## Usage

``` r
# S3 method for class 'bada'
plot(
  x,
  dims = c(1, 2),
  labels = FALSE,
  show_barycenters = TRUE,
  col = NULL,
  pch = 19,
  cex = 1,
  main = "BaDA Score Plot",
  ...
)
```

## Arguments

- x:

  A fitted \`bada\` object.

- dims:

  Integer vector of length 2 specifying which components to plot.

- labels:

  Logical; if \`TRUE\`, label points.

- show_barycenters:

  Logical; if \`TRUE\`, show class centroids.

- col:

  Point colors. If \`NULL\`, colors by class label.

- pch:

  Point character.

- cex:

  Point size.

- main:

  Plot title.

- ...:

  Additional arguments passed to \[plot()\].

## Value

Invisibly returns the plotted scores.
