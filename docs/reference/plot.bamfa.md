# Plot BaMFA Results

Displays global scores from a BaMFA analysis.

## Usage

``` r
# S3 method for class 'bamfa'
plot(
  x,
  dims = c(1, 2),
  block = "mean",
  labels = FALSE,
  col = NULL,
  pch = 19,
  cex = 1,
  main = "BaMFA Global Scores",
  ...
)
```

## Arguments

- x:

  A fitted \`bamfa\` object.

- dims:

  Integer vector of length 2 specifying which components to plot.

- block:

  Which block's scores to show, or \`"mean"\` for averaged scores.

- labels:

  Logical; if \`TRUE\`, label points.

- col:

  Point colors.

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
