# Plot COVSTATIS Results

Displays ROI scores or subject G-scores from a COVSTATIS analysis.

## Usage

``` r
# S3 method for class 'covstatis'
plot(
  x,
  dims = c(1, 2),
  type = c("roi", "subjects"),
  labels = FALSE,
  col = NULL,
  pch = 19,
  cex = 1,
  main = NULL,
  ...
)
```

## Arguments

- x:

  A fitted \`covstatis\` object.

- dims:

  Integer vector of length 2 specifying which components to plot.

- type:

  \`"roi"\` for ROI-level scores, \`"subjects"\` for subject G-scores.

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
