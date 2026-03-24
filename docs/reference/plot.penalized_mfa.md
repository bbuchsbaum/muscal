# Plot Penalized MFA Results

Displays consensus scores from a Penalized MFA analysis.

## Usage

``` r
# S3 method for class 'penalized_mfa'
plot(
  x,
  dims = c(1, 2),
  block = "consensus",
  labels = FALSE,
  col = NULL,
  pch = 19,
  cex = 1,
  main = "Penalized MFA Scores",
  ...
)
```

## Arguments

- x:

  A fitted \`penalized_mfa\` object.

- dims:

  Integer vector of length 2 specifying which components to plot.

- block:

  Which block's scores to show, or \`"consensus"\` for averaged.

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
