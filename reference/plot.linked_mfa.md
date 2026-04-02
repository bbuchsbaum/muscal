# Plot Linked MFA Scores

Displays observations in the Linked MFA shared score space, optionally
indicating block coverage.

## Usage

``` r
# S3 method for class 'linked_mfa'
plot(
  x,
  dims = c(1, 2),
  show_coverage = TRUE,
  labels = FALSE,
  col = NULL,
  pch = 19,
  cex = 1,
  main = "Linked MFA Score Plot",
  ...
)
```

## Arguments

- x:

  An object of class \`linked_mfa\`.

- dims:

  Integer vector of length 2 specifying which components to plot.

- show_coverage:

  Logical; if \`TRUE\`, color points by number of blocks contributing to
  each observation.

- labels:

  Logical; if \`TRUE\`, label points with row names or indices.

- col:

  Point colors. Overrides \`show_coverage\` coloring if provided.

- pch:

  Point character.

- cex:

  Point size.

- main:

  Plot title.

- ...:

  Additional arguments passed to \[plot()\].

## Value

Invisibly returns the score matrix (for the plotted dimensions).
