# Plot Class Barycenters for BaDA

Visualizes class centroids with optional confidence ellipses.

## Usage

``` r
plot_barycenters(x, dims = c(1, 2), ellipses = TRUE, level = 0.95, ...)
```

## Arguments

- x:

  A fitted \`bada\` object.

- dims:

  Integer vector of length 2 for component axes.

- ellipses:

  Logical; if \`TRUE\`, draw confidence ellipses around classes.

- level:

  Confidence level for ellipses (default 0.95).

- ...:

  Additional arguments.

## Value

A ggplot object or base R plot.
