# Plot Variable Loadings

Displays variable contributions to components, optionally as a
correlation circle or bar chart.

## Usage

``` r
plot_loadings(
  x,
  dims = c(1, 2),
  type = c("circle", "bar"),
  component = 1,
  color_by = "block",
  top_n = 20,
  ...
)
```

## Arguments

- x:

  A fitted \`mfa\` or \`linked_mfa\` object.

- dims:

  Integer vector of length 2 for component axes.

- type:

  \`"circle"\` for correlation circle, \`"bar"\` for bar chart of
  loadings on a single component.

- component:

  For \`type = "bar"\`, which component to show.

- color_by:

  For MFA, color variables by \`"block"\` or \`NULL\`.

- top_n:

  For \`type = "bar"\`, show only top N variables by absolute loading.

- ...:

  Additional arguments.

## Value

A ggplot object or base R plot.
