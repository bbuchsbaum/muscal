# Plot Loading Alignment for Penalized MFA

Visualizes how well block-specific loadings align with each other,
showing the effect of the penalty parameter.

## Usage

``` r
plot_loading_alignment(x, dims = c(1, 2), top_n = 10, ...)
```

## Arguments

- x:

  A fitted \`penalized_mfa\` object.

- dims:

  Integer vector of length 2 for component axes.

- top_n:

  Number of top variables to show per block.

- ...:

  Additional arguments.

## Value

A ggplot object or base R plot.
