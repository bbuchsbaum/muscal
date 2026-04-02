# Plot Variance Explained (Scree Plot)

Shows the proportion of variance captured by each component.

## Usage

``` r
plot_variance(x, type = c("bar", "line"), ...)

plot_variance.covstatis(x, type = c("bar", "line"), ...)

plot_variance.bada(x, type = c("bar", "line"), ...)

plot_variance.penalized_mfa(x, type = c("bar", "line"), ...)
```

## Arguments

- x:

  A fitted \`mfa\` or \`linked_mfa\` object.

- type:

  \`"bar"\` for bar chart, \`"line"\` for scree plot with cumulative.

- ...:

  Additional arguments.

## Value

A ggplot object or base R plot.
