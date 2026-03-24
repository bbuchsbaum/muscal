# Plot Convergence Trace

Shows the objective function value over iterations for iterative
methods.

## Usage

``` r
plot_convergence(x, log_scale = FALSE, ...)

# S3 method for class 'linked_mfa'
plot_convergence(x, log_scale = FALSE, ...)

# S3 method for class 'bamfa'
plot_convergence(x, log_scale = FALSE, ...)

# S3 method for class 'penalized_mfa'
plot_convergence(x, log_scale = FALSE, ...)
```

## Arguments

- x:

  A fitted object with an \`objective_trace\` (e.g., \`linked_mfa\`,
  \`bamfa\`, \`penalized_mfa\`).

- log_scale:

  Logical; if \`TRUE\`, use log scale for y-axis.

- ...:

  Additional arguments.

## Value

A ggplot object or base R plot.
