# Plot a Design-Effect Component Trajectory

Plot a Design-Effect Component Trajectory

## Usage

``` r
# S3 method for class 'design_mfa'
plot(x, term = NULL, component = 1L, xvar = NULL, ...)
```

## Arguments

- x:

  A fitted \`design_mfa\` object.

- term:

  Optional fitted design term. By default, the first estimable effect
  with a retained component is plotted.

- component:

  Component index.

- xvar:

  Optional reference-grid variable for the horizontal axis.

- ...:

  Additional arguments passed to \[graphics::plot()\].

## Value

The plotted data, invisibly.

## Examples

``` r
example_data <- design_mfa_example_data()
example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency)
plot(example_fit)
```
