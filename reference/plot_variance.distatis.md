# DISTATIS Inertia Plot

DISTATIS Inertia Plot

## Usage

``` r
plot_variance.distatis(
  x,
  type = c("bar", "line"),
  axes = c("objects", "subjects"),
  ...
)
```

## Arguments

- x:

  A fitted `distatis` model.

- type:

  `"bar"` or `"line"` (with cumulative inertia).

- axes:

  Plot `"objects"` or `"subjects"` inertia.

- ...:

  Additional arguments passed to
  [`plot_variance()`](https://bbuchsbaum.github.io/muscal/reference/plot_variance.md).

## Value

A ggplot of the complete positive spectrum, including axes not retained
in the score maps. Call through `plot_variance(fit, axes = "subjects")`.
