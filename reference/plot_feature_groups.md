# Plot Feature Group Loadings

For Linked MFA models with feature groups, shows how grouped features
align in the loading space.

## Usage

``` r
plot_feature_groups(x, dims = c(1, 2), groups = NULL, show_centers = TRUE, ...)
```

## Arguments

- x:

  A fitted \`linked_mfa\` object with \`feature_groups\`.

- dims:

  Integer vector of length 2 for component axes.

- groups:

  Which groups to show. \`NULL\` means all groups.

- show_centers:

  Logical; if \`TRUE\`, show group centers.

- ...:

  Additional arguments.

## Value

A ggplot object or base R plot.
