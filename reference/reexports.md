# Extract Design Effects

These objects are imported from other packages. Follow the links below
to see their documentation.

- multivarious:

  [`components`](https://bbuchsbaum.github.io/multivarious/reference/components.html)

- stats:

  [`effects`](https://rdrr.io/r/stats/effects.html)

## Arguments

- x:

  A fitted \`design_mfa\` object.

- term:

  A fitted design term.

- ...:

  Unused.

## Value

A list of effects, or one \`design_mfa_effect\` when \`term\` is
supplied.

A data frame of singular values and fractions of total weighted
reference-grid effect energy. Fractions need not sum to one when
\`rank\` truncates the effect decomposition.

## Examples

``` r
example_data <- design_mfa_example_data()
example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency + probe)
effects(example_fit, "saliency")
#> Design effect: saliency
#>   Status: estimable
#>   Estimable df: 4
#>   Components: 2
#>   Partial R2: 0.9995
```
