# Summarize Components for a Design-MFA Effect

Summarize Components for a Design-MFA Effect

## Usage

``` r
# S3 method for class 'design_mfa'
components(x, term, ...)
```

## Arguments

- x:

  A fitted \`design_mfa\` object.

- term:

  A fitted design term.

- ...:

  Unused.

## Value

A data frame of singular values and fractions of total weighted
reference-grid effect energy. Fractions need not sum to one when
\`rank\` truncates the effect decomposition.

## Examples

``` r
example_data <- design_mfa_example_data()
example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency)
components(example_fit, "saliency")
#>   component singular_value    explained
#> 1         1   2.6187385119 1.000000e+00
#> 2         2   0.0003246124 1.536547e-08
```
