# Project New Fields onto a Design-Effect Basis

Project New Fields onto a Design-Effect Basis

## Usage

``` r
# S3 method for class 'design_mfa_effect'
project(x, new_data, ...)
```

## Arguments

- x:

  A \`design_mfa_effect\` object.

- new_data:

  A numeric row-by-feature matrix.

- ...:

  Unused.

## Value

Component scores for the new fields.

## Examples

``` r
example_data <- design_mfa_example_data()
example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency)
project(effects(example_fit, "saliency"), example_data$Y[[1L]][1:3, ])
#>     Component_1 Component_2
#> old -3.59257972   -1.061929
#> old -1.75263447   -1.106770
#> old  0.08139509   -1.095468
```
