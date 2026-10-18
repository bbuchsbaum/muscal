# Extract Design Effects from a Design-MFA Fit

Extract Design Effects from a Design-MFA Fit

## Usage

``` r
# S3 method for class 'design_mfa'
effects(object, term = NULL, ...)
```

## Arguments

- object:

  A fitted \`design_mfa\` object.

- term:

  Optional design term.

- ...:

  Unused.

## Value

A list of effects, or one \`design_mfa_effect\` when \`term\` is
supplied.

## Examples

``` r
example_data <- design_mfa_example_data()
example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency + probe)
effects(example_fit, "probe")
#> Design effect: probe
#>   Status: estimable
#>   Estimable df: 2
#>   Components: 1
#>   Partial R2: 0.9952
```
