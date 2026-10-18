# Extract Subject-Specific Design Effects

Extract Subject-Specific Design Effects

## Usage

``` r
subject_effects(x, term = NULL)
```

## Arguments

- x:

  A fitted \`design_mfa\` object.

- term:

  Optional term.

## Value

Subject deviations when fitted by a partial-pooling model.

## Examples

``` r
example_data <- design_mfa_example_data()
pooled <- design_mfa(
  example_data$Y, example_data$design, ~ saliency,
  subject_model = "ridge_slopes", random_terms = "saliency"
)
lapply(subject_effects(pooled, "saliency"), dim)
#> $S1
#> [1] 4 8
#> 
#> $S2
#> [1] 4 8
#> 
#> $S3
#> [1] 4 8
#> 
#> $S4
#> [1] 4 8
#> 
#> $S5
#> [1] 4 8
#> 
#> $S6
#> [1] 4 8
#> 
```
