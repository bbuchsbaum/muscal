# Report Design-MFA Estimability

Report Design-MFA Estimability

## Usage

``` r
estimability(x, ...)

# S3 method for class 'design_mfa'
estimability(x, ...)
```

## Arguments

- x:

  A fitted object.

- ...:

  Unused.

## Value

A data frame with term-level rank and reference-grid diagnostics.
\`effect_nonestimable_grid_rows\` concerns the tested term contrast,
whereas \`prediction_nonestimable_grid_rows\` concerns the complete
fitted response, always against the global model and independently of
\`ss_type\`. These can differ when an unrelated between-subject
covariate is aliased with subject fixed effects.

## Examples

``` r
example_data <- design_mfa_example_data()
ages <- seq(20, 70, length.out = length(example_data$design))
adjusted_design <- Map(function(d, age) transform(d, age = age),
                       example_data$design, ages)
names(adjusted_design) <- names(example_data$design)
adjusted <- design_mfa(
  example_data$Y, adjusted_design, ~ saliency + age,
  contrasts = list(saliency = "poly"), terms = "saliency"
)
estimability(adjusted)
#>       term    status estimable_df grid_rows effect_nonestimable_grid_rows
#> 1 saliency estimable            4        30                             0
#>   prediction_nonestimable_grid_rows estimability_scope partial_ss partial_r2
#> 1                                30        effect_only   6.857792  0.9009939
#>   objective_delta objective_fraction  comparison component_rank
#> 1        6.857792          0.9009939 partial_sse              2
```
