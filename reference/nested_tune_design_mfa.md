# Nested Held-Out-Subject Evaluation for Design-MFA Tuning

Runs \[tune_design_mfa()\] entirely within each outer analysis set and
scores the selected model against an independently estimated design
effect in the outer assessment subjects. The resulting outer scores
evaluate the complete tuning procedure without reusing an outer subject
for candidate selection.

## Usage

``` r
nested_tune_design_mfa(
  Y,
  design,
  formula,
  subject = NULL,
  terms = NULL,
  grid = NULL,
  outer_folds = 3L,
  inner_folds = 3L,
  seed = 1L,
  refit = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- Y:

  A named list of subject trial-by-feature matrices.

- design:

  A matching list of trial-design data frames.

- formula:

  A one-sided design formula.

- subject:

  Optional unique subject identifiers.

- terms:

  Optional terms scored in inner and outer folds.

- grid:

  Candidate grid passed to \[tune_design_mfa()\].

- outer_folds:

  Number of outer subject folds.

- inner_folds:

  Number of inner subject folds. This is reduced when necessary to
  retain at least two assessment subjects per fold.

- seed:

  Seed controlling both named-subject fold levels. The caller's global
  random-number state is restored on return, including error exits.

- refit:

  Run a final tuning pass on all subjects for deployment.

- verbose:

  Report inner fitting failures.

- ...:

  Fixed \[design_mfa()\] arguments passed through both levels.

## Value

A \`nested_design_mfa_tuning\` object with outer scores and folds, every
inner tuning result, and an optional final full-data tuning/refit.

## Examples

``` r
# \donttest{
example_data <- design_mfa_example_data()
nested <- nested_tune_design_mfa(
  example_data$Y, example_data$design, ~ saliency,
  terms = "saliency", grid = data.frame(rank = 1L),
  outer_folds = 3, inner_folds = 2, seed = 12, refit = FALSE
)
nested$outer_results
#>   fold        score n_terms rank spatial_penalty subject_penalty
#> 1    1 2.236171e-05       1    1               0               1
#> 2    2 1.819925e-05       1    1               0               1
#> 3    3 2.193891e-05       1    1               0               1
#>   alignment_penalty
#> 1                 0
#> 2                 0
#> 3                 0
# }
```
