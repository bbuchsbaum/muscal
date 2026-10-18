# Tune Design-MFA Components with Held-Out Subjects

Selects component rank, spatial smoothness, fixed partial-pooling
penalties, and soft-alignment stabilization using subject-disjoint
folds. For each fold, a candidate effect is learned from analysis
subjects and compared with the same design-effect matrix estimated
independently from assessment subjects.

## Usage

``` r
tune_design_mfa(
  Y,
  design,
  formula,
  subject = NULL,
  terms = NULL,
  grid = NULL,
  n_folds = 3L,
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

  A one-sided design formula passed to \[design_mfa()\].

- subject:

  Optional unique subject identifiers.

- terms:

  Optional terms scored during tuning. The default scores every
  requested formula term.

- grid:

  Candidate data frame or named list with \`rank\`, \`spatial_penalty\`,
  \`subject_penalty\`, and \`alignment_penalty\` fields.

- n_folds:

  Number of subject-disjoint folds. Every assessment fold must contain
  at least two subjects.

- seed:

  Seed for deterministic named-subject fold construction. The caller's
  global random-number state is restored on successful and error exits.

- refit:

  Refit \[design_mfa()\] on all subjects with the selected values.

- verbose:

  Report failed candidate-fold fits as they occur.

- ...:

  Fixed arguments passed to \[design_mfa()\]. Tuned fields and \`store\`
  are controlled by this function.

## Value

A \`design_mfa_tuning\` object containing candidate summaries,
fold-level scores, retained failure diagnostics, subject membership,
selected parameters, and the full refit when \`refit = TRUE\`.

## Details

The score is the mean, over estimable requested terms, of
reference-weighted squared disagreement divided by the held-out effect
energy. This measures recovery of a separately estimated design effect,
not trial-response prediction. Each analysis fold derives its reference
grid/weights without assessment responses, and that same target
grid/weights is used for its assessment fit. Spatial generalized
components are reconstructed with their dual metric loadings. Every
candidate must succeed on every fold to be eligible. Failed
candidate-fold fits remain visible in \`fold_scores\` and are aggregated
in \`failure_summary\`; an otherwise successful tuning run warns when
any candidate was excluded by a fit failure.

This procedure selects hyperparameters; its minimum score is not an
unbiased estimate of final performance. For outer-subject performance
assessment, run \`tune_design_mfa()\` inside an outer subject-level
resampling loop.

## Examples

``` r
# \donttest{
example_data <- design_mfa_example_data()
tuned <- tune_design_mfa(
  example_data$Y, example_data$design, ~ saliency,
  terms = "saliency", grid = data.frame(rank = c(1L, 2L)),
  n_folds = 3, seed = 11, refit = FALSE
)
tuned$best_params
#> $rank
#> [1] 1
#> 
#> $spatial_penalty
#> [1] 0
#> 
#> $subject_penalty
#> [1] 1
#> 
#> $alignment_penalty
#> [1] 0
#> 
# }
```
