# Repeated-measures COVSTATIS

Estimate the common covariance geometry and the covariance patterns that
separate conditions in a balanced, one-factor repeated-measures design.

## Usage

``` r
covstatis_rm(
  data,
  subject,
  condition,
  ncomp = 2L,
  ncomp_condition = NULL,
  norm_method = c("none", "frobenius", "mfa"),
  dcenter = FALSE,
  weights = NULL
)

# S3 method for class 'covstatis_rm'
project_subjects(x, new_data, subject, condition, ...)
```

## Arguments

- data:

  List of covariance matrices, one per subject-condition cell. All
  matrices must be finite, symmetric, and describe the same variables in
  the same order. Supply consistent variable names throughout or use
  positional matching throughout. Positive semidefiniteness is assumed.

- subject, condition:

  Vectors identifying the subject and condition of each matrix, in input
  order. Each subject must have one matrix per condition. Fitting
  requires at least two subjects and two conditions. Observed factor
  levels specify condition order. Other labels use order of first
  appearance.

- ncomp:

  Maximum number of common variable axes to retain.

- ncomp_condition:

  Maximum number of condition axes to retain. The default keeps all
  numerically supported axes, up to one fewer than the number of
  conditions.

- norm_method:

  Per-cell scaling. The default, `"none"`, preserves covariance
  magnitude. `"frobenius"` divides each matrix by its Frobenius norm and
  `"mfa"` divides by its largest eigenvalue.

- dcenter:

  Apply Gower double centering to each cell before scaling. The default,
  `FALSE`, preserves covariance along the all-ones direction.

- weights:

  Optional vector of nonnegative finite variable weights, with at least
  one positive entry. Unnamed weights follow matrix order. Named weights
  are matched to variable names. Each covariance entry (i, j) is
  multiplied by sqrt(w\[i\] \* w\[j\]) before double centering and
  normalization. The fitted model stores the weights and applies them to
  raw projection matrices in the training variable order. The default
  gives unit weights.

- x:

  A fitted \`covstatis_rm\` model.

- new_data:

  Covariance matrices for one or more new subjects, with all training
  conditions observed once per subject. Supply \`subject\` and
  \`condition\` as for the training data.

- ...:

  Reserved for methods.

## Value

A `covstatis_rm` list with the following fields:

- `common` is the fitted COVSTATIS object on subject means.

- `weights` contains the variable weights used for raw cells at fit and
  projection.

- `condition_effects` contains one signed mean deviation matrix per
  condition.

- `condition_maps` contains the retained Frobenius-orthonormal loading
  matrices.

- `condition_scores` gives condition centroids on the retained axes.

- `paired_scores` gives subject, condition, and component scores in
  input order.

- `explained` gives the fraction of total condition-effect energy on
  each retained axis. `sums_of_squares` partitions total variation about
  the equal-cell grand mean into subject, condition, and residual terms.
  The residual includes subject-by-condition variation and estimation
  error.

A zero condition effect yields an empty list of maps and zero score
columns.
[`project_subjects()`](https://bbuchsbaum.github.io/muscal/reference/project_subjects.md)
applies the fitted maps to complete panels of new subjects.

For \`project_subjects()\`, a list with \`paired_scores\` in input order
and \`common_scores\`, one partial variable map per new subject.

## Details

The common map is an ordinary
[`covstatis()`](https://bbuchsbaum.github.io/muscal/reference/covstatis.md)
fit on each subject's average preprocessed covariance, with no further
preprocessing. Variable weights apply to every cell in both analyses.
The estimated STATIS table weights apply only to the common map.

For the condition analysis, subtract each subject's mean matrix from
their cells, then average these paired deviations across subjects with
equal weights. Decompose the resulting condition means using the
Frobenius inner product. For each unit-norm loading matrix, the cell
score is its inner product with that cell's paired deviation. Component
signs are arbitrary. A tied set of components identifies a subspace up
to rotation.

For two conditions, subtract the first fitted condition effect from the
second to obtain the mean paired difference in the preprocessed units.
Eigenvectors of this signed difference give variable combinations along
which projected variance increased or decreased.

## References

Madssen TS et al. (2021). Repeated measures ASCA+ for analysis of
longitudinal intervention studies with multivariate outcome data.
[doi:10.1371/journal.pcbi.1009585](https://doi.org/10.1371/journal.pcbi.1009585)
.

## See also

[`covstatis()`](https://bbuchsbaum.github.io/muscal/reference/covstatis.md),
[`design_mfa()`](https://bbuchsbaum.github.io/muscal/reference/design_mfa.md)

## Examples

``` r
subject <- rep(c("s1", "s2", "s3"), each = 2)
condition <- rep(c("rest", "task"), 3)
cells <- lapply(seq_along(subject), function(i) {
  S <- diag(c(5 + match(subject[i], unique(subject)), 2, 2))
  S[2, 3] <- S[3, 2] <- if (condition[i] == "task") 0.5 else 0
  S
})
fit <- covstatis_rm(cells, subject, condition)
fit$common$roi_scores
#>          [,1]    [,2]
#> Obs_1 2.66117 0.00000
#> Obs_2 0.00000 1.06066
#> Obs_3 0.00000 1.06066
fit$paired_scores
#>   subject condition Component_1
#> 1      s1      rest  -0.3535534
#> 2      s1      task   0.3535534
#> 3      s2      rest  -0.3535534
#> 4      s2      task   0.3535534
#> 5      s3      rest  -0.3535534
#> 6      s3      task   0.3535534
delta <- fit$condition_effects$task - fit$condition_effects$rest
eigen(delta, symmetric = TRUE)$values
#> [1]  0.5  0.0 -0.5
```
