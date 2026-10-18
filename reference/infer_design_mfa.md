# Subject-Level Inference for a Design-MFA Effect

Resamples subjects as intact clusters for bootstrap uncertainty.
Permutation inference is available only when the caller supplies an
explicit \`exchangeability\` callback; there is deliberately no generic
trial shuffle.

Bootstrap loading bases are Procrustes-aligned to the fitted reference
basis. Nearly tied components should be interpreted through the returned
subspace distances rather than voxelwise loading intervals.

## Usage

``` r
infer_design_mfa(
  x,
  term,
  method = c("bootstrap", "permutation"),
  statistics = NULL,
  nrep = 200L,
  alpha = 0.05,
  seed = NULL,
  exchangeability = NULL,
  verbose = FALSE,
  use_future = FALSE,
  ...
)
```

## Arguments

- x:

  A fitted \[design_mfa()\] object.

- term:

  One estimable design term.

- method:

  \`"bootstrap"\` or \`"permutation"\`.

- statistics:

  Summaries to compute. The method-dependent default is
  \`c("partial_ss", "singular_values", "subspace")\` for the bootstrap
  and \`c("partial_ss", "singular_values")\` for permutation inference.
  Loading intervals are optional because they can be large.

- nrep:

  Number of resampling replicates.

- alpha:

  Total error probability for central percentile intervals.

- seed:

  Optional random seed. The caller's global random-number state is
  restored on successful and error exits.

- exchangeability:

  For permutation, a function accepting \`data\`, \`term\`, and
  \`replicate\` and returning a complete refit-data list. The callback
  owns the scientific validity of the null transformation. Return
  aligned named \`Y\`, \`design\`, \`trial_weights\`, \`block_weight\`,
  \`reference_grid\`, and \`reference_weight\` fields (plus
  \`feature_map\` for soft alignment). Preserve the frozen reference
  grid and weights. Malformed refits are failures.

- verbose:

  Emit per-replicate progress. Failed replicates always warn.

- use_future:

  If \`TRUE\` and the multifer path is unavailable, run independent
  native refits with \[furrr::future_map()\]. Multifer permutation and
  bootstrap stay sequential.

- ...:

  Unused.

## Value

A \`design_mfa_inference\` object containing observed statistics,
replicate values, failures, and bootstrap summaries or a permutation
p-value.

## Details

The target is conditional on the fitted reference grid, its cell
weights, feature maps, contrasts and chosen penalties/rank. Subject
resampling keeps each subject's trials, design and trial weights
together and preserves the fitted reference grid/weights; it does not
estimate a new population grid. Exchangeability of subjects or a
caller-specified null transformation is a scientific assumption, not
established by callback schema validation. No unconditional
rejection-rate or interval-coverage calibration is claimed.

Bootstrap percentile intervals on \`partial_ss\` (and other nonnegative
statistics) have a positive null expectation and must not be read as a
hypothesis test. Use permutation with an explicit \`exchangeability\`
callback to test a term.

Permutation p-values use successful replicates only: \`(1 + b) /
(n_success + 1)\`. Failed refits are not treated as null draws. If
\`n_success \< nrep\`, the resolution is \`1 / (n_success + 1)\` and a
warning is issued.

When a multifer build with \`effect_significance\` is available, serial
permutation calls \`multifer::infer()\`. When the statistic-stability
sidecar is also available, bootstrap calls
\`multifer::bootstrap_fits()\` and maps percentile intervals from that
sidecar. \`use_future = TRUE\` no longer starts a local \`furrr\` loop
on those paths.

## Examples

``` r
# \donttest{
example_data <- design_mfa_example_data()
example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency)
uncertainty <- infer_design_mfa(
  example_fit, "saliency", statistics = "partial_ss",
  nrep = 5, seed = 10
)
uncertainty$summary$partial_ss
#>   component      label observed    mean          sd    lower    upper
#> 1         1 partial_ss 6.857792 6.85631 0.003164841 6.852106 6.859263
# }
```
