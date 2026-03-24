# Tune Tied Penalties for iPCA via 1D Alpha Search

Performs a lightweight 1D search over \`alpha\` and ties block penalties
as \`lambda_k = alpha / p_k\` (or related schemes), then selects the
value with the lowest held-out entry reconstruction MSE.

## Usage

``` r
ipca_tune_alpha(
  data,
  preproc = multivarious::center(),
  ncomp = 2,
  alpha_grid = 10^seq(-3, 3, by = 1),
  tie = c("inv_p", "inv_pbar", "equal"),
  holdout_frac = 0.05,
  n_masks = 1,
  warm_start = TRUE,
  method = c("auto", "gram", "dense"),
  max_iter = 100,
  tol = 1e-06,
  normalize_trace = TRUE,
  use_future = FALSE,
  eig_solver = c("auto", "full", "truncated"),
  eig_rank = NULL,
  eig_trunc_min_n = 400,
  seed = 1,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  A list of matrices/data.frames or a \`multiblock\` object. Blocks must
  share rows.

- preproc:

  A \`multivarious\` preprocessing pipeline.

- ncomp:

  Integer; number of joint components to evaluate.

- alpha_grid:

  Positive numeric vector of candidate alpha values.

- tie:

  Penalty tying rule: \`"inv_p"\` (\`alpha / p_k\`), \`"inv_pbar"\`
  (\`alpha \* mean(p) / p_k\`), or \`"equal"\` (\`alpha\` for all
  blocks).

- holdout_frac:

  Fraction of entries to hold out in each block.

- n_masks:

  Integer; number of random holdout masks to average over.

- warm_start:

  Logical; if \`TRUE\`, alpha candidates are fit in ascending order and
  each fit is initialized from the previous alpha within a mask.

- method:

  One of \`"auto"\`, \`"gram"\`, or \`"dense"\` for iPCA fits.

- max_iter:

  Maximum number of iPCA iterations per alpha.

- tol:

  Convergence tolerance for iPCA fits.

- normalize_trace:

  Logical; passed to \[ipca()\].

- use_future:

  Logical; passed to \[ipca()\].

- eig_solver:

  Eigensolver policy passed to \[ipca()\].

- eig_rank:

  Optional truncated eigensolver rank passed to \[ipca()\].

- eig_trunc_min_n:

  Minimum \`n\` at which \`eig_solver = "auto"\` can switch to truncated
  eigendecomposition in dense mode.

- seed:

  Integer random seed for holdout mask generation.

- verbose:

  Logical; print per-alpha progress.

- ...:

  Additional arguments passed to \[ipca()\].

## Value

A list with:

- best_alpha:

  Selected alpha value.

- best_lambda:

  Selected lambda vector of length K.

- results:

  Data frame of alpha candidates and held-out MSE.

- fit:

  iPCA fit refit on full data with \`best_alpha\`.
