# Multifer Inference for Aligned RRR

Run \`multifer\` inference for a fitted \[aligned_rrr()\] model. The
adapter uses multifer's opaque predictive geometry so the model can keep
its native multiblock predictor/response data shape.

## Usage

``` r
infer_aligned_rrr(
  x,
  data = NULL,
  targets = "default",
  B = 1000L,
  R = 500L,
  alpha = 0.05,
  seed = NULL,
  strict = TRUE,
  ...
)
```

## Arguments

- x:

  A fitted \`aligned_rrr\` object.

- data:

  Optional original data list with elements \`X\` and \`Y\`. When
  \`NULL\`, the data stored in \`x\$fit_spec\$refit\$data\` are used.

- targets:

  Inference targets passed to \[multifer::infer()\]. Defaults to all
  targets declared by the adapter: component significance plus loading,
  score, and subspace stability.

- B:

  Number of Monte Carlo draws per component-significance rung.

- R:

  Number of bootstrap replicates for stability targets.

- alpha:

  Significance level passed to \[multifer::infer()\].

- seed:

  Optional random seed.

- strict:

  Logical; passed to \[multifer::infer()\].

- ...:

  Additional arguments passed to \[multifer::infer()\].

## Value

A \`multifer::infer_result\` object, or a \`multifer\` bundle when
\`return_bundle = TRUE\` is passed through \`...\`.
