# Multifer Inference for BaDA

Run \`multifer\` inference for a fitted BaDA model using subject-level
class-barycenter blocks as the multiblock data representation.

## Usage

``` r
infer_bada(
  x,
  data,
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

  A fitted \`bada\` object.

- data:

  The original \`multidesign\` data used to fit \`x\`.

- targets:

  Inference targets passed to \[multifer::infer()\].

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

A \`multifer::infer_result\` object.
