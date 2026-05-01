# Multifer Component Tests for COVSTATIS

Runs adapter-owned \`multifer\` component-significance tests for either
the COVSTATIS interstructure/RV axes across input tables or the
compromise ROI eigencomponents.

## Usage

``` r
infer_covstatis(
  x,
  data,
  axes = c("compromise", "interstructure", "both"),
  B = 1000L,
  alpha = 0.05,
  seed = NULL,
  strict = TRUE,
  ...
)
```

## Arguments

- x:

  A fitted \`covstatis\` object.

- data:

  Original list of covariance/correlation matrices used to fit \`x\`.

- axes:

  Which component family to test: \`"interstructure"\` tests axes of the
  table-by-table RV/interstructure matrix, \`"compromise"\` tests ROI
  eigencomponents of the weighted compromise matrix, and \`"both"\` runs
  both tests and returns a named list.

- B:

  Number of Monte Carlo draws per ladder rung.

- alpha:

  Significance threshold passed to \[multifer::infer()\].

- seed:

  Optional random seed.

- strict:

  Logical; passed to \[multifer::infer()\].

- ...:

  Additional arguments passed to \[multifer::infer()\].

## Value

A \`multifer::infer_result\` object, or a named list of two such objects
when \`axes = "both"\`.
