# Standard Resampling and Inference for Muscal Fits

Runs bootstrap or permutation inference for \`muscal\` fits that expose
a standard refit contract through \`fit_spec\$refit\`. This keeps
resampling logic generic across method families while still allowing
each method to define its own data reconstruction and null-generation
strategy.

## Usage

``` r
infer_muscal(
  object,
  method = c("bootstrap", "permutation"),
  statistic = c("sdev", "loadings"),
  statistic_fn = NULL,
  nrep = 100,
  alpha = 0.05,
  seed = NULL,
  alternative = c("greater", "less", "two.sided"),
  refit = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted \`muscal\` model.

- method:

  One of \`"bootstrap"\` or \`"permutation"\`.

- statistic:

  One of \`"sdev"\` or \`"loadings"\`. Ignored if \`statistic_fn\` is
  supplied.

- statistic_fn:

  Optional extractor \`function(fit)\`, used for custom statistics. For
  permutation inference this must return a numeric vector.

- nrep:

  Number of resampling replicates.

- alpha:

  Tail probability used for bootstrap intervals and permutation
  reference quantiles.

- seed:

  Optional RNG seed for reproducible resampling.

- alternative:

  Alternative hypothesis for permutation p-values.

- refit:

  Optional refit specification overriding \`object\$fit_spec\$refit\`.

- verbose:

  Logical; if \`TRUE\`, emit progress and failure counts.

- ...:

  Reserved for future extensions.

## Value

A list inheriting from \`muscal_inference_result\` and either
\`muscal_bootstrap_result\` or \`muscal_permutation_result\`.

## Details

The default statistics currently supported are: - \`"sdev"\`:
component-wise standard deviations / singular values - \`"loadings"\`:
the concatenated loading matrix \`fit\$v\` (bootstrap only)

Methods without stored refit metadata will error unless an explicit
\`refit\` specification is supplied.
