# Resolve a muscal fit to a multifer adapter

Returns the \`multifer\` adapter for a fitted \`muscal\` method.
Design-MFA exposes \`effect_significance\` with a caller-owned null.
Other methods currently expose feature-importance and
component-significance hooks that can be passed to
\`multifer::feature_importance_from_fit()\`,
\`multifer::check_feature_importance_adapter()\`, and, when the fit/data
expose refitting hooks, \`multifer::feature_importance_pvalues()\`.

## Usage

``` r
as_multifer_adapter(x, ...)
```

## Arguments

- x:

  A fitted \`muscal\` model.

- ...:

  Additional arguments passed to methods.

## Value

A \`multifer_adapter\`.
