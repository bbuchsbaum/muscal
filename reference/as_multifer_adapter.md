# Resolve a muscal fit to a multifer adapter

Returns the \`multifer\` adapter used to expose feature-importance hooks
for fitted \`muscal\` methods. The returned adapter can be passed
directly to \`multifer::feature_importance_from_fit()\`,
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
