# Muscal Fit Object Contract

CV-capable \`muscal\` methods should expose a small common contract so
they can participate in shared cross-validation and inference workflows.

## Usage

``` r
muscal_fit_contract(
  method,
  task,
  oos_types,
  fit_call = NULL,
  refit_supported = FALSE,
  prediction_target = NULL,
  refit = NULL
)
```

## Arguments

- method:

  Character scalar naming the fitting method.

- task:

  Character scalar naming the primary evaluation task.

- oos_types:

  Character vector naming the supported out-of-sample prediction types.

- fit_call:

  Optional matched call captured at fit time.

- refit_supported:

  Logical; whether the object currently stores enough metadata for a
  generic refit path.

- prediction_target:

  Optional character scalar naming the high-level prediction target.

- refit:

  Optional refit specification created by \`.muscal_make_refit_spec()\`.

## Value

A named list storing contract metadata.

## Details

The current contract stores: - \`task\`: the primary evaluation task -
\`oos_types\`: the supported \`predict(type = )\` outputs for
out-of-sample use - \`fit_spec\`: a small metadata record describing the
fitting method and whether generic refitting is currently supported -
\`refit\`: optional metadata and callbacks for standard resampling
workflows

The contract does not guarantee refitting for every method. Methods that
support standard resampling attach a \`refit\` record describing how to
rebuild the fit from stored training data and how to generate bootstrap
or permutation replicates.
