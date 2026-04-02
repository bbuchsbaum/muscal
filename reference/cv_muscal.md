# Cross-Validate Muscal Models with Multidesign Folds

A thin cross-validation wrapper that reuses \`multidesign\` fold objects
and the \`muscal\` fit contract. Users provide: - \`fit_fn(analysis)\`
to fit a model on the analysis split, - \`estimate_fn(model, assessment,
...)\` to produce predictions or held-out estimates, and -
\`truth_fn(assessment, ...)\` to extract the corresponding truth object.

## Usage

``` r
cv_muscal(
  folds,
  fit_fn,
  estimate_fn,
  truth_fn,
  task = NULL,
  metrics = NULL,
  performance_args = list(),
  ...
)
```

## Arguments

- folds:

  A \`foldlist\` object, typically created by
  \`multidesign::fold_over()\` or \`multidesign::cv_rows()\`.

- fit_fn:

  Function with signature \`fit_fn(analysis)\` returning a fitted
  \`muscal\` model.

- estimate_fn:

  Function that produces held-out estimates from \`estimate_fn(model,
  assessment, ...)\`.

- truth_fn:

  Function that extracts the truth object from \`truth_fn(assessment,
  ...)\`.

- task:

  Optional task override. If \`NULL\`, the task is taken from
  \`model\$task\`.

- metrics:

  Optional metric vector passed to \[performance_metrics()\].

- performance_args:

  Optional named list of additional arguments forwarded to
  \[performance_metrics()\] on every fold.

- ...:

  Reserved for future extensions.

## Value

An object inheriting from \`cv_result\` with additional fields
\`folds\`, \`foldframe\`, \`task\`, \`metrics\`, and \`call\`.

## Details

The returned metrics are computed through \[performance_metrics()\], so
fold-wise outputs have a stable, task-aware tibble shape that works with
\`multidesign::cross_validate()\`.

This function deliberately does not define a new split abstraction. Use
\`multidesign::fold_over()\` or \`multidesign::cv_rows()\` to construct
\`folds\`.
