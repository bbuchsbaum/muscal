# Compute Task-Aware Performance Metrics

Computes a one-row tibble of performance metrics for a declared
evaluation task. This provides a stable metric interface for fold-wise
cross-validation and held-out evaluation.

## Usage

``` r
performance_metrics(
  task,
  truth,
  estimate,
  metrics = NULL,
  k = 5L,
  truth_ids = NULL,
  by_column = FALSE,
  ...
)
```

## Arguments

- task:

  Task name.

- truth:

  Truth object for the task. For reconstruction and response prediction
  this is a numeric matrix/data.frame. For retrieval/alignment this is a
  numeric matrix/data.frame of query-side feature vectors.

- estimate:

  Estimated or predicted object for the task.

- metrics:

  Optional character vector of metrics. If \`NULL\`, task-specific
  defaults are used.

- k:

  Integer \`k\` used by top-k retrieval metrics and recall@k.

- truth_ids:

  Optional vector of true ids for retrieval metrics such as
  \`recall_at_k\` and \`mrr\`.

- by_column:

  Logical; forwarded to reconstruction-style \`r2\` calculations.

- ...:

  Reserved for future task-specific options.

## Value

A one-row tibble of metric values.

## Details

Supported task families currently include: - \`"reconstruction"\` -
\`"response_prediction"\` - \`"retrieval_alignment"\`

The alias \`"row_alignment"\` maps to \`"retrieval_alignment"\`.

For \`"retrieval_alignment"\`, \`estimate\` may be: - a numeric
matrix/data.frame interpreted as a single top-1 retrieved feature vector
per query, or - a list containing \`retrieved_features\` (a list of
ranked feature matrices), plus optional \`retrieved_ids\`,
\`oracle_similarity\`, or \`oracle_features\`.
