# Metric Registry for Muscal Tasks

Returns the task-aware performance metric registry used by
\[performance_metrics()\]. The registry records which metrics belong to
which evaluation task, whether larger values are better, and a short
description of the metric.

## Usage

``` r
metric_registry(task = NULL)
```

## Arguments

- task:

  Optional task name. If supplied, returns only the metrics for that
  task (after alias normalization).

## Value

A tibble with one row per supported metric.

## Details

Supported task families currently include: - \`"reconstruction"\` -
\`"response_prediction"\` - \`"retrieval_alignment"\`

The alias \`"row_alignment"\` maps to \`"retrieval_alignment"\`.
