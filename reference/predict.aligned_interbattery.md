# Predict from an Aligned Interbattery Fit

Predict from an Aligned Interbattery Fit

## Usage

``` r
# S3 method for class 'aligned_interbattery'
predict(
  object,
  new_data,
  block = NULL,
  from = c("x", "y", "both"),
  type = c("prediction", "reconstruction", "common_scores", "side_scores"),
  target_block = NULL,
  preprocess = TRUE,
  x_row_map_new = NULL,
  y_row_map_new = NULL,
  row_graph_new = NULL,
  row_graph_form = NULL,
  ...
)
```

## Arguments

- object:

  A fitted \`aligned_interbattery\` object.

- new_data:

  Observed out-of-sample data. For \`from = "x"\` or \`"y"\`, this may
  be a single matrix/data.frame or a list of observed blocks from that
  side. For \`from = "both"\`, supply a list with optional elements
  \`x\` and \`y\`, each containing a single block or block bundle for
  that side.

- block:

  Optional source block identifier(s). For single-matrix input, this
  identifies the source loading type or training block. For unnamed
  block bundles, provide one identifier per observed block. For \`from =
  "both"\`, \`block\` may be \`NULL\` or a list with optional
  \`x\`/\`y\` entries matching the supplied side bundles.

- from:

  One of \`"x"\`, \`"y"\`, or \`"both"\`.

- type:

  One of \`"prediction"\`, \`"reconstruction"\`, \`"common_scores"\`, or
  \`"side_scores"\`.

- target_block:

  Optional target loading type or training block name. For one-sided
  prediction this targets the opposite side. For \`from = "both"\`, it
  may identify either side; ambiguous names can be disambiguated via
  \`x:\<name\>\` or \`y:\<name\>\`. When omitted, multi-target calls
  return named lists.

- preprocess:

  Logical; if \`TRUE\` (default), apply the fitted shared preprocessing
  pipeline(s) before projection.

- x_row_map_new:

  Optional integer row-map structure for new \`x\`-side blocks. Missing
  maps imply identity mapping and therefore require aligned row counts
  within the supplied \`x\` bundle.

- y_row_map_new:

  Optional integer row-map structure for new \`y\`-side blocks.

- row_graph_new:

  Optional graph/Laplacian on the new rows used to refine the
  common-score solve. This is never required.

- row_graph_form:

  Interpretation of \`row_graph_new\`. When \`NULL\` (default), the
  fitted object's \`row_graph_form\` is used.

- ...:

  Unused.

## Value

A numeric matrix or a named list of matrices for multi-target prediction
/ completion.
