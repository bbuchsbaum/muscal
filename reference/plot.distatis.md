# Plot DISTATIS Object and Assessor Maps

Plot DISTATIS Object and Assessor Maps

## Usage

``` r
# S3 method for class 'distatis'
plot(
  x,
  dims = c(1, 2),
  type = c("objects", "subjects", "roi"),
  labels = FALSE,
  col = NULL,
  pch = 19,
  cex = 1,
  main = NULL,
  ...
)

# S3 method for class 'distatis'
autoplot(
  object,
  dims = c(1, 2),
  type = c("objects", "subjects", "roi"),
  labels = FALSE,
  color = NULL,
  alpha = 0.8,
  size = 2,
  ...
)

# S3 method for class 'distatis'
plot_partial_scores(
  x,
  dims = c(1, 2),
  blocks = NULL,
  connect = TRUE,
  show_consensus = TRUE,
  alpha = 0.6,
  ...
)

# S3 method for class 'distatis'
plot_block_weights(x, ...)

# S3 method for class 'distatis'
plot_block_similarity(x, method = c("RV", "RV2"), ...)
```

## Arguments

- x, object:

  A fitted `distatis` model.

- dims:

  Two distinct retained axes in the chosen map. Each map has its own
  rank.

- type:

  `"objects"` (alias `"roi"`) for the compromise map, or `"subjects"`
  for the assessor principal coordinates in RV space.

- labels:

  Logical; show object or assessor labels?

- col, pch, cex, main:

  Base graphics point color, symbol, size, and title.

- ...:

  Further arguments passed to base graphics, or unused by `autoplot`.

- color, alpha, size:

  ggplot point colors, opacity, and size.

- blocks:

  Optional integer indices of assessors included in the partial map.

- connect, show_consensus:

  Connect partial points to the compromise and show compromise points?

- method:

  Similarity measure. Only `"RV"` is supported for DISTATIS.

## Value

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) invisibly
returns the displayed coordinates; `autoplot()` returns a ggplot. Axis
percentages use the complete spectrum of the chosen map, including axes
not retained in the scores.

## See also

[`distatis()`](https://bbuchsbaum.github.io/muscal/reference/distatis.md),
[`plot_partial_scores()`](https://bbuchsbaum.github.io/muscal/reference/plot_partial_scores.md),
[`plot_block_weights()`](https://bbuchsbaum.github.io/muscal/reference/plot_block_weights.md),
[`plot_block_similarity()`](https://bbuchsbaum.github.io/muscal/reference/plot_block_similarity.md),
[`plot_variance()`](https://bbuchsbaum.github.io/muscal/reference/plot_variance.md)
