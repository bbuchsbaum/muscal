# Plot Block Similarity Matrix

Heatmap of pairwise RV coefficients between blocks, showing which data
sources tell similar stories.

## Usage

``` r
plot_block_similarity(x, method = c("RV", "RV2"), ...)
```

## Arguments

- x:

  A fitted \`mfa\` object.

- method:

  Similarity method: \`"RV"\` or \`"RV2"\`.

- ...:

  Additional arguments.

## Value

A ggplot object or base R image plot.
