# Plot Block Reconstruction Error

Bar chart showing per-block reconstruction error, indicating which
blocks are well-captured by the shared scores.

## Usage

``` r
plot_block_fit(x, normalize = TRUE, ...)
```

## Arguments

- x:

  A fitted \`linked_mfa\` object.

- normalize:

  Logical; if \`TRUE\`, show relative (percentage) error.

- ...:

  Additional arguments.

## Value

A ggplot object or base R plot.
