# Plot Partial Factor Scores

Displays how each block "sees" the observations. For MFA, this shows
block-specific projections. For Linked MFA, connects auxiliary block
observations to their reference positions.

## Usage

``` r
plot_partial_scores(
  x,
  dims = c(1, 2),
  blocks = NULL,
  connect = TRUE,
  show_consensus = TRUE,
  alpha = 0.6,
  ...
)
```

## Arguments

- x:

  A fitted \`mfa\` or \`linked_mfa\` object.

- dims:

  Integer vector of length 2 for component axes.

- blocks:

  Which blocks to include. \`NULL\` means all blocks.

- connect:

  Logical; if \`TRUE\`, draw lines connecting the same observation
  across blocks.

- show_consensus:

  Logical; if \`TRUE\`, show the global scores as larger points.

- alpha:

  Transparency for partial score points.

- ...:

  Additional arguments passed to plotting functions.

## Value

A ggplot object (if ggplot2 available) or base R plot.

## Details

For MFA: Each block's partial scores are computed by projecting that
block alone onto the compromise space. Observations where partial scores
diverge indicate block disagreement.

For Linked MFA the auxiliary block observations are shown at their
projected positions, with lines connecting to the reference score.
