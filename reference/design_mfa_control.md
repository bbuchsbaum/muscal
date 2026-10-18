# Control Numerical Design-MFA Calculations

Control Numerical Design-MFA Calculations

## Usage

``` r
design_mfa_control(
  rank_tol = sqrt(.Machine$double.eps),
  estimability_tol = 1e-08,
  max_numeric_grid = 20L,
  max_grid_rows = 5000L,
  feature_chunk_size = 25000L,
  multilevel_max_iter = 5000L,
  multilevel_tol = 1e-08,
  soft_max_iter = 500L,
  soft_tol = 1e-08
)
```

## Arguments

- rank_tol:

  Relative eigenvalue threshold for numerical design rank and
  effect-component rank.

- estimability_tol:

  Relative tolerance for reference-grid estimability.

- max_numeric_grid:

  Maximum number of distinct numeric values included automatically in a
  reference grid. Variables with more distinct values use this many
  bounded sample quantiles.

- max_grid_rows:

  Maximum allowed reference-grid size.

- feature_chunk_size:

  Maximum number of response features processed at once while
  accumulating shared-grid sufficient statistics.

- multilevel_max_iter:

  Maximum block-coordinate iterations for ridge subject deviations. The
  default leaves headroom for weak-pooling fits, which can converge
  substantially more slowly than strongly pooled fits.

- multilevel_tol:

  Relative penalized-objective convergence tolerance.

- soft_max_iter:

  Maximum conjugate-gradient iterations for soft alignment.

- soft_tol:

  Relative normal-equation residual tolerance for soft alignment.

## Value

A named control list for \[design_mfa()\].

## Examples

``` r
control <- design_mfa_control(feature_chunk_size = 1000L, soft_tol = 1e-7)
control[c("feature_chunk_size", "soft_tol")]
#> $feature_chunk_size
#> [1] 1000
#> 
#> $soft_tol
#> [1] 1e-07
#> 
```
