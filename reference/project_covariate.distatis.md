# Project Covariates onto DISTATIS Maps

Project Covariates onto DISTATIS Maps

## Usage

``` r
# S3 method for class 'distatis'
project_covariate(
  x,
  y,
  side = c("table", "feature"),
  what = "dimension",
  scale = c("cosine", "beta"),
  ...
)
```

## Arguments

- x:

  A fitted `distatis` model.

- y:

  Finite numeric covariate, one value per assessor or object.

- side:

  `"table"` for assessors or `"feature"` for objects.

- what:

  `"dimension"` returns uncentered columnwise cosines or regression
  slopes in the selected map. For table covariates, `"feature"` (alias
  `"observation"`) returns an object-by-component pattern by projecting
  the partial scores across assessors onto `y`.

- scale:

  `"cosine"` or `"beta"`. For dimension projections the slope regresses
  `y` on each map axis. For feature patterns it regresses each
  partial-score coordinate on `y`. No intercept or automatic centering
  is used.

- ...:

  Unused arguments, rejected.

## Value

A named vector for dimension projections, otherwise an object score
matrix. Assessor and object projections may have different axis counts.
