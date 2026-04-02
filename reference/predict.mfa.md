# Predict from a Multiple Factor Analysis Fit

Predict from a Multiple Factor Analysis Fit

## Usage

``` r
# S3 method for class 'mfa'
predict(object, new_data = NULL, type = c("scores", "reconstruction"), ...)
```

## Arguments

- object:

  A fitted \`mfa\` object.

- new_data:

  Optional matrix/data.frame of new observations in the full
  concatenated feature space.

- type:

  One of \`"scores"\` or \`"reconstruction"\`.

- ...:

  Additional arguments passed to the underlying projector helpers.

## Value

A numeric matrix of projected scores or reconstructed observations.
