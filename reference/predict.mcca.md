# Predict from an MCCA Fit

Predict from an MCCA Fit

## Usage

``` r
# S3 method for class 'mcca'
predict(object, new_data = NULL, type = c("scores", "reconstruction"), ...)
```

## Arguments

- object:

  A fitted \`mcca\` object.

- new_data:

  Optional matrix/data.frame or list of blocks with the same structure
  as the training data.

- type:

  One of \`"scores"\` or \`"reconstruction"\`.

- ...:

  Additional arguments passed to the underlying projection helpers.

## Value

A numeric matrix of projected scores or reconstructed observations.
