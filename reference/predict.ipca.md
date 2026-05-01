# Predict from an iPCA Fit

Predict from an iPCA Fit

## Usage

``` r
# S3 method for class 'ipca'
predict(object, new_data = NULL, type = c("scores", "reconstruction"), ...)
```

## Arguments

- object:

  A fitted \`ipca\` object.

- new_data:

  Optional matrix/data.frame or list of blocks with the same structure
  as the training data.

- type:

  One of \`"scores"\` or \`"reconstruction"\`.

- ...:

  Additional arguments passed to the underlying projection helpers.

## Value

A numeric matrix of projected scores or reconstructed observations.
