# Predict method for BaMFA objects

Reconstructs data blocks from a fitted BaMFA model.

## Usage

``` r
# S3 method for class 'bamfa'
predict(object, new_data = NULL, ...)
```

## Arguments

- object:

  A fitted \`bamfa\` object.

- new_data:

  Optional new data to reconstruct. If NULL (default), reconstructs the
  training data.

- ...:

  Additional arguments (unused).

## Value

A list of reconstructed data matrices.
