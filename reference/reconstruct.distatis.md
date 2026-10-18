# Reconstruct the Normalized DISTATIS Compromise

Reconstruct the Normalized DISTATIS Compromise

## Usage

``` r
# S3 method for class 'distatis'
reconstruct(
  x,
  comp = seq_len(ncol(x$v)),
  what = c("gram", "distance"),
  squared = FALSE,
  ...
)
```

## Arguments

- x:

  A fitted `distatis` model.

- comp:

  Nonempty vector of distinct retained object component indices.

- what:

  Return the centered `"gram"` matrix (default) or `"distance"`.

- squared:

  For `what = "distance"`, return squared distances?

- ...:

  Unused arguments, rejected.

## Value

A square matrix in normalized compromise units, not the original
assessor-specific distance scales. Distance reconstruction uses the Gram
diagonal identity and optionally its elementwise square root.
