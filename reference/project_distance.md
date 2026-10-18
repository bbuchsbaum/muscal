# Project a Supplementary Assessor's Distances

Project a Supplementary Assessor's Distances

## Usage

``` r
project_distance(x, new_data, squared = x$squared, ...)

# S3 method for class 'distatis'
project(x, new_data, ...)

# S3 method for class 'distatis'
project_cov(x, new_data, ...)
```

## Arguments

- x:

  A fitted `distatis` model.

- new_data:

  A distance matrix or `dist` object on the fitted objects in their
  original order. This projects a new assessor, not new objects.

- squared:

  Whether new distances are already squared. Defaults to the setting
  used when fitting `x`.

- ...:

  Unused arguments, rejected.

## Value

An object-by-component partial score matrix. The fitted compromise and
assessor weights remain fixed.

## See also

[`project_subjects()`](https://bbuchsbaum.github.io/muscal/reference/project_subjects.md),
[`distatis()`](https://bbuchsbaum.github.io/muscal/reference/distatis.md)
