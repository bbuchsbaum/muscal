# Convert Sorting Partitions to Distance Matrices

Convert Sorting Partitions to Distance Matrices

## Usage

``` r
sorting_distances(data)
```

## Arguments

- data:

  An object-by-assessor matrix or data frame of categorical group
  identifiers. Identifiers may be numeric, character, logical, or
  factors. Every object must have exactly one nonmissing group
  identifier per assessor.

## Value

A named list of binary distance matrices. A pair has distance zero when
the assessor grouped it together and one otherwise. Row/column names
identify objects. Renaming group identifiers leaves the distances
unchanged.

## See also

[`distatis()`](https://bbuchsbaum.github.io/muscal/reference/distatis.md)
