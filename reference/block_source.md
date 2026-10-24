# Out-of-Core Block Sources

A block source describes a list of data blocks without holding them in
memory. Fitting functions that accept a block source (currently
\[aligned_mfa()\]) read it one column chunk of one block at a time, so
peak memory is one chunk plus the method's small summary statistics.

## Usage

``` r
block_source(
  loader,
  nrow,
  ncol,
  names = NULL,
  colnames = NULL,
  chunk_size = 5000L
)

as_block_source(x, ...)

# S3 method for class 'muscal_block_source'
as_block_source(x, ...)

# S3 method for class 'list'
as_block_source(x, ...)

# Default S3 method
as_block_source(x, ...)

read_block(x, k, cols = NULL)
```

## Arguments

- loader:

  Function \`function(k, cols)\` returning block \`k\` restricted to the
  integer column indices \`cols\`, as a numeric matrix with
  \`nrow\[k\]\` rows and \`length(cols)\` columns. Every call for block
  \`k\` must return the rows in the same order. Missing and infinite
  values are not supported.

- nrow, ncol:

  Integer vectors giving the number of rows and columns of each block.

- names:

  Optional block names; defaults to \`names(nrow)\` or \`X1, X2, ...\`.

- colnames:

  Optional list of column-name vectors, one per block (each \`NULL\` or
  of length \`ncol\[k\]\`). Needed for \`feature_groups = "colnames"\`
  or for feature groups that refer to features by name.

- chunk_size:

  Maximum number of columns read per call to \`loader\`. Peak memory per
  read is about \`8 \* nrow\[k\] \* chunk_size\` bytes.

- x:

  An object to convert: a list of matrices/data frames, or an existing
  block source.

- ...:

  Passed to \`block_source()\` (e.g. \`chunk_size\`).

- k:

  Block index or name.

- cols:

  Integer column indices to read; defaults to all columns.

## Value

An object of class \`muscal_block_source\`.

For \`read_block()\`, the requested numeric matrix.

## Details

\`block_source()\` wraps a loader function, which is the general
mechanism: any storage (HDF5, file-backed matrices, NIfTI volumes, a
database) can be exposed by writing \`loader(k, cols)\`.
\`as_block_source()\` converts an in-memory list of matrices, which is
mainly useful for testing the out-of-core code paths against the
in-memory ones.

## Examples

``` r
X <- list(A = matrix(rnorm(40), 10), B = matrix(rnorm(60), 10))
src <- block_source(
  loader = function(k, cols) X[[k]][, cols, drop = FALSE],
  nrow = c(10, 10), ncol = c(4, 6), names = c("A", "B"), chunk_size = 3
)
src
#> <muscal_block_source> 2 blocks, chunk_size = 3
#>  block nrow ncol
#>      A   10    4
#>      B   10    6
```
