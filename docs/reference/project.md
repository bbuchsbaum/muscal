# Re‑export selected generics from \*\*multivarious\*\*

The package extends the \`project()\` and \`reprocess()\` generics
defined in multivarious. Re‑exporting them avoids forcing users to load
multivarious explicitly.

## Usage

``` r
project(x, new_data, ...)

reprocess(x, new_data, colind, ...)

# S3 method for class 'ipca'
project(x, new_data, ...)
```

## Arguments

- x:

  A fitted model object.

- new_data:

  New data to project or reprocess.

- ...:

  Additional arguments passed to methods.

- colind:

  Optional column indices for subsetting.
