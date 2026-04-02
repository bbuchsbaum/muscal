# Barycentric Discriminant Analysis (generic)

A generic function for fitting a \*\*Barycentric Discriminant
Analysis\*\* (BaDA) model to multi–subject multivariate data. Concrete
methods (e.g. for objects of class \`multidesign\`) implement the actual
estimation procedure and return an object that inherits from class
\``bada`\`.

Perform bootstrap resampling on a multivariate bada model to estimate
the variability of components and scores.

This function transforms new data using the preprocessing pipeline from
a fitted BaDA model. It handles different preprocessing scenarios based
on the provided parameters.

This function projects new data onto a previously fitted BaDA model,
returning the scores of the new data in the space of the original model.

## Usage

``` r
bada(data, y, subject, preproc, ncomp, ...)

# S3 method for class 'bada'
bootstrap(x, nboot = 500, ...)

# S3 method for class 'multidesign'
bada(
  data,
  y,
  subject,
  preproc = multivarious::center(),
  ncomp = 2,
  resdim = 20,
  rescomp = 0,
  ...
)

# S3 method for class 'bada'
reprocess(x, new_data, colind = NULL, ...)

# S3 method for class 'bada'
project(x, new_data, ...)
```

## Arguments

- data:

  A data object for which a BaDA method is defined (see individual
  methods for details).

- y:

  A factor or bare column name giving the grouping variable.

- subject:

  A factor or bare column name identifying repeated‑measures, i.e. the
  subject/block variable.

- preproc:

  A preprocessing pipeline from the multivarious package. Defaults to
  \`multivarious::center()\` in the concrete methods.

- ncomp:

  Integer; the number of discriminant components to compute.

- ...:

  Additional arguments:

  block

  :   An optional character string specifying which block's
      preprocessing to apply before projection. If not provided, the
      generic method is used.

- x:

  A fitted BaDA model object.

- nboot:

  An integer specifying the number of bootstrap resamples to perform
  (default is 500).

- resdim:

  pca dimensionality for residual analysis (only relevant if \`rescomp\`
  \> 0)

- rescomp:

  number of final residual components (default = 0, no residual
  aanalysis)

- new_data:

  A numeric matrix of new data to be projected.

- colind:

  An optional integer vector specifying column indices to use within
  blocks. If NULL and block is also NULL, all blocks are used. If NULL
  but block is provided, all columns in the specified block are used.

## Value

An object that inherits from class \`bada\`. See the help for the
corresponding method (e.g. \[`bada.multidesign`\]) for the exact
structure.

A preprocessed numeric matrix with the same number of rows as the input
data.

A numeric matrix of projected scores.

## Details

The function returns a list containing the summarized bootstrap
resampled components and scores for the model. The returned list
contains eight elements: \* \`boot_scores_mean\`: A matrix which is the
mean of all bootstrapped scores matrices. \* \`boot_scores_sd\`: A
matrix which is the standard deviation of all bootstrapped scores
matrices. \* \`boot_scores_upper\`: A matrix which is the upper alpha
percentile of all bootstrapped scores matrices. \*
\`boot_scores_lower\`: A matrix which is the lower alpha percentile of
all bootstrapped scores matrices. \* \`boot_lds_mean\`: A matrix which
is the mean of all bootstrapped loadings matrices. \* \`boot_lds_sd\`: A
matrix which is the standard deviation of all bootstrapped loadings
matrices. \* \`boot_lds_upper\`: A matrix which is the upper alpha
percentile of all bootstrapped loadings matrices. \* \`boot_lds_lower\`:
A matrix which is the lower alpha percentile of all bootstrapped
loadings matrices. The dimensions of each matrix in the list correspond
to the dimensions of the respective matrices in the input data.

The function handles three scenarios: 1. When both colind and block are
NULL: The function applies preprocessing for each block and averages the
results. 2. When block is provided: The function applies preprocessing
specific to the named block. If colind is also provided, it acts as a
relative subset within the block. 3. When only colind is provided: The
function applies preprocessing for each block using the specified column
indices and averages the results.

When a specific block is provided, the function first reprocesses the
data using that block's preprocessing pipeline, then projects it onto
the model space. If no block is specified, it delegates to the default
project method.

## Documentation strategy

\* This block documents \*\*all\*\* BaDA S3 methods via the shared
\`@rdname bada\` tag. Method‐specific files should therefore
\*\*omit\*\* the \`@return\` field and simply add \`@inheritParams
bada\` (or \`@inherit bada\`) plus any method‑specific \`@details\` or
\`@note\` sections.

## References

Abdi, H., Williams, L. J., & Bera, M. (2017). \*Barycentric discriminant
analysis\*. In \*\*Encyclopedia of Social Network Analysis and
Mining\*\* (pp. 1–20).
