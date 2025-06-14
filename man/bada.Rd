% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/bada.R
\name{bada}
\alias{bada}
\alias{bootstrap.bada}
\alias{bada.multidesign}
\alias{reprocess.bada}
\alias{project.bada}
\title{Barycentric Discriminant Analysis (generic)}
\usage{
bada(data, y, subject, preproc, ncomp, ...)

\method{bootstrap}{bada}(x, data, nboot = 500, alpha = 0.05, verbose = FALSE, ...)

\method{bada}{multidesign}(
  data,
  y,
  subject,
  preproc = multivarious::center(),
  ncomp = 2,
  resdim = 20,
  rescomp = 0,
  ...
)

\method{reprocess}{bada}(x, new_data, colind = NULL, block = NULL)

\method{project}{bada}(x, new_data, block)
}
\arguments{
\item{data}{The dataset on which the bootstrap resampling should be performed.}

\item{y}{A factor or bare column name giving the grouping variable.}

\item{subject}{A factor or bare column name identifying repeated‑measures, i.e.
the subject/block variable.}

\item{preproc}{A \link[multivarious]{preprocessing pipeline}.  Defaults to
`multivarious::center()` in the concrete methods.}

\item{ncomp}{Integer; the number of discriminant components to compute.}

\item{...}{Additional arguments passed to methods (currently unused).}

\item{x}{A fitted BaDA model object.}

\item{nboot}{An integer specifying the number of bootstrap resamples to perform (default is 500).}

\item{alpha}{The percentile level for the computation of the lower and upper percentiles (default is 0.05).}

\item{verbose}{Logical; if `TRUE`, progress information is printed during
resampling. Defaults to `FALSE`.}

\item{resdim}{pca dimensionality for residual analysis (only relevant if `rescomp` > 0)}

\item{rescomp}{number of final residual components (default = 0, no residual aanalysis)}

\item{new_data}{A numeric matrix of new data to be projected.}

\item{colind}{An optional integer vector specifying column indices to use within blocks.
If NULL and block is also NULL, all blocks are used. If NULL but block is provided,
all columns in the specified block are used.}

\item{block}{An optional character string specifying which block's preprocessing
to apply before projection. If missing, the generic method is used.}
}
\value{
An object that inherits from class `bada`.  See the help for the
  corresponding method (e.g. [\code{bada.multidesign}]) for the exact
  structure.

A preprocessed numeric matrix with the same number of rows as the input data.

A numeric matrix of projected scores.
}
\description{
A generic function for fitting a **Barycentric Discriminant Analysis** (BaDA) model
to multi–subject multivariate data.  Concrete methods (e.g. for objects of class
`multidesign`) implement the actual estimation procedure and return an object that
inherits from class `\code{bada}`.

Perform bootstrap resampling on a multivariate bada model to estimate the variability 
of components and scores.

This function transforms new data using the preprocessing pipeline from a fitted BaDA model.
It handles different preprocessing scenarios based on the provided parameters.

This function projects new data onto a previously fitted BaDA model,
returning the scores of the new data in the space of the original model.
}
\details{
The function returns a list containing the summarized bootstrap resampled components and scores for the model. 
The returned list contains eight elements:
   * `boot_scores_mean`: A matrix which is the mean of all bootstrapped scores matrices.
   * `boot_scores_sd`: A matrix which is the standard deviation of all bootstrapped scores matrices.
   * `boot_scores_upper`: A matrix which is the upper alpha percentile of all bootstrapped scores matrices.
   * `boot_scores_lower`: A matrix which is the lower alpha percentile of all bootstrapped scores matrices.
   * `boot_lds_mean`: A matrix which is the mean of all bootstrapped loadings matrices.
   * `boot_lds_sd`: A matrix which is the standard deviation of all bootstrapped loadings matrices.
   * `boot_lds_upper`: A matrix which is the upper alpha percentile of all bootstrapped loadings matrices.
   * `boot_lds_lower`: A matrix which is the lower alpha percentile of all bootstrapped loadings matrices.
The dimensions of each matrix in the list correspond to the dimensions of the respective matrices in the input data.

The function handles three scenarios:
1. When both colind and block are NULL: The function applies preprocessing for each block
   and averages the results.
2. When block is provided: The function applies preprocessing specific to the named block.
   If colind is also provided, it acts as a relative subset within the block.
3. When only colind is provided: The function applies preprocessing for each block
   using the specified column indices and averages the results.

When a specific block is provided, the function first reprocesses the data using
that block's preprocessing pipeline, then projects it onto the model space.
If no block is specified, it delegates to the default project method.
}
\section{Documentation strategy}{

  * This block documents **all** BaDA S3 methods via the shared
    `@rdname bada` tag.  Method‐specific files should therefore **omit** the
    `@return` field and simply add `@inheritParams bada` (or `@inherit bada`)
    plus any method‑specific `@details` or `@note` sections.
}

\references{
Abdi, H., Williams, L. J., & Bera, M. (2017). *Barycentric discriminant
analysis*. In **Encyclopedia of Social Network Analysis and Mining** (pp.
1–20).
}
