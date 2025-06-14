% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/penalized_mfa.R
\name{penalized_mfa}
\alias{penalized_mfa}
\alias{penalized_mfa.list}
\alias{penalized_mfa.multiblock}
\alias{penalized_mfa.multidesign}
\title{Penalized Multiple Factor Analysis (MFA)}
\usage{
penalized_mfa(data, ...)

\method{penalized_mfa}{list}(
  data,
  ncomp = 2,
  lambda = 1,
  penalty_method = c("projection", "pairwise", "global_mean"),
  max_iter = 10,
  nsteps_inner = 5,
  learning_rate = 0.01,
  optimizer = c("adam", "gradient"),
  preproc = multivarious::center(),
  beta1 = 0.9,
  beta2 = 0.999,
  adam_epsilon = 1e-08,
  tol_obj = 1e-07,
  tol_inner = NULL,
  compute_consensus = FALSE,
  verbose = FALSE,
  ...
)

\method{penalized_mfa}{multiblock}(data, ...)

\method{penalized_mfa}{multidesign}(data, subject, ...)
}
\arguments{
\item{data}{A list of matrices, a `multiblock` object, or a `multidesign` object.}

\item{...}{Additional arguments passed to methods.}

\item{ncomp}{Integer number of latent components.}

\item{lambda}{Non-negative scalar controlling the penalty strength.}

\item{penalty_method}{One of `"pairwise"`, `"global_mean"`, or `"projection"`.}

\item{max_iter}{Maximum number of outer BCD iterations.}

\item{nsteps_inner}{Number of gradient updates per block in each iteration.}

\item{learning_rate}{Step size for the optimizer.}

\item{optimizer}{One of `"gradient"` (fixed step size) or `"adam"` (adaptive).}

\item{preproc}{A `pre_processor` object from `multivarious` (e.g., `center()`), a list of such objects, or `NULL` for no preprocessing.}

\item{beta1, beta2}{Adam hyperparameters.}

\item{adam_epsilon}{Small constant for Adam's denominator.}

\item{tol_obj}{Numeric tolerance for the relative change in the objective function to determine convergence.}

\item{tol_inner}{Optional tolerance for stopping inner loop updates based on the norm of the change in `V`.}

\item{compute_consensus}{If `TRUE`, computes a consensus (average) loading matrix. Uses Procrustes alignment via the `vegan` package if installed.}

\item{verbose}{If `TRUE`, prints iteration logs using the `cli` package.}

\item{subject}{Required for `multidesign` method: the name of the subject variable.}
}
\value{
A `multiblock_projector` object of class `penalized_mfa`, containing the concatenated loadings (`v`), block indices, and preprocessor.
  Additional results are stored as attributes:
  \item{V_list}{The list of final, orthonormal loading matrices for each block.}
  \item{obj_values}{A vector of objective function values at each iteration.}
  \item{consensus}{The consensus loading matrix, if `compute_consensus=TRUE`.}
  \item{lambda, penalty_method, iterations_run}{Model and convergence information.}
}
\description{
This function implements a penalized MFA-like decomposition via block-coordinate
descent (BCD). Each block \(\mathbf{X}_i \in \mathbb{R}^{n \times p_i}\) has a loading
matrix \(\mathbf{V}_i \in \mathbb{R}^{p_i \times k}\) with orthonormal columns.
It assumes data blocks `X_i` are pre-processed (e.g., column-centered) as needed.

We solve:
\deqn{
  \min_{\{V_i\}} \sum_{i=1}^S \|X_i - X_i V_i V_i^\top\|_F^2
  \;+\; \lambda \,\mathrm{Penalty}(\{V_i\}),
}
where the penalty is one of:
\itemize{
  \item \code{"pairwise"}: \(\sum_{i < j} \|V_i - V_j\|_F^2\). Penalizes Euclidean distance between loading matrices. Not rotation-invariant.
  \item \code{"global_mean"}: \(\sum_{i} \|V_i - \bar{V}\|_F^2\). Similar to "pairwise" but computationally simpler. Not rotation-invariant.
  \item \code{"projection"}: \(\sum_{i} \|P_i - \bar{P}\|_F^2\), where \(P_i = V_i V_i^\top\) is the projection matrix for block `i`. This penalty is rotation-invariant.
}
}
\section{Input types}{

The function supports several input types for `data`: `list`, `multiblock`, and `multidesign`.
}

\examples{
\dontrun{
# Example with a list of matrices
data_list <- lapply(1:3, function(i) scale(matrix(rnorm(100), 10, 10), scale=FALSE))
res <- penalized_mfa(data_list, ncomp=2, lambda=1, penalty_method="projection",
                     optimizer="adam", max_iter=50, verbose=TRUE)

# Example with a multiblock object
mb <- multiblock(data_list)
res_mb <- penalized_mfa(mb, ncomp=2, lambda=1)
}

}
