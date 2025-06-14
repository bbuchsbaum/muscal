% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synthdat.R
\name{synthetic_multiblock}
\alias{synthetic_multiblock}
\title{Generate synthetic multiblock data}
\usage{
synthetic_multiblock(
  S = 5,
  n = 100,
  p = 200,
  r = 3,
  sigma = 0.1,
  sphere = FALSE,
  k_nn = 6,
  seed = 1
)
}
\arguments{
\item{S}{Number of blocks/subjects.}

\item{n}{Number of rows (observations) per block.}

\item{p}{Number of columns (variables) per block. Can be a single value or a
vector of length `S`.}

\item{r}{Rank (number of latent components).}

\item{sigma}{Standard deviation of the Gaussian noise.}

\item{sphere}{Logical; if `TRUE`, generate coordinates on the unit sphere and
build a k-nearest-neighbour graph.}

\item{k_nn}{Number of neighbours for the graph when `sphere = TRUE`.}

\item{seed}{Integer seed for reproducible generation.}
}
\value{
A list with the following elements:
\itemize{
  \item `data_list` — centred data matrices `X_s` of size `n` \times `p_s`;
  \item `coords_list` — 3-D coordinates for each block when `sphere = TRUE`;
  \item `V_true` — list of true loading matrices (`p_s` \times `r`);
  \item `F_true` — matrix of shared factor scores (`n` \times `r`);
  \item `Sadj` — sparse Laplacian matrix for smoothness penalties.
}
}
\description{
This helper function creates several blocks of multivariate data that share
a common set of latent factor scores.  Optionally the variables can be placed
on the unit sphere to yield spatial coordinates and a sparse k-nearest-neighbour
graph.  The returned object also contains the ground-truth loadings and scores
used for simulation.
}
