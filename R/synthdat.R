#' Generate synthetic multiblock data
#'
#' This helper function creates several blocks of multivariate data that share
#' a common set of latent factor scores.  Optionally the variables can be placed
#' on the unit sphere to yield spatial coordinates and a sparse k-nearest-neighbour
#' graph.  The returned object also contains the ground-truth loadings and scores
#' used for simulation.
#'
#' @param S Number of subjects/blocks to generate.
#' @param n Number of observations (rows) per block.
#' @param p_list A list of integers specifying the number of variables (columns) for each block.
#' @param ncomp The rank of the shared component structure.
#' @param g_list A list specifying the shared structure for each block, where each element is a matrix of size \eqn{p_i \times ncomp}.
#' @param u A matrix of size \eqn{n \times ncomp} representing the subject-specific scores on the shared components.
#' @param sigma The standard deviation of the noise added to the data.
#' @param orthogonal_space Logical; if TRUE, the block-specific scores are orthogonal to the shared scores.
#' @param A A matrix of size \eqn{p_i \times p_i} used to induce spatial correlation (if `spatial=TRUE`).
#' @param ... extra arguments.
#' @param spatial whether to add spatial correlations.
#' @param k_nn number of nearest neighbors for spatial correlations.
#'
#' @return A list of matrices, where each matrix is a data block.
#' @export
synthetic_multiblock <- function(S       = 5,
                                 n       = 100,
                                 p       = 200,
                                 r       = 3,
                                 sigma   = 0.1,
                                 sphere  = FALSE,
                                 k_nn    = 6,
                                 seed    = 1)
{
  if (length(p) == 1) p <- rep(p, S)
  stopifnot(length(p) == S, r < min(p))

  set.seed(seed)

  # ------------------------------------------------ shared factors
  F_true <- matrix(rnorm(n * r), n, r)                    # N(0,1)
  F_true <- scale(F_true, scale = FALSE)                  # centre
  F_true <- qr.Q(qr(F_true))                              # orthogonal

  # ------------------------------------------------ per‑block loadings
  V_true   <- vector("list", S)
  data_lst <- vector("list", S)

  for (s in seq_len(S)) {
    # start from a common orthogonal basis ...
    B  <- svd(matrix(rnorm(p[s] * r), p[s], r))$u
    # ... then apply a random orthogonal rotation (makes blocks non‑identical)
    Q  <- svd(matrix(rnorm(r * r), r, r))$u
    V_true[[s]] <- B %*% Q                               # p_s × r   orthonormal
    # synthesise data
    Xs <- F_true %*% t(V_true[[s]]) + sigma * matrix(rnorm(n * p[s]), n, p[s])
    # column‑centre (as assumed by MFA code)
    data_lst[[s]] <- scale(Xs, scale = FALSE)
  }

  # ------------------------------------------------ optional spatial coords
  coords_lst <- NULL
  Sadj       <- NULL
  if (sphere) {
    if (!requireNamespace("RANN", quietly = TRUE) ||
        !requireNamespace("Matrix", quietly = TRUE))
      stop("Install packages 'RANN' and 'Matrix' for the spatial variant.")

    # stack block‑wise coords on the unit sphere
    coords_lst <- lapply(p, function(ps) {
      # generate uniformly on sphere via normalisation
      mat <- matrix(rnorm(ps * 3), ps, 3)
      mat / sqrt(rowSums(mat^2))
    })
    coords_all <- do.call(rbind, coords_lst)

    # k‑NN adjacency (sparse)
    nn   <- RANN::nn2(coords_all, k = k_nn + 1)$nn.idx[, -1]  # drop self
    i    <- rep(seq_len(nrow(coords_all)), each = k_nn)
    j    <- as.vector(t(nn))
    w    <- rep(1, length(i))
    A    <- Matrix::sparseMatrix(i = i, j = j, x = w,
                                 dims = c(nrow(coords_all), nrow(coords_all)))
    A    <- Matrix::forceSymmetric(A, uplo = "U")             # make undirected
    deg  <- Matrix::rowSums(A)
    Sadj <- Matrix::Diagonal(x = deg) - A                     # Laplacian
  }

  list(data_list   = data_lst,
       coords_list = coords_lst,
       V_true      = V_true,
       F_true      = F_true,
       Sadj        = Sadj)
}