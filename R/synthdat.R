#' Generate synthetic multiblock data
#'
#' This helper function creates several blocks of multivariate data that share
#' a common set of latent factor scores.  Optionally the variables can be placed
#' on the unit sphere to yield spatial coordinates and a sparse k-nearest-neighbour
#' graph.  The returned object also contains the ground-truth loadings and scores
#' used for simulation.
#'
#' @param S Number of subjects/blocks to generate (default 5).
#' @param n Number of observations (rows) per block (default 100).
#' @param p Number of variables (columns) per block, or a vector of length S
#'   specifying different dimensions per block (default 200).
#' @param r The rank of the shared component structure (default 3).
#' @param sigma The standard deviation of the noise added to the data (default 0.1).
#' @param sphere Logical; if TRUE, variables are placed on the unit sphere and
#'   a k-nearest-neighbour graph is computed (default FALSE).
#' @param k_nn Number of nearest neighbors for spatial correlations when
#'   `sphere = TRUE` (default 6).
#' @param seed Random seed for reproducibility (default 1).
#'
#' @return A list containing:
#'   \describe{
#'     \item{data_list}{A list of data matrices, one per block.}
#'     \item{coords_list}{Coordinates on the unit sphere (if `sphere = TRUE`), otherwise NULL.}
#'     \item{V_true}{List of ground-truth loading matrices.}
#'     \item{F_true}{Ground-truth factor score matrix.}
#'     \item{Sadj}{Spatial adjacency Laplacian matrix (if `sphere = TRUE`), otherwise NULL.}
#'   }
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
    # column‑centre without adding scale() attributes
    Xs <- sweep(Xs, 2, colMeans(Xs), "-")
    data_lst[[s]] <- Xs
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
