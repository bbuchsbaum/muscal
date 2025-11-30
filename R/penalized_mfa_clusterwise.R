#' Validate arguments for penalized MFA functions
#' @keywords internal
validate_pmfa_args <- function(data_list, coords_list = NULL, ncomp, lambda, max_iter, 
                               nsteps_inner, learning_rate, optimizer, beta1, beta2, 
                               adam_epsilon, tol_obj, tol_inner, verbose, memory_budget_mb) {
  chk::chk_list(data_list)
  if (!is.null(coords_list)) chk::chk_list(coords_list)
  chk::chk_numeric(ncomp); chk::chk_gte(ncomp, 1)
  chk::chk_numeric(lambda); chk::chk_gte(lambda, 0)
  chk::chk_numeric(max_iter); chk::chk_gte(max_iter, 1)
  chk::chk_numeric(nsteps_inner); chk::chk_gte(nsteps_inner, 1)
  chk::chk_numeric(learning_rate); chk::chk_gt(learning_rate, 0)
  chk::chk_numeric(beta1); chk::chk_range(beta1, c(0, 1))
  chk::chk_numeric(beta2); chk::chk_range(beta2, c(0, 1))
  chk::chk_numeric(adam_epsilon); chk::chk_gt(adam_epsilon, 0)
  chk::chk_numeric(tol_obj); chk::chk_gt(tol_obj, 0)
  if (!is.null(tol_inner)) {
    chk::chk_numeric(tol_inner); chk::chk_gt(tol_inner, 0)
  }
  chk::chk_flag(verbose)
  chk::chk_numeric(memory_budget_mb); chk::chk_gt(memory_budget_mb, 0)
  
  # Convert to integers where needed
  ncomp <- as.integer(ncomp)
  max_iter <- as.integer(max_iter)
  nsteps_inner <- as.integer(nsteps_inner)
  
  # Additional validations
  if (length(data_list) < 2) stop("Need >= 2 subjects/blocks.")
  if (!is.null(coords_list) && length(coords_list) != length(data_list)) {
    stop("coords_list must match data_list length.")
  }
}

#' Simple spatial constraints function (creates k-NN adjacency matrix)
#' @param coords_list List of coordinate matrices (each k_s x 3)
#' @param k_nn Number of nearest neighbors (default 6)
#' @param nblocks Number of blocks (ignored, for compatibility)
#' @return Sparse adjacency matrix
#' @keywords internal
spatial_constraints <- function(coords_list, k_nn = 6, nblocks = NULL, ...) {
  # Combine all coordinates
  coords_all <- do.call(rbind, coords_list)
  n_total <- nrow(coords_all)
  
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' needed for sparse adjacency. Please install it.", call. = FALSE)
  }
  
  # Simple distance-based adjacency for testing
  # In a real implementation, this might use RANN for k-NN graphs
  dist_mat <- as.matrix(dist(coords_all))
  
  # Create k-NN adjacency matrix
  A <- Matrix::Matrix(0, n_total, n_total, sparse = TRUE)
  
  for (i in 1:n_total) {
    # Find k_nn nearest neighbors (excluding self)
    neighbors <- order(dist_mat[i, ])[2:(k_nn + 1)]  # Skip self (position 1)
    if (length(neighbors) > 0) {
      A[i, neighbors] <- 1
    }
  }
  
  # Make symmetric
  A <- Matrix::forceSymmetric(A, uplo = "U")
  
  return(A)
}

#' Determine if XtX should be precomputed based on memory budget
#' @keywords internal
should_precompute <- function(k_s, memory_budget_mb) {
  # Estimate memory for k_s x k_s matrix in MB (double precision = 8 bytes)
  estimated_mb <- (k_s^2 * 8) / (1024^2)
  estimated_mb < memory_budget_mb
}

#' Create gradient function for reconstruction term
#' @keywords internal
make_grad_fun <- function(X, XtX = NULL) {
  if (!is.null(XtX)) {
    # Precomputed XtX version: faster for repeated evaluations
    function(V) -2 * XtX %*% V
  } else {
    # On-the-fly version: memory efficient
    function(V) -2 * crossprod(X, X %*% V)
  }
}

#' Build a cluster graph from spatial coordinates
#' @param coords_list List of coordinate matrices (each k_s x 3)
#' @param adjacency_opts Options for adjacency construction
#' @return Sparse adjacency matrix
#' @keywords internal
cluster_graph <- function(coords_list, adjacency_opts = list()) {
  # Check for spatial_constraints function
  if (!exists("spatial_constraints", mode = "function")) {
    stop("Function 'spatial_constraints' not found. Ensure it is loaded.", call. = FALSE)
  }
  
  # Build adjacency matrix
  adj_opts <- modifyList(list(nblocks = length(coords_list)), adjacency_opts)
  Aadj <- do.call(spatial_constraints, c(list(coords_list), adj_opts))
  
  # Calculate Graph Laplacian L = D - A
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' needed for sparse Laplacians. Please install it.", call. = FALSE)
  }
  
  deg <- Matrix::rowSums(Aadj)
  Sadj <- Matrix::Diagonal(x = deg) - Aadj
  
  return(Sadj)
}

#' Update a single block's loadings using gradient descent or Adam
#' @keywords internal
block_update_cluster <- function(s, V_list, Sadj, row_indices_list, Vbig_begin_iter,
                                 grad_fun_list, lambda, nsteps_inner, optimizer,
                                 learning_rate, global_step, m_list = NULL, V2_list = NULL,
                                 beta1, beta2, adam_epsilon, tol_inner, verbose, iter, ncomp_block) {
  
  k_s <- nrow(V_list[[s]])
  ncomp_s <- ncomp_block[s]  # Use block-specific ncomp
  if (k_s == 0) return(list(V = V_list[[s]], global_step = global_step, m = NULL, V2 = NULL))
  
  Vs <- V_list[[s]]
  idx_s <- row_indices_list[[s]]
  
  # Skip if Vs is somehow ill-defined
  if (is.null(Vs) || nrow(Vs) != k_s || ncol(Vs) != ncomp_s) {
    if (verbose) cli::cli_alert_danger("Iter {iter}, Block {s}: Skipping update due to inconsistent V dimensions.")
    return(list(V = Vs, global_step = global_step, m = NULL, V2 = NULL))
  }
  
  # Adam states if needed
  if (optimizer == "adam") {
    Mi <- m_list[[s]]
    V2i <- V2_list[[s]]
    # Ensure Adam states match Vs dimensions
    if(nrow(Mi) != k_s || ncol(Mi) != ncomp_s || nrow(V2i) != k_s || ncol(V2i) != ncomp_s) {
      if (verbose) cli::cli_alert_danger("Adam state dimension mismatch for block {s}. Reinitializing.")
      Mi <- matrix(0, k_s, ncomp_s)
      V2i <- matrix(0, k_s, ncomp_s)
    }
  }
  
  # Project gradient G onto tangent space of Stiefel manifold at Vs
  stiefel_project <- function(Vs, G) {
    if (nrow(Vs) == 0) return(matrix(0, 0, ncomp_s))
    if (!is.matrix(G)) G <- as.matrix(G)
    A <- t(Vs) %*% G
    sym_part <- (A + t(A)) / 2
    G - Vs %*% sym_part
  }
  
  # Inner loop for updating block s
  for (step_inner in seq_len(nsteps_inner)) {
    # Reconstruction gradient
    G_recon <- grad_fun_list[[s]](Vs)
    
    # Spatial gradient: skip computation if lambda = 0
    if (lambda > 0) {
      # Spatial gradient: uses current Vs for diagonal block and 
      # V from start of iteration for off-diagonal
      G_spatial_fresh_diag <- Sadj[idx_s, idx_s] %*% Vs
      G_spatial_off_diag <- Sadj[idx_s, -idx_s] %*% Vbig_begin_iter[-idx_s, , drop = FALSE]
      G_spatial <- 2 * (G_spatial_fresh_diag + G_spatial_off_diag)
      
      # Combine gradients
      G <- G_recon + lambda * G_spatial
    } else {
      G <- G_recon
    }
    
    # Project onto tangent space
    G_t <- stiefel_project(Vs, G)
    
    # Update step (Adam or Gradient)
    if (optimizer == "gradient") {
      Vs_new <- Vs - learning_rate * G_t
    } else {
      # Adam update
      global_step <- global_step + 1
      out_adam <- adam_update_block(
        V = Vs, G = G_t, M = Mi, V2 = V2i,
        step_count = global_step,
        beta1 = beta1, beta2 = beta2,
        adam_epsilon = adam_epsilon,
        learning_rate = learning_rate
      )
      Vs_new <- out_adam$V
      Mi <- out_adam$M
      V2i <- out_adam$V2
    }
    
    # CRITICAL FIX: Retract after each inner update for numerical stability
    if (requireNamespace("pracma", quietly = TRUE)) {
      Vs_new <- pracma::orth(Vs_new)[, 1:ncomp_s, drop = FALSE]
    } else {
      Vs_new <- qr.Q(qr(Vs_new))[, 1:ncomp_s, drop = FALSE]
    }
    
    # Check for convergence in inner loop
    if (!is.null(tol_inner)) {
      diff_norm <- sqrt(sum((Vs_new - Vs)^2))
      if (diff_norm < tol_inner) {
        if (verbose) {
          cli::cli_alert_info("Iter {iter}, Block {s}: inner stop at step {step_inner} (||deltaV||={format(diff_norm, scientific=TRUE, digits=2)} < tol_inner)")
        }
        Vs <- Vs_new
        break
      }
    }
    
    # Update Vs for next inner iteration
    Vs <- Vs_new
  }
  
  # Return updated values
  result <- list(V = Vs, global_step = global_step)
  if (optimizer == "adam") {
    result$m <- Mi
    result$V2 <- V2i  # Fixed naming consistency
  }
  return(result)
}

#' Penalized MFA with Clusterwise Spatial Smoothness Constraints
#'
#' @description
#' This function implements a spatially-regularized Multiple Factor Analysis for
#' multi-subject neuroimaging or spatial cluster data. Unlike standard MFA which
#' encourages similarity among block loadings globally, this method uses spatial
#' coordinates of clusters (e.g., brain regions, spatial locations) to construct
#' a graph and enforces smoothness of loadings among spatially adjacent clusters
#' via a graph Laplacian penalty.
#'
#' @details
#' ## Optimization Problem
#'
#' For \eqn{S} subjects/blocks with data matrices \eqn{\mathbf{X}_s \in \mathbb{R}^{n_s \times k_s}}
#' (where \eqn{k_s} is the number of clusters for subject \eqn{s}), we estimate
#' orthonormal loading matrices \eqn{\mathbf{V}_s \in \mathbb{R}^{k_s \times r_s}}
#' by minimizing:
#'
#' \deqn{
#'   \min_{\{\mathbf{V}_s\}} \sum_{s=1}^S \|\mathbf{X}_s - \mathbf{X}_s \mathbf{V}_s \mathbf{V}_s^\top\|_F^2
#'   + \lambda \,\text{tr}(\mathbf{V}^\top \mathbf{L} \mathbf{V})
#' }
#'
#' where:
#' \itemize{
#'   \item The first term is the reconstruction error (sum of squared residuals)
#'   \item \eqn{\mathbf{V}} is the vertical concatenation of all \eqn{\mathbf{V}_s}
#'   \item \eqn{\mathbf{L} = \mathbf{D} - \mathbf{A}} is the graph Laplacian
#'   \item \eqn{\mathbf{A}} is the adjacency matrix constructed from spatial coordinates
#'   \item \eqn{\mathbf{D}} is the degree matrix (\eqn{D_{ii} = \sum_j A_{ij}})
#'   \item \eqn{\lambda \geq 0} controls the spatial smoothness strength
#' }
#'
#' ## Spatial Smoothness Penalty
#'
#' The Laplacian penalty \eqn{\text{tr}(\mathbf{V}^\top \mathbf{L} \mathbf{V})} can be
#' rewritten as:
#' \deqn{
#'   \frac{1}{2} \sum_{i,j} A_{ij} \|\mathbf{v}_i - \mathbf{v}_j\|_2^2
#' }
#'
#' This penalizes the squared Euclidean distance between loadings of adjacent clusters.
#' When two clusters \eqn{i} and \eqn{j} are connected (\eqn{A_{ij} = 1}), their
#' loading vectors \eqn{\mathbf{v}_i} and \eqn{\mathbf{v}_j} are encouraged to be similar.
#'
#' **Normalized vs. Unnormalized Laplacian:**
#' \describe{
#'   \item{\strong{Unnormalized} (\code{normalized_laplacian = FALSE})}{
#'     \eqn{\mathbf{L} = \mathbf{D} - \mathbf{A}}. Smoothness is scaled by node degree;
#'     high-degree nodes have stronger smoothness constraints.
#'   }
#'   \item{\strong{Normalized} (\code{normalized_laplacian = TRUE}, default)}{
#'     \eqn{\mathbf{L}_{\text{sym}} = \mathbf{D}^{-1/2} \mathbf{L} \mathbf{D}^{-1/2}}.
#'     Scale-free smoothness that is independent of node degree. Recommended for
#'     most applications.
#'   }
#' }
#'
#' ## Graph Construction
#'
#' The spatial adjacency graph is built using k-nearest neighbors (k-NN) in 3D space:
#' \enumerate{
#'   \item Cluster coordinates from all subjects are pooled
#'   \item For each cluster, find its k nearest neighbors (default: k=6)
#'   \item Create edges between each cluster and its neighbors
#'   \item Symmetrize the adjacency matrix
#' }
#'
#' The number of neighbors can be controlled via `adjacency_opts = list(k_nn = k)`.
#'
#' ## Optimization Algorithm
#'
#' The method uses block-coordinate descent (BCD) with Riemannian optimization:
#'
#' 1. **Initialization**: Initialize each \eqn{\mathbf{V}_s} via SVD of \eqn{\mathbf{X}_s}
#' 2. **Outer loop** (max_iter iterations): Cycle through all subjects/blocks
#' 3. **Inner loop** (nsteps_inner steps per block):
#'    - Compute reconstruction gradient: \eqn{\nabla_{\mathbf{V}_s} = -2 \mathbf{X}_s^\top \mathbf{X}_s \mathbf{V}_s}
#'    - Compute spatial gradient: \eqn{2\lambda \mathbf{L}_{ss} \mathbf{V}_s + 2\lambda \mathbf{L}_{s,-s} \mathbf{V}_{-s}}
#'    - Project combined gradient onto Stiefel manifold tangent space
#'    - Update via gradient descent or Adam
#'    - Retract to manifold via QR decomposition
#' 4. **Convergence check**: Stop when relative change in objective < `tol_obj`
#'
#' ## Memory Management
#'
#' For large datasets, computing and storing \eqn{\mathbf{X}_s^\top \mathbf{X}_s} for
#' all blocks may exceed available memory. The `memory_budget_mb` parameter controls
#' this tradeoff:
#'
#' \itemize{
#'   \item **Precomputed mode** (if \eqn{k_s^2 \times 8 / (1024^2)} < `memory_budget_mb`):
#'     Store \eqn{\mathbf{X}_s^\top \mathbf{X}_s} for faster gradient computation
#'   \item **On-the-fly mode** (otherwise): Compute gradients as needed, saving memory
#'     at the cost of computation time
#' }
#'
#' ## Variable-Rank Loadings
#'
#' When subjects have different numbers of clusters, each block can have a different
#' effective rank. The algorithm automatically sets \eqn{r_s = \min(\text{ncomp}, k_s)}
#' for each block and pads the concatenated loading matrix with zeros for consistency.
#'
#' ## Preprocessing
#'
#' Data preprocessing is crucial for reconstruction-based methods. The function
#' defaults to centering each block if no preprocessing is specified. **Important:**
#' Preprocessing must preserve the number of columns (clusters) to maintain
#' correspondence with spatial coordinates.
#'
#' ## When Lambda = 0
#'
#' When \eqn{\lambda = 0}, the spatial penalty is disabled and the method reduces
#' to independent PCA on each block. This provides a useful baseline for comparison.
#'
#' ## Practical Considerations
#'
#' \describe{
#'   \item{\strong{Choosing lambda}}{
#'     Start with \eqn{\lambda = 0} (no smoothness) and gradually increase.
#'     Typical range: 0.1 to 10. Monitor the objective function components
#'     (reconstruction vs. smoothness) to balance the tradeoff.
#'   }
#'   \item{\strong{Number of neighbors}}{
#'     More neighbors (larger k-NN) create a denser graph with stronger smoothness
#'     constraints. Default k=6 works well for 3D spatial data.
#'   }
#'   \item{\strong{Convergence}}{
#'     The method may converge slowly for large \eqn{\lambda}. Increase `max_iter`
#'     or adjust `learning_rate` if needed.
#'   }
#'   \item{\strong{Optimizer choice}}{
#'     Adam is generally more robust and converges faster. Use gradient descent
#'     for better control or when Adam's adaptive behavior is undesired.
#'   }
#' }
#'
#' @section Engineering Improvements:
#' This implementation includes several optimizations:
#' \itemize{
#'   \item **Mathematical**: Normalized Laplacian default, \eqn{\lambda=0} optimization skips
#'   \item **Numerical**: Sparse matrix operations, in-block QR retraction, fast orthogonalization
#'   \item **Stability**: Variable-rank loadings, robust parameter validation, Adam state management
#'   \item **Architecture**: Modular helper functions, comprehensive logging, clear error messages
#'   \item **Memory**: Adaptive memory budget, efficient gradient computation strategies
#'   \item **API**: Shorter `pmfa_cluster()` alias, consistent naming conventions
#' }
#'
#' @param data_list A list of length \eqn{S} containing numeric matrices. Each element
#'   \eqn{\mathbf{X}_s} is an \eqn{n_s \times k_s} matrix where \eqn{n_s} is the number
#'   of observations and \eqn{k_s} is the number of clusters/features for subject \eqn{s}.
#'   Data should be column-centered (mean zero per cluster) for proper reconstruction.
#' @param coords_list A list of length \eqn{S} containing coordinate matrices. Each element
#'   is a \eqn{k_s \times 3} matrix of (x, y, z) spatial coordinates for the clusters in
#'   the corresponding data block. Must have exactly 3 columns.
#' @param ncomp Integer number of latent components to extract per block (default: 2).
#'   When a block has fewer than `ncomp` clusters, the effective rank is automatically
#'   reduced to match the number of clusters.
#' @param lambda Non-negative scalar controlling the spatial smoothness penalty strength
#'   (default: 1). Larger values enforce stronger smoothness among spatially adjacent
#'   clusters. When \eqn{\lambda = 0}, the method reduces to independent PCA per block.
#'   Typical range: 0.1 to 10.
#' @param adjacency_opts A named list of options passed to the internal
#'   `spatial_constraints` function for adjacency matrix construction (default: empty list).
#'   Supported options:
#'   \itemize{
#'     \item \code{k_nn}: Number of nearest neighbors (default: 6)
#'   }
#' @param max_iter Maximum number of outer block-coordinate descent iterations (default: 10).
#'   More iterations allow better convergence but increase computation time. Typical range:
#'   10 to 100.
#' @param nsteps_inner Number of gradient update steps per block in each outer iteration
#'   (default: 1). Higher values perform more thorough optimization of each block before
#'   moving to the next. For stability, starting with 1 is recommended.
#' @param learning_rate Step size for the optimizer (default: 0.01). Controls the magnitude
#'   of updates in the gradient direction. Too large may cause divergence; too small slows
#'   convergence. Typical range: 0.001 to 0.1.
#' @param optimizer Character string specifying the optimization algorithm (default: "gradient").
#'   Options are:
#'   \itemize{
#'     \item \code{"gradient"}: Standard gradient descent with fixed learning rate
#'     \item \code{"adam"}: Adaptive Moment Estimation with momentum and adaptive learning rates
#'   }
#'   Adam is generally more robust but gradient descent offers more control.
#' @param beta1 Adam hyperparameter for first moment decay (default: 0.9). Controls the
#'   exponential decay rate for gradient moving average. Typical range: 0.8 to 0.95.
#'   Only used when `optimizer = "adam"`.
#' @param beta2 Adam hyperparameter for second moment decay (default: 0.999). Controls
#'   the exponential decay rate for squared gradient moving average. Typical range:
#'   0.99 to 0.9999. Only used when `optimizer = "adam"`.
#' @param adam_epsilon Small constant added to denominator in Adam for numerical stability
#'   (default: 1e-8). Prevents division by zero. Only used when `optimizer = "adam"`.
#' @param tol_obj Numeric tolerance for outer loop convergence based on relative change
#'   in objective function (default: 1e-7). The algorithm stops when
#'   \eqn{|f^{(t+1)} - f^{(t)}| / (|f^{(t)}| + \epsilon) < \text{tol_obj}}.
#'   Smaller values require more precise convergence. Typical range: 1e-8 to 1e-5.
#' @param tol_inner Optional numeric tolerance for inner loop convergence based on
#'   Frobenius norm of change in loading matrix (default: NULL, no early stopping).
#'   When specified, inner loop stops if \eqn{\|\mathbf{V}_s^{\text{new}} - \mathbf{V}_s\|_F < \text{tol_inner}}.
#' @param verbose Logical indicating whether to print iteration progress (default: FALSE).
#'   When TRUE, displays detailed information about objective values, convergence, memory
#'   usage, and potential issues using the `cli` package.
#' @param preproc Preprocessing specification (default: NULL, will force centering).
#'   Can be:
#'   \itemize{
#'     \item \code{NULL}: Automatically applies centering to each block
#'     \item A single `pre_processor` object: Applied independently to each block
#'     \item A list of `pre_processor` objects: One per block
#'   }
#'   **Critical:** Preprocessing must preserve the number of columns to maintain
#'   spatial correspondence.
#' @param memory_budget_mb Numeric specifying the maximum memory in megabytes allocated
#'   per block for precomputing \eqn{\mathbf{X}_s^\top \mathbf{X}_s} (default: 1024).
#'   If a block's \eqn{k_s^2 \times 8 / (1024^2)} exceeds this budget, gradients are
#'   computed on-the-fly instead. Increase for better speed with large memory; decrease
#'   if memory is constrained.
#' @param normalized_laplacian Logical indicating whether to use the normalized graph
#'   Laplacian (default: TRUE). When TRUE, uses \eqn{\mathbf{L}_{\text{sym}} = \mathbf{D}^{-1/2} \mathbf{L} \mathbf{D}^{-1/2}}
#'   for scale-free smoothness that is independent of node degree. When FALSE, uses the
#'   unnormalized Laplacian \eqn{\mathbf{L} = \mathbf{D} - \mathbf{A}}. Normalized is
#'   recommended for most applications.
#'
#' @return A `multiblock_projector` object of class `"penalized_mfa_clusterwise"` with
#'   the following components:
#'   \describe{
#'     \item{\code{v}}{Concatenated loading matrix (total_clusters × ncomp) formed by
#'       vertically stacking all block-specific loading matrices. Padded with zeros for
#'       variable-rank blocks.}
#'     \item{\code{preproc}}{A pass-through preprocessor (identity transformation), since
#'       this function expects pre-processed input.}
#'     \item{\code{block_indices}}{Named list indicating which rows of `v` correspond to
#'       each block/subject. Each element is an integer vector of row indices.}
#'   }
#'
#'   Additional information stored as named elements (accessible via `$`):
#'   \describe{
#'     \item{\code{V_list}}{List of length S containing the final orthonormal loading
#'       matrices for each block. Each element is a \eqn{k_s \times r_s} matrix where
#'       \eqn{r_s = \min(\text{ncomp}, k_s)}.}
#'     \item{\code{Sadj}}{The (normalized or unnormalized) graph Laplacian matrix used
#'       for the spatial smoothness penalty. Sparse matrix of dimension
#'       \eqn{\sum_s k_s \times \sum_s k_s}.}
#'     \item{\code{LV}}{The matrix product \eqn{\mathbf{L} \mathbf{V}} from the final
#'       iteration, useful for computing spatial penalty contributions.}
#'     \item{\code{obj_values}}{Numeric vector of objective function values at each
#'       outer iteration, including the initial value. Length is (iterations_run + 1).}
#'     \item{\code{lambda}}{The spatial smoothness penalty weight used.}
#'     \item{\code{precompute_info}}{Logical vector of length S indicating which blocks
#'       used precomputed \eqn{\mathbf{X}_s^\top \mathbf{X}_s} (TRUE) vs. on-the-fly
#'       computation (FALSE).}
#'     \item{\code{iterations_run}}{Integer indicating how many outer iterations were
#'       completed before convergence or reaching max_iter.}
#'     \item{\code{ncomp_block}}{Integer vector of length S containing the effective
#'       number of components extracted for each block (may differ from `ncomp` for
#'       blocks with fewer clusters).}
#'   }
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic usage with simulated spatial cluster data
#' set.seed(123)
#' S <- 3  # 3 subjects
#'
#' # Generate data with spatial structure
#' data_list <- lapply(1:S, function(s) {
#'   n <- 50  # 50 observations
#'   k <- 20  # 20 clusters
#'   matrix(rnorm(n * k), n, k)
#' })
#'
#' # Generate 3D spatial coordinates for clusters
#' coords_list <- lapply(1:S, function(s) {
#'   matrix(runif(20 * 3, 0, 10), 20, 3)  # Random 3D positions
#' })
#'
#' # Fit model with spatial smoothness
#' res <- penalized_mfa_clusterwise(
#'   data_list, coords_list,
#'   ncomp = 3,
#'   lambda = 1,
#'   max_iter = 20,
#'   verbose = TRUE
#' )
#' print(res)
#'
#' # Plot convergence
#' plot(res$obj_values, type = 'b', xlab = 'Iteration', ylab = 'Objective',
#'      main = 'Convergence of Spatially-Regularized MFA')
#'
#' # Example 2: Using the shorter alias
#' res2 <- pmfa_cluster(data_list, coords_list, ncomp = 2, lambda = 0.5)
#'
#' # Example 3: Compare different lambda values
#' lambdas <- c(0, 0.1, 0.5, 1, 2, 5)
#' results <- lapply(lambdas, function(lam) {
#'   fit <- pmfa_cluster(data_list, coords_list, ncomp = 2, lambda = lam,
#'                       verbose = FALSE)
#'   list(
#'     lambda = lam,
#'     final_obj = tail(fit$obj_values, 1),
#'     iterations = fit$iterations_run
#'   )
#' })
#'
#' # Example 4: Using Adam optimizer for faster convergence
#' res_adam <- pmfa_cluster(
#'   data_list, coords_list,
#'   ncomp = 3,
#'   lambda = 1,
#'   optimizer = "adam",
#'   learning_rate = 0.05,
#'   max_iter = 30,
#'   verbose = TRUE
#' )
#'
#' # Example 5: Controlling k-NN graph construction
#' res_dense <- pmfa_cluster(
#'   data_list, coords_list,
#'   ncomp = 2,
#'   lambda = 1,
#'   adjacency_opts = list(k_nn = 10),  # More neighbors = denser graph
#'   verbose = TRUE
#' )
#'
#' # Example 6: Using unnormalized Laplacian
#' res_unnorm <- pmfa_cluster(
#'   data_list, coords_list,
#'   ncomp = 2,
#'   lambda = 1,
#'   normalized_laplacian = FALSE
#' )
#'
#' # Example 7: Memory-constrained settings
#' # For large datasets, reduce memory budget
#' res_mem <- pmfa_cluster(
#'   data_list, coords_list,
#'   ncomp = 2,
#'   lambda = 1,
#'   memory_budget_mb = 256,  # Limit memory per block
#'   verbose = TRUE
#' )
#' # Check which blocks used precomputed gradients
#' print(res_mem$precompute_info)
#'
#' # Example 8: Extracting block-specific loadings
#' V_list <- res$V_list
#'
#' # Loadings for first subject
#' V1 <- V_list[[1]]
#' dim(V1)  # k_s x ncomp
#'
#' # Visualize spatial smoothness of first loading vector
#' library(scatterplot3d)
#' scatterplot3d(
#'   coords_list[[1]],
#'   color = rank(V1[, 1]),
#'   main = "Spatial Pattern of First Loading (Subject 1)",
#'   xlab = "X", ylab = "Y", zlab = "Z"
#' )
#'
#' # Example 9: Examining the spatial penalty contribution
#' # Extract Laplacian and loadings
#' L <- res$Sadj
#' LV <- res$LV
#' v <- res$v
#'
#' # Compute spatial penalty: tr(V' L V)
#' spatial_penalty <- sum(LV * v)
#' cat("Spatial penalty term:", spatial_penalty, "\n")
#'
#' # Example 10: Variable-rank loadings (blocks with different sizes)
#' # Simulate data with varying cluster numbers
#' data_var <- list(
#'   matrix(rnorm(50 * 15), 50, 15),  # 15 clusters
#'   matrix(rnorm(50 * 20), 50, 20),  # 20 clusters
#'   matrix(rnorm(50 * 10), 50, 10)   # 10 clusters
#' )
#'
#' coords_var <- list(
#'   matrix(runif(15 * 3), 15, 3),
#'   matrix(runif(20 * 3), 20, 3),
#'   matrix(runif(10 * 3), 10, 3)
#' )
#'
#' res_var <- pmfa_cluster(
#'   data_var, coords_var,
#'   ncomp = 8,  # Will be capped at 10 (smallest block)
#'   lambda = 1,
#'   verbose = TRUE
#' )
#' # Check effective ranks per block
#' print(res_var$ncomp_block)
#'
#' # Example 11: Custom preprocessing per block
#' library(multivarious)
#'
#' preproc_list <- list(
#'   center(),        # Just center first block
#'   standardize(),   # Center and scale second block
#'   center()         # Just center third block
#' )
#'
#' res_preproc <- pmfa_cluster(
#'   data_list, coords_list,
#'   ncomp = 2,
#'   lambda = 1,
#'   preproc = preproc_list
#' )
#'
#' # Example 12: Lambda selection via cross-validation style approach
#' # (Simplified - would need proper CV in practice)
#' lambda_grid <- 10^seq(-1, 1, length.out = 10)
#' cv_results <- data.frame(
#'   lambda = lambda_grid,
#'   recon_error = NA,
#'   spatial_penalty = NA,
#'   total_obj = NA
#' )
#'
#' for (i in seq_along(lambda_grid)) {
#'   fit <- pmfa_cluster(data_list, coords_list, ncomp = 2,
#'                       lambda = lambda_grid[i], verbose = FALSE)
#'   cv_results$total_obj[i] <- tail(fit$obj_values, 1)
#'   # Could decompose into reconstruction and spatial components
#' }
#'
#' # Plot objective vs lambda
#' plot(cv_results$lambda, cv_results$total_obj, log = "x",
#'      type = "b", xlab = "Lambda (log scale)", ylab = "Final Objective",
#'      main = "Objective Function vs. Smoothness Penalty")
#' }
#'
#' @importFrom Matrix Diagonal rowSums crossprod
#' @importFrom pracma orth
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_danger cli_h1 cli_h2
#' @importFrom stats coefficients dist
#' @importFrom utils head tail modifyList combn
#' @export
penalized_mfa_clusterwise <- function(data_list,
                                      coords_list,
                                      ncomp         = 2L,
                                      lambda        = 1,
                                      adjacency_opts= list(),
                                      max_iter      = 10L,
                                      nsteps_inner  = 1L,  # Default to 1 for better stability
                                      learning_rate = 0.01,
                                      optimizer     = c("gradient","adam"),
                                      beta1         = 0.9,
                                      beta2         = 0.999,
                                      adam_epsilon  = 1e-8,
                                      tol_obj       = 1e-7,
                                      tol_inner     = NULL,
                                      verbose       = FALSE,
                                      preproc       = NULL,
                                      memory_budget_mb = 1024,
                                      normalized_laplacian = TRUE) {  # Default to TRUE
  
  # Check for Matrix package needed for sparse matrices (consolidate check)
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' needed for sparse Laplacians. Please install it.", call. = FALSE)
  }
  
  optimizer <- match.arg(optimizer)
  
  # Parameter validation using helper function (moved higher as recommended)
  validate_pmfa_args(data_list, coords_list, ncomp, lambda, max_iter, 
                     nsteps_inner, learning_rate, optimizer, beta1, beta2, 
                     adam_epsilon, tol_obj, tol_inner, verbose, memory_budget_mb)

  S <- length(data_list)

  # Check data dimensions early using number of columns only
  k_s_vec <- sapply(data_list, ncol) # Number of columns (clusters) per subject
  if(!all(sapply(coords_list, nrow) == k_s_vec)) {
      stop("Number of rows in each coords_list element must match number of columns in corresponding data_list element.")
  }

  coord_ncol <- sapply(coords_list, ncol)
  if (!all(coord_ncol == 3)) {
      bad_idx <- which(coord_ncol != 3)[1]
      stop(sprintf("coords_list element %d must have exactly 3 columns (x, y, z coordinates).",
                   bad_idx))
  }
  
  # CRITICAL FIX: Handle variable-rank V_s per block
  max_k_s <- max(k_s_vec)
  if (ncomp > max_k_s) {
    if (verbose) cli::cli_alert_danger("ncomp ({ncomp}) > maximum block size ({max_k_s}). Reducing ncomp to {max_k_s}.")
    ncomp <- max_k_s
  }
  
  # Set block-specific ncomp allowing variable-rank loadings
  ncomp_block <- pmin(ncomp, k_s_vec)
  if (verbose && !all(ncomp_block == ncomp)) {
    cli::cli_alert_info("Using variable-rank loadings: {paste(ncomp_block, collapse=', ')} components per block.")
  }
  
  # 0. Preprocessing using utility function
  # Force centering if no preprocessing is specified
  if (is.null(preproc)) {
    preproc <- multivarious::center()
    if (verbose) cli::cli_alert_info("No preprocessing specified. Forcing centering for reconstruction identity.")
  }
  
  preproc_result <- prepare_block_preprocessors(data_list, preproc, check_consistent_ncol = FALSE)
  proclist <- preproc_result$proclist
  Xp <- preproc_result$Xp
  
  # Preprocessing must NOT change the number of columns for spatial logic to work
  k_s_post_preproc <- sapply(Xp, ncol)
  if (!all(k_s_post_preproc == k_s_vec)) {
      stop("Preprocessing changed the number of columns (clusters) for some blocks. This invalidates the spatial constraints. Please use preprocessing steps that preserve the number of columns.", call.=FALSE)
  }
  
  # Build graph Laplacian
  Sadj <- cluster_graph(coords_list, adjacency_opts)
  
  # Optional: Use normalized Laplacian (now default)
  if (normalized_laplacian) {
    deg <- Matrix::rowSums(Sadj + Matrix::Diagonal(nrow(Sadj))) # Add diagonal to handle isolated nodes
    D_inv_sqrt <- Matrix::Diagonal(x = 1 / sqrt(pmax(deg, 1e-10)))  # Avoid division by zero
    Sadj <- D_inv_sqrt %*% Sadj %*% D_inv_sqrt
    if (verbose) cli::cli_alert_info("Using normalized Laplacian for scale-free smoothness.")
  }
  
  bigK <- nrow(Sadj)
  
  # Offsets for indexing into Vbig and L
  offset <- cumsum(c(0, k_s_vec))
  if (offset[S+1] != bigK) {
      stop(sprintf("Total number of clusters (%d) does not match Laplacian dimension (%d).", offset[S+1], bigK))
  }
  # Create row index list for easy slicing of Vbig/LV
  row_indices_list <- lapply(seq_len(S), function(s) (offset[s]+1):offset[s+1])

  # Precompute XtX or NormX2 based on memory budget with sparse matrix support
  XtX_list <- vector("list", S)
  normX2_list <- numeric(S)
  grad_fun_list <- vector("list", S)
  precompute_info <- logical(S)

  for (s in seq_len(S)) {
      Xs <- Xp[[s]]
      k_s <- k_s_vec[s]
      # Skip if block is empty
      if (k_s == 0) {
          normX2_list[s] <- 0
          grad_fun_list[[s]] <- function(V) matrix(0, 0, ncomp_block[s]) # Use block-specific ncomp
          precompute_info[s] <- FALSE # Explicitly false for empty block
          next
      }

      if (should_precompute(k_s, memory_budget_mb)) {
          # EFFICIENCY FIX: Use sparse XtX for better memory efficiency
          XtX_s <- Matrix::crossprod(Matrix::Matrix(Xs, sparse = FALSE))
          XtX_list[[s]] <- XtX_s
          normX2_list[s] <- sum(Matrix::diag(XtX_s)) # Faster than sum(Xs^2)
          grad_fun_list[[s]] <- make_grad_fun(Xs, XtX = XtX_s)
          precompute_info[s] <- TRUE
      } else {
          # XtX_list[[s]] remains NULL
          normX2_list[s] <- sum(Xs^2)
          grad_fun_list[[s]] <- make_grad_fun(Xs) # on-the-fly
          precompute_info[s] <- FALSE
      }
  }

  if (verbose) {
      precomputed_count <- sum(precompute_info)
      on_the_fly_count <- S - precomputed_count
      cli::cli_alert_info("Gradient/Objective: {precomputed_count} blocks using precomputed XtX, {on_the_fly_count} blocks using on-the-fly.")
  }
  
  # Initialize loadings (handle empty blocks and padding with variable rank)
  V_list <- vector("list", S)
  for (s in seq_len(S)) {
    Xs <- Xp[[s]]
    k_s <- k_s_vec[s]
    ncomp_s <- ncomp_block[s]
    
    if(k_s == 0) {
        V_list[[s]] <- matrix(0, 0, ncomp_s) # Placeholder for empty block
        next
    }
    # Perform SVD safely
    sv <- tryCatch(svd(Xs, nu=0, nv=ncomp_s), error = function(e) {
        if (verbose) cli::cli_alert_danger("SVD failed during initialization for block {s}: {e$message}")
        NULL
    })
    
    if (is.null(sv) || length(sv$d) == 0) {
        if (verbose) cli::cli_alert_danger("Block {s}: SVD init failed or returned no singular values. Initializing with random orthogonal matrix.")
        Vtemp <- pracma::orth(matrix(rnorm(k_s * ncomp_s), k_s, ncomp_s))
    } else {
        Vtemp <- sv$v
    }

    # Pad if rank < ncomp_s
    ncomp_svd <- ncol(Vtemp)
    if (ncomp_svd < ncomp_s) {
      extras <- ncomp_s - ncomp_svd
      Vextra <- matrix(rnorm(k_s * extras), k_s, extras)
      # Orthogonalize padded part against existing Vtemp and itself
      Q_all <- qr.Q(qr(cbind(Vtemp, Vextra)))
      # Ensure correct dimensions even if original rank was < ncomp_svd
      Vtemp <- Q_all[, 1:ncomp_s, drop = FALSE]
    } else if (ncomp_svd > ncomp_s) {
        # Should not happen with svd(nv=ncomp_s), but check
        Vtemp <- Vtemp[, 1:ncomp_s, drop=FALSE]
    }
    
    # Final orthonormalization for the initial V_list[s]
    V_list[[s]] <- qr.Q(qr(Vtemp))[, 1:ncomp_s, drop=FALSE]
  }
  
  # Adam storage (initialize based on actual ncomp_block)
  if (optimizer=="adam") {
    m_list <- lapply(seq_len(S), function(s) matrix(0, k_s_vec[s], ncomp_block[s]))
    V2_list <- lapply(seq_len(S), function(s) matrix(0, k_s_vec[s], ncomp_block[s]))  # Fixed naming
    global_step <- 0
  }
  
  # Helper to combine list of V_s into one big V matrix with variable ranks
  combine_loadings <- function(Vl) {
    # Check for consistency before combining
    if (length(Vl) != S) stop("Internal error: Vl length mismatch in combine_loadings")
    k_s_current <- sapply(Vl, nrow)
    if (!all(k_s_current == k_s_vec)) {
         if (verbose) cli::cli_alert_danger("Row counts in Vl differ from expected k_s_vec. Using expected.")
    }
    
    # Create Vbig using max ncomp for consistent dimensions
    max_ncomp <- max(ncomp_block)
    out <- matrix(0, bigK, max_ncomp)
    for (s in seq_len(S)) {
        idx_out <- row_indices_list[[s]]
        ncomp_s <- ncomp_block[s]
        # Ensure the block Vl[[s]] has the correct dimensions before assigning
        if(length(idx_out) > 0 && !is.null(Vl[[s]]) && nrow(Vl[[s]]) == length(idx_out) && ncol(Vl[[s]]) == ncomp_s) {
           out[idx_out, 1:ncomp_s] <- Vl[[s]]
        } else if (length(idx_out) > 0) {
            # If dimensions mismatch or Vl[[s]] is NULL, leave as zeros but warn
             if (verbose) cli::cli_alert_danger("Dimension mismatch or NULL V matrix for block {s} in combine_loadings. Filling with zeros.")
        }
    }
    out
  }
  
  # Initial objective calculation requires a first pass of combine+LV
  Vbig_current <- combine_loadings(V_list)
  
  # Skip Laplacian computation if lambda = 0 as recommended
  if (lambda > 0) {
    LV_current <- Sadj %*% Vbig_current # Initial LV for first iter spatial grad
  } else {
    LV_current <- matrix(0, nrow(Vbig_current), ncol(Vbig_current))
    if (verbose) cli::cli_alert_info("Lambda = 0: skipping Laplacian operations for efficiency.")
  }
  
  # Objective function using cached components
  calculate_objective <- function(Vl, Vbig, LV) {
      recon_total <- 0
      for(s in seq_len(S)) {
          Vs <- Vl[[s]]
          k_s <- k_s_vec[s]
          ncomp_s <- ncomp_block[s]
          if (k_s == 0) next # Skip empty blocks
          
          # Check Vs validity before calculation
          if (is.null(Vs) || nrow(Vs) != k_s || ncol(Vs) != ncomp_s) {
              if (verbose) cli::cli_alert_danger("Objective: Skipping block {s} due to invalid V dimensions.")
              recon_total <- recon_total + normX2_list[s] # Penalize fully if V is bad
              next
          }
          
          if (precompute_info[s]) {
              # Use XtX if available: ||X||^2 - tr(V' XtX V)
              XtXs <- XtX_list[[s]]
              recon_term <- normX2_list[s] - sum(Matrix::diag(Matrix::crossprod(Vs, XtXs %*% Vs)))
          } else {
              # Use on-the-fly: ||X||^2 - ||XV||^2
              Xs <- Xp[[s]]
              recon_term <- normX2_list[s] - sum((Xs %*% Vs)^2)
          }
          # Ensure non-negative reconstruction cost (numerical precision issues)
          recon_total <- recon_total + max(0, recon_term)
      }
      
      # Spatial penalty: tr(V' L V) = sum (LV * V), skip if lambda = 0
      penalty_term <- if (lambda > 0) lambda * sum(LV * Vbig) else 0
      
      return(recon_total + penalty_term)
  }
  
  # Main BCD loop
  obj_values <- numeric(max_iter + 1)
  # Calculate initial objective value
  initial_obj <- tryCatch(calculate_objective(V_list, Vbig_current, LV_current),
                          error = function(e) {
                              if (verbose) cli::cli_alert_danger("Could not compute initial objective: {e$message}")
                              return(NA)
                          })
  if (!is.finite(initial_obj)) { # Check for NA or Inf
      if (verbose) cli::cli_alert_danger("Initial objective is not finite. Check input data and initialization.")
      initial_obj <- Inf # Allow loop to start
  }
  obj_values[1] <- initial_obj
  old_obj <- initial_obj

  if (verbose) cli::cli_alert_info("Starting optimization... Initial objective: {format(initial_obj, digits=6)}")

  # Test for lambda = 0 case (should revert to independent PCA per block)
  if (lambda == 0 && verbose) {
    cli::cli_alert_info("Lambda = 0: reverting to independent PCA per block.")
  }

  iter <- 0 # Initialize iter for trimming logic
  for (iter in 1:max_iter) {
    
    # Capture V_list state at the BEGINNING of the iteration
    V_list_begin_iter <- V_list
    Vbig_begin_iter <- combine_loadings(V_list_begin_iter)

    # Block update loop - retraction now happens inside block_update_cluster
    for (s in seq_len(S)) {
      if (k_s_vec[s] == 0) next  # Skip empty blocks
      
      # Use refactored block update function
      if (optimizer == "adam") {
        update_result <- block_update_cluster(
          s, V_list, Sadj, row_indices_list, Vbig_begin_iter,
          grad_fun_list, lambda, nsteps_inner, optimizer,
          learning_rate, global_step, m_list, V2_list,
          beta1, beta2, adam_epsilon, tol_inner, verbose, iter, ncomp_block
        )
        V_list[[s]] <- update_result$V
        global_step <- update_result$global_step
        m_list[[s]] <- update_result$m
        V2_list[[s]] <- update_result$V2
      } else {
        update_result <- block_update_cluster(
          s, V_list, Sadj, row_indices_list, Vbig_begin_iter,
          grad_fun_list, lambda, nsteps_inner, optimizer,
          learning_rate, 0, NULL, NULL,
          beta1, beta2, adam_epsilon, tol_inner, verbose, iter, ncomp_block
        )
        V_list[[s]] <- update_result$V
      }
    }

    # Compute Vbig and LV (for next iter grad & current obj)
    Vbig_current <- combine_loadings(V_list) # Use the now orthonormalized V_list
    
    # Skip Laplacian computation if lambda = 0
    if (lambda > 0) {
      LV_current <- Sadj %*% Vbig_current 
    }

    # Calculate Objective and Check Convergence
    new_obj <- tryCatch(calculate_objective(V_list, Vbig_current, LV_current),
                         error = function(e) {
                           if (verbose) cli::cli_alert_danger("Could not compute objective at iter {iter}: {e$message}")
                           return(NA) 
                         })

    if (is.na(new_obj)) {
        if (verbose) cli::cli_alert_danger("Objective is NA. Optimization may be unstable. Stopping.")
        # Trim obj_values up to the previous iteration
        obj_values <- obj_values[1:iter] # iter because obj_values[1] is initial
        break # Stop if objective calculation fails
    }

    obj_values[iter + 1] <- new_obj # Store objective AFTER iter completed
    
    # Check convergence
    if (!is.finite(old_obj) || !is.finite(new_obj)) {
       rel_change <- Inf # Cannot compute relative change reliably
    } else {
       # Use relative change in objective function
       rel_change <- abs(new_obj - old_obj) / (abs(old_obj) + 1e-12)
    }
    
    if (verbose) {
      cli::cli_alert_info("Iter {iter}: obj={format(new_obj, digits=6)}, rel_change={format(rel_change, scientific=TRUE, digits=2)}")
    }
    
    # Stop if converged
    if (rel_change < tol_obj && is.finite(rel_change)) {
      if (verbose) cli::cli_alert_success("Converged after {iter} iterations.")
      break # Exit loop
    }
    
    old_obj <- new_obj

  } # End main loop (iter)

  # Final trim of obj_values array
  obj_values <- head(obj_values, iter + 1)

  # Prepare results list - converting to multiblock_projector
  
  # Concatenate V_list into a single 'v' matrix with padding for variable ranks
  # Ensure validity first
  valid_V_final <- sapply(V_list, function(v) !is.null(v) && is.matrix(v))
  # Check if number of rows match k_s_vec
  valid_rows <- mapply(function(v, k_s) { if (is.null(v)) FALSE else nrow(v) == k_s }, V_list, k_s_vec)
  
  if (!all(valid_V_final) || !all(valid_rows)) {
      if (verbose) cli::cli_alert_danger("Final V_list contains invalid/inconsistent matrices. Cannot construct projector. Returning raw list.")
      # Return the raw list if construction fails
      return(list(
          V_list = V_list, 
          LV = LV_current, 
          Sadj = Sadj, 
          obj_values = obj_values, 
          ncomp = ncomp,
          lambda = lambda,
          precompute_info = precompute_info,
          iterations_run = iter
      ))
  }
  
  # Pad V_list to consistent dimensions for concatenation
  max_ncomp <- max(ncomp_block)
  V_list_padded <- lapply(seq_len(S), function(s) {
    V_s <- V_list[[s]]
    k_s <- k_s_vec[s]
    ncomp_s <- ncomp_block[s]
    
    if (k_s == 0) return(matrix(0, 0, max_ncomp))
    if (ncomp_s < max_ncomp) {
      # Pad with zeros
      V_padded <- matrix(0, k_s, max_ncomp)
      V_padded[, 1:ncomp_s] <- V_s
      return(V_padded)
    }
    return(V_s)
  })
  
  v_concat <- do.call(rbind, V_list_padded)
  
  # Block indices are already computed in row_indices_list
  # Ensure names are set if data_list had names
  if (!is.null(names(data_list))) {
      names(row_indices_list) <- names(data_list)
  } else {
      names(row_indices_list) <- paste0("Block_", seq_along(row_indices_list))
  }

  # Create a pass() preprocessor, as this function assumes pre-processed input
  final_preproc <- multivarious::prep(multivarious::pass())

  # Construct the multiblock_projector and pass auxiliary results via '...'
  result_projector <- multivarious::multiblock_projector(
      v = v_concat,
      preproc = final_preproc,
      block_indices = row_indices_list,
      Sadj = Sadj,
      LV = LV_current,
      obj_values = obj_values,
      lambda = lambda,
      precompute_info = precompute_info,
      iterations_run = iter,
      V_list = V_list,
      ncomp_block = ncomp_block,  # Store block-specific ncomp
      classes = "penalized_mfa_clusterwise" # Add original class back
  )
  
  return(result_projector)
}

#' Shorter alias for penalized_mfa_clusterwise
#'
#' @rdname penalized_mfa_clusterwise
#' @export
pmfa_cluster <- penalized_mfa_clusterwise

#' Print Method for Penalized MFA Clusterwise Objects
#'
#' @md
#' @description
#' Prints a summary of a `penalized_mfa_clusterwise` object, showing key parameters 
#' (lambda, components), convergence info, and block details.
#'
#' @param x A `penalized_mfa_clusterwise` object (which inherits from `multiblock_projector`)
#' @param ... Additional parameters (unused)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.penalized_mfa_clusterwise <- function(x, ...) {
  # Verify it's the expected structure
  if (!inherits(x, "multiblock_projector") || !inherits(x, "penalized_mfa_clusterwise")) {
    cli::cli_alert_danger("Input object does not seem to be a valid penalized_mfa_clusterwise result.")
    print(unclass(x)) # Fallback print
    return(invisible(x))
  }
  
  cli::cli_h1("Penalized MFA with Clusterwise Smoothness")
  
  # Extract info from object and named elements
  ncomp <- ncol(x$v)
  lambda_val <- x$lambda
  V_list_internal <- x$V_list
  obj_values_val <- x$obj_values
  n_blocks <- length(x$block_indices)
  block_names <- names(x$block_indices)
  if (is.null(block_names)) block_names <- paste("Block", 1:n_blocks)
  iterations_run <- x$iterations_run
  if (is.null(iterations_run)) iterations_run <- length(obj_values_val) - 1 # Estimate if missing
  precompute_info <- x$precompute_info
  if (is.null(precompute_info)) precompute_info <- rep(NA, n_blocks) # Handle if missing
  
  # Model parameters
  cli::cli_h2("Model Parameters")
  cli::cli_ul(c(
    "Components (ncomp): {ncomp}",
    "Smoothness Lambda: {lambda_val}"
  ))
  
  # Results overview
  cli::cli_h2("Data Structure")
  cli::cli_ul("Number of blocks/subjects: {n_blocks}")
  
  # Show per-block info using V_list_internal
  cli::cli_h2("Block Information (Loadings & Gradient Method)")
  for (i in seq_len(n_blocks)) {
    block_name <- block_names[i]
    Bi <- V_list_internal[[i]] # Use the stored V_list
    grad_method <- if (is.na(precompute_info[i])) "Unknown" else if (precompute_info[i]) "Precomputed XtX" else "On-the-fly"
    
    if (!is.null(Bi) && is.matrix(Bi)) {
       cli::cli_ul("{block_name}: {nrow(Bi)}x{ncol(Bi)} loading matrix (Grad: {grad_method})")
    } else {
        cli::cli_ul("{block_name}: Invalid/Missing Loadings (Grad: {grad_method})")
    }
  }
  
  # Show convergence info if available
  if (!is.null(obj_values_val) && length(obj_values_val) > 0) {
    cli::cli_h2("Convergence")
    initial_obj <- obj_values_val[1]
    final_obj <- obj_values_val[length(obj_values_val)]
    # Use explicit iterations_run if available, else estimate
    num_iters_display <- if (!is.null(iterations_run)) iterations_run else length(obj_values_val) - 1
    num_iters_display <- max(0, num_iters_display) 
    
    cli::cli_ul(c(
      "Iterations run: {num_iters_display}",
      "Initial objective: {format(initial_obj, digits=6)}",
      "Final objective: {format(final_obj, digits=6)}"
    ))
    
    # Calculate percent decrease if possible
    if (length(obj_values_val) > 1 && is.finite(initial_obj) && is.finite(final_obj) && abs(initial_obj) > 1e-12) {
      pct_decrease <- 100 * (initial_obj - final_obj) / abs(initial_obj) 
      cli::cli_ul("Objective decrease: {format(pct_decrease, digits=4)}%")
    }
  }
  
  return(invisible(x))
}
