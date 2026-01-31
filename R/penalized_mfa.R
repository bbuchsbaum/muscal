#' Perform a single Adam update step on a block's loadings.
#'
#' This function updates the loadings matrix \code{V} using one step of the Adam
#' optimizer. It maintains and returns the updated first (\code{M}) and second
#' (\code{V2}) moment estimates. The update is performed in the ambient space.
#'
#' @param V Current loading matrix (\eqn{p \times k}).
#' @param G The gradient of the loss with respect to V (\eqn{p \times k}).
#' @param M,V2 Current first and second moment estimate matrices.
#' @param step_count The global step number, used for bias correction.
#' @param beta1,beta2 Adam hyperparameters (e.g., 0.9, 0.999).
#' @param adam_epsilon A small constant to prevent division by zero (e.g., 1e-8).
#' @param learning_rate The base step size for the update.
#' @return A list containing the updated `V`, `M`, and `V2` matrices.
#' @keywords internal
adam_update_block <- function(V, G, M, V2, step_count,
                              beta1, beta2, adam_epsilon,
                              learning_rate) {
  # Update biased first and second moment estimates
  M <- beta1 * M + (1 - beta1) * G
  V2 <- beta2 * V2 + (1 - beta2) * (G * G)
  
  # Compute bias-corrected moment estimates
  M_hat <- M / (1 - beta1^step_count)
  V_hat <- V2 / (1 - beta2^step_count)
  
  # Compute the update step
  step_mat <- learning_rate * M_hat / (sqrt(V_hat) + adam_epsilon)
  
  # Apply the update in the ambient space
  V_new <- V - step_mat
  list(V = V_new, M = M, V2 = V2)
}

#' Prepare preprocessors and data for Penalized MFA
#'
#' @param data list of data blocks
#' @param preproc a preprocessor object or list of them
#' @param check_consistent_ncol whether to check for consistent number of columns post-preprocessing.
#' @return a list containing the preprocessed data `Xp`, the fitted preprocessors `proclist`,
#'         and the post-preprocessing dimension `p_post`.
#' @keywords internal
#' @noRd
prepare_block_preprocessors <- function(data, preproc, check_consistent_ncol = TRUE) {
  proclist <- if (is.null(preproc)) {
    # If no preproc, create a pass-through preprocessor for each block
    lapply(seq_along(data), function(i) multivarious::fit(multivarious::pass(), data[[i]]))
  } else if (inherits(preproc, "pre_processor") || inherits(preproc, "prepper")) {
    # If single preproc, create a fresh copy for each block
    lapply(seq_along(data), function(i) multivarious::fit(multivarious::fresh(preproc), data[[i]]))
  } else if (is.list(preproc)) {
    # If list of preprocs, prepare each one on its corresponding block
    chk::chk_equal(length(preproc), length(data))
    lapply(seq_along(data), function(i) {
      p <- preproc[[i]]
      if (!inherits(p, "pre_processor") && !inherits(p, "prepper")) {
        stop("Each element of 'preproc' list must be a pre_processor or prepper.")
      }
      multivarious::fit(multivarious::fresh(p), data[[i]])
    })
  } else {
    stop("`preproc` must be NULL, a pre_processor/prepper object, or a list of them.")
  }

  # Apply preprocessing
  Xp <- lapply(seq_along(data), function(i) multivarious::transform(proclist[[i]], data[[i]]))
  
  # Check for consistent dimensions after preprocessing
  dims_post <- sapply(Xp, dim)
  p_post <- dims_post[2, 1]
  
  if (check_consistent_ncol && !all(dims_post[2, ] == p_post)) {
    stop("All blocks must have the same number of columns after preprocessing.")
  }
  
  list(Xp = Xp, proclist = proclist, p_post = p_post)
}

#' @importFrom chk chk_numeric chk_list chk_integer chk_gte chk_lte chk_flag
#' @importFrom purrr map
#' @importFrom multivarious fresh prep init_transform center concat_pre_processors pass
#' @importFrom dplyr pull select
#' @importFrom rlang enquo
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_danger
#' @importFrom utils head tail
NULL

#' Penalized Multiple Factor Analysis (MFA)
#'
#' @description
#' This function implements a penalized Multiple Factor Analysis decomposition that
#' encourages similarity among block-specific loading matrices. The method uses
#' block-coordinate descent (BCD) with Riemannian optimization on the Stiefel manifold
#' to estimate orthonormal loading matrices for each data block while applying a
#' regularization penalty that promotes structural similarity across blocks.
#'
#' @details
#' ## Optimization Problem
#'
#' For \eqn{S} data blocks \eqn{\mathbf{X}_i \in \mathbb{R}^{n \times p_i}} (where
#' \eqn{i = 1, \ldots, S}), the method estimates loading matrices
#' \eqn{\mathbf{V}_i \in \mathbb{R}^{p_i \times k}} with orthonormal columns
#' (\eqn{\mathbf{V}_i^\top \mathbf{V}_i = \mathbf{I}_k}) by minimizing:
#'
#' \deqn{
#'   \min_{\{V_i\}} \sum_{i=1}^S \|\mathbf{X}_i - \mathbf{X}_i \mathbf{V}_i \mathbf{V}_i^\top\|_F^2
#'   \;+\; \lambda \,\text{Penalty}(\{\mathbf{V}_i\})
#' }
#'
#' The first term is the reconstruction error (sum of squared residuals) for each block,
#' and the second term is a penalty that encourages similarity among the loading matrices.
#'
#' ## Penalty Options
#'
#' Three penalty formulations are available:
#'
#' \describe{
#'   \item{\strong{Projection} (default, rotation-invariant):}{
#'     \deqn{\sum_{i=1}^S \|\mathbf{P}_i - \bar{\mathbf{P}}\|_F^2}
#'     where \eqn{\mathbf{P}_i = \mathbf{V}_i \mathbf{V}_i^\top} is the projection
#'     matrix for block \eqn{i} and \eqn{\bar{\mathbf{P}} = \frac{1}{S}\sum_{i=1}^S \mathbf{P}_i}.
#'     This penalty is invariant to rotations of the loading matrices and is recommended
#'     when the orientation of the latent space is arbitrary.
#'   }
#'   \item{\strong{Global Mean} (not rotation-invariant):}{
#'     \deqn{\sum_{i=1}^S \|\mathbf{V}_i - \bar{\mathbf{V}}\|_F^2}
#'     where \eqn{\bar{\mathbf{V}} = \frac{1}{S}\sum_{i=1}^S \mathbf{V}_i}. This directly
#'     penalizes Euclidean distance from the mean loading matrix. Computationally simpler
#'     than pairwise but not rotation-invariant.
#'   }
#'   \item{\strong{Pairwise} (not rotation-invariant):}{
#'     \deqn{\sum_{i < j} \|\mathbf{V}_i - \mathbf{V}_j\|_F^2}
#'     This penalizes all pairwise Euclidean distances between loading matrices.
#'     Equivalent to global mean up to a scaling factor but conceptually different.
#'   }
#' }
#'
#' ## Optimization Algorithm
#'
#' The algorithm uses block-coordinate descent with Riemannian gradient descent on
#' the Stiefel manifold:
#'
#' 1. **Outer loop** (max_iter iterations): Cycles through all blocks
#' 2. **Inner loop** (nsteps_inner steps per block): Updates each block's loading matrix
#' 3. **Gradient computation**: Combines reconstruction gradient and penalty gradient
#' 4. **Tangent space projection**: Projects gradient onto Stiefel manifold tangent space
#' 5. **Update step**: Uses either gradient descent or Adam optimizer
#' 6. **Retraction**: Re-orthonormalizes via QR decomposition to maintain manifold constraint
#'
#' The Stiefel manifold constraint ensures \eqn{\mathbf{V}_i^\top \mathbf{V}_i = \mathbf{I}_k}
#' at all iterations, which is critical for identifiability and interpretability.
#'
#' ## Convergence Criteria
#'
#' The outer loop stops when the relative change in the objective function falls below
#' `tol_obj`:
#' \deqn{\frac{|f^{(t+1)} - f^{(t)}|}{|f^{(t)}| + \epsilon} < \text{tol_obj}}
#'
#' Optionally, the inner loop can stop early if the Frobenius norm of the change in
#' \eqn{\mathbf{V}_i} falls below `tol_inner`.
#'
#' ## Preprocessing
#'
#' Data preprocessing is handled via the `multivarious` package. Common preprocessing
#' steps include centering (default) and scaling. Each block can have its own
#' preprocessor, or a single preprocessor can be applied to all blocks. Preprocessing
#' is "learned" from the input data and stored in the result object for later use on
#' new data.
#'
#' ## Consensus Loading Matrix
#'
#' When `compute_consensus = TRUE` and blocks have equal feature dimensions, a consensus
#' loading matrix is computed by orthonormalizing the sum of block-specific loadings:
#' \deqn{\mathbf{V}_{\text{consensus}} = \text{orth}\left(\sum_{i=1}^S \mathbf{V}_i\right)}
#'
#' This provides a single "average" loading matrix that summarizes the common structure.
#'
#' @section Input types:
#' The function supports three input types via S3 methods:
#' \describe{
#'   \item{\code{list}}{A list of numeric matrices, each representing a data block}
#'   \item{\code{multiblock}}{A multiblock object from the \code{multivarious} package}
#'   \item{\code{multidesign}}{A multidesign object where each subject is treated as a separate block}
#' }
#'
#' @param data A list of matrices, a `multiblock` object, or a `multidesign` object.
#'   For lists, each element should be a numeric matrix with the same number of rows
#'   (observations). Columns represent features and can differ across blocks.
#' @param ncomp Integer number of latent components to extract (default: 2). This is
#'   automatically capped at the minimum number of columns across all blocks after
#'   preprocessing. Must be at least 1.
#' @param lambda Non-negative scalar controlling the strength of the penalty (default: 1).
#'   Larger values encourage more similarity among block loadings. When `lambda = 0`,
#'   the method reduces to independent PCA on each block. Typical range: 0 to 10.
#' @param penalty_method Character string specifying the penalty type (default: "projection").
#'   Options are:
#'   \itemize{
#'     \item \code{"projection"}: Rotation-invariant penalty on projection matrices (recommended)
#'     \item \code{"pairwise"}: Pairwise Euclidean distances between loading matrices
#'     \item \code{"global_mean"}: Distance from mean loading matrix
#'   }
#' @param max_iter Maximum number of outer block-coordinate descent iterations (default: 10).
#'   Typical range: 10 to 100. More iterations allow better convergence but increase
#'   computation time.
#' @param nsteps_inner Number of gradient update steps per block in each outer iteration
#'   (default: 5). Higher values perform more thorough optimization of each block before
#'   moving to the next. Typical range: 1 to 20.
#' @param learning_rate Step size for the optimizer (default: 0.01). Controls how far
#'   to move in the gradient direction. Too large may cause divergence; too small slows
#'   convergence. Typical range for gradient descent: 0.001 to 0.1. Adam is less
#'   sensitive to this parameter.
#' @param optimizer Character string specifying the optimization algorithm (default: "adam").
#'   Options are:
#'   \itemize{
#'     \item \code{"adam"}: Adaptive Moment Estimation with momentum and adaptive learning rates
#'     \item \code{"gradient"}: Standard gradient descent with fixed learning rate
#'   }
#'   Adam is generally more robust and converges faster with default settings.
#' @param preproc A preprocessing specification (default: `multivarious::center()`).
#'   Can be:
#'   \itemize{
#'     \item A single `pre_processor` object applied to all blocks
#'     \item A list of `pre_processor` objects (one per block)
#'     \item \code{NULL} for no preprocessing (identity transformation)
#'   }
#'   Common options: `center()`, `standardize()`, `pass()` (no-op).
#' @param beta1 Adam hyperparameter for first moment decay (default: 0.9). Controls
#'   the exponential decay rate for gradient moving average. Typical range: 0.8 to 0.95.
#'   Only used when `optimizer = "adam"`.
#' @param beta2 Adam hyperparameter for second moment decay (default: 0.999). Controls
#'   the exponential decay rate for squared gradient moving average. Typical range:
#'   0.99 to 0.9999. Only used when `optimizer = "adam"`.
#' @param adam_epsilon Small constant added to denominator in Adam for numerical
#'   stability (default: 1e-8). Prevents division by zero. Only used when
#'   `optimizer = "adam"`.
#' @param tol_obj Numeric tolerance for outer loop convergence based on relative
#'   change in objective function (default: 1e-7). Smaller values require more precise
#'   convergence. Typical range: 1e-8 to 1e-5.
#' @param tol_inner Optional numeric tolerance for inner loop convergence based on
#'   Frobenius norm of change in loading matrix (default: NULL, no early stopping).
#'   When specified, inner loop stops if \eqn{\|\mathbf{V}_i^{\text{new}} - \mathbf{V}_i\|_F < \text{tol_inner}}.
#' @param compute_consensus Logical indicating whether to compute a consensus loading
#'   matrix (default: FALSE). Only computed when all blocks have the same number of
#'   features after preprocessing. The consensus is the orthonormalized sum of
#'   block-specific loadings.
#' @param verbose Logical indicating whether to print iteration progress (default: FALSE).
#'   When TRUE, displays objective values and relative changes at each iteration using
#'   the `cli` package.
#' @param subject Required for `multidesign` method: the name (unquoted) of the subject
#'   variable used to split data into blocks. Each subject's data becomes a separate
#'   block in the analysis.
#' @param ... Additional arguments passed to methods. Currently unused but included
#'   for S3 method consistency.
#'
#' @return A `multiblock_projector` object of class `penalized_mfa` with the following components:
#'   \describe{
#'     \item{\code{v}}{Concatenated loading matrix (total_features × ncomp) formed by
#'       vertically stacking all block-specific loading matrices.}
#'     \item{\code{preproc}}{A `concat_pre_processor` object containing the preprocessing
#'       transformations applied to each block. Can be used to transform new data.}
#'     \item{\code{block_indices}}{A named list indicating which rows of `v` correspond
#'       to each block. Each element is an integer vector of row indices.}
#'   }
#'
#'   Additional information is stored as attributes:
#'   \describe{
#'     \item{\code{V_list}}{List of length S containing the final orthonormal loading
#'       matrices for each block. Each element is a (p_i × ncomp) matrix.}
#'     \item{\code{obj_values}}{Numeric vector of objective function values at each
#'       iteration, including the initial value. Length is (iterations_run + 1).}
#'     \item{\code{consensus}}{The consensus loading matrix (p\_post × ncomp) if
#'       `compute_consensus=TRUE`, otherwise `NULL`. This is the orthonormalized
#'       average of block loadings over blocks with matching feature dimension.}
#'     \item{\code{lambda}}{The penalty strength used in the optimization.}
#'     \item{\code{penalty_method}}{Character string indicating the penalty type used.
#'       May be "none" if penalty was disabled due to inconsistent block dimensions.}
#'     \item{\code{iterations_run}}{Integer indicating how many outer iterations were
#'       completed before convergence or reaching max_iter.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic usage with simulated data
#' set.seed(123)
#' data_list <- lapply(1:3, function(i) {
#'   matrix(rnorm(100), 10, 10)
#' })
#' res <- penalized_mfa(data_list, ncomp=2, lambda=1, penalty_method="projection",
#'                      optimizer="adam", max_iter=50, verbose=TRUE)
#' print(res)
#'
#' # Access block-specific loadings
#' V_list <- attr(res, "V_list")
#' print(V_list[[1]])  # Loadings for first block
#'
#' # Plot convergence
#' obj_vals <- attr(res, "obj_values")
#' plot(obj_vals, type='b', xlab='Iteration', ylab='Objective',
#'      main='Convergence of Penalized MFA')
#'
#' # Example 2: With consensus and custom preprocessing
#' library(multivarious)
#' res2 <- penalized_mfa(data_list, ncomp=3, lambda=2,
#'                       preproc=standardize(),  # Center and scale
#'                       compute_consensus=TRUE,
#'                       verbose=FALSE)
#' consensus_loadings <- attr(res2, "consensus")
#'
#' # Example 3: Comparing penalty methods
#' res_proj <- penalized_mfa(data_list, ncomp=2, lambda=1, penalty_method="projection")
#' res_mean <- penalized_mfa(data_list, ncomp=2, lambda=1, penalty_method="global_mean")
#' res_pair <- penalized_mfa(data_list, ncomp=2, lambda=1, penalty_method="pairwise")
#'
#' # Example 4: Lambda selection via objective values
#' lambdas <- c(0, 0.1, 0.5, 1, 2, 5, 10)
#' results <- lapply(lambdas, function(lam) {
#'   fit <- penalized_mfa(data_list, ncomp=2, lambda=lam, verbose=FALSE)
#'   list(lambda=lam, final_obj=tail(attr(fit, "obj_values"), 1))
#' })
#'
#' # Example 5: Using with multiblock object
#' mb <- multiblock(data_list)
#' res_mb <- penalized_mfa(mb, ncomp=2, lambda=1)
#'
#' # Example 6: Different preprocessors per block
#' preproc_list <- list(
#'   center(),
#'   standardize(),
#'   pass()  # No preprocessing for third block
#' )
#' res_custom <- penalized_mfa(data_list, ncomp=2, lambda=1,
#'                             preproc=preproc_list)
#'
#' # Example 7: Gradient descent instead of Adam
#' res_gd <- penalized_mfa(data_list, ncomp=2, lambda=1,
#'                         optimizer="gradient",
#'                         learning_rate=0.001,
#'                         max_iter=100)
#' }
#'
#' @export
penalized_mfa <- function(data, ...) {
  UseMethod("penalized_mfa")
}

#' @rdname penalized_mfa
#' @export
penalized_mfa.list <- function(data,
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
                               adam_epsilon = 1e-8,
                               tol_obj = 1e-7,
                               tol_inner = NULL,
                               compute_consensus = FALSE,
                               verbose = FALSE,
                               ...
                               ) {
  
  penalty_method <- match.arg(penalty_method)
  optimizer <- match.arg(optimizer)
  
  chk::chk_list(data)
  if (length(data) < 2) stop("penalized_mfa requires at least 2 data blocks.")
  
  prep_res <- prepare_block_preprocessors(data, preproc, check_consistent_ncol = FALSE)
  Xp <- prep_res$Xp
  proclist <- prep_res$proclist
  p_dims <- sapply(Xp, ncol)
  ncomp <- min(ncomp, min(p_dims))
  if (ncomp < 1) stop("ncomp must be at least 1 after preprocessing.")
  
  engine_res <- penalized_mfa_engine(
    Xp, ncomp, lambda, penalty_method, max_iter, nsteps_inner,
    learning_rate, optimizer, beta1, beta2, adam_epsilon,
    tol_obj, tol_inner, verbose, compute_consensus
  )
  
  V_list <- engine_res$V_list
  v_concat <- do.call(rbind, V_list)

  # Per-block scores in the learned subspaces (n x k each); define a simple
  # "consensus" score matrix as the mean across blocks for plotting.
  scores_list <- lapply(seq_along(Xp), function(i) Xp[[i]] %*% V_list[[i]])
  names(scores_list) <- names(Xp)
  s_consensus <- Reduce(`+`, scores_list) / length(scores_list)
  sdev <- apply(s_consensus, 2, stats::sd)

  cor_loadings <- tryCatch(
    {
      Xp_concat <- do.call(cbind, Xp)
      C <- stats::cor(Xp_concat, s_consensus)
      C[!is.finite(C)] <- 0
      C
    },
    error = function(e) NULL
  )
  
  block_indices <- list()
  current_start <- 1
  for(i in 1:length(Xp)) {
      p_i <- ncol(Xp[[i]])
      block_indices[[i]] <- current_start:(current_start + p_i - 1)
      current_start <- current_start + p_i
  }
  names(block_indices) <- names(Xp)
  
  final_preproc <- multivarious::concat_pre_processors(proclist, block_indices)
  
  result_projector <- multivarious::multiblock_biprojector(
    v = v_concat,
    s = s_consensus,
    sdev = sdev,
    preproc = final_preproc,
    block_indices = block_indices,
    # extra fields for plotting / diagnostics
    V_list = V_list,
    scores_list = scores_list,
    objective_trace = engine_res$obj_values,
    cor_loadings = cor_loadings,
    classes = "penalized_mfa"
  )
  
  attr(result_projector, "V_list") <- V_list
  attr(result_projector, "obj_values") <- engine_res$obj_values
  attr(result_projector, "lambda") <- lambda
  attr(result_projector, "penalty_method") <- attr(engine_res, "penalty_method_used")
  attr(result_projector, "iterations_run") <- engine_res$iterations_run
  attr(result_projector, "consensus") <- engine_res$consensus
  
  result_projector
}

#' @rdname penalized_mfa
#' @export
penalized_mfa.multiblock <- function(data, ...) {
  # The list method can handle the multiblock object by converting it.
  penalized_mfa.list(as.list(data), ...)
}

#' @rdname penalized_mfa
#' @export
penalized_mfa.multidesign <- function(data, subject, ...) {
  # This method now treats each subject as a "block" in the penalized analysis.
  # This is a change from the original implementation which iterated over subjects.
  # The penalty is now applied across subjects' data.
  
  # Check for multidesign package
  if (!requireNamespace("multidesign", quietly = TRUE)) {
    stop("Package 'multidesign' is required for this method.")
  }
  
  subject_quo <- rlang::enquo(subject)
  
  # `split` creates a list of multidesign objects, one per subject.
  sdat_multiblocks <- split(data, !!subject_quo)
  
  # For this cross-subject analysis, we extract the full data matrix from each subject.
  # This assumes each subject's multiblock should be treated as a single observation matrix.
  sdat_matrices <- lapply(sdat_multiblocks, function(mb) mb$x)
  names(sdat_matrices) <- names(sdat_multiblocks)
  
  # Now, call the list method on the list of subject-level matrices.
  penalized_mfa.list(sdat_matrices, ...)
}

#' Print Method for Penalized MFA Objects
#'
#' @param x A `penalized_mfa` object
#' @param ... Additional parameters (unused)
#' @return Invisibly returns the input object
#' @export
print.penalized_mfa <- function(x, ...) {
  is_projector <- inherits(x, "multiblock_projector")
  
  cli::cli_h1("Penalized Multiple Factor Analysis (MFA)")
  
  if (!is_projector) {
    cli::cli_alert_danger("Object is not a valid 'penalized_mfa' projector. Printing basic info.")
    print(unclass(x))
    return(invisible(x))
  }
  
  lambda_val <- attr(x, "lambda")
  penalty_method_val <- attr(x, "penalty_method")
  V_list_internal <- attr(x, "V_list")
  obj_values_val <- attr(x, "obj_values")
  consensus_val <- attr(x, "consensus")
  n_blocks <- length(x$block_indices)
  iterations_run <- attr(x, "iterations_run")
  
  cli::cli_h2("Model Parameters")
  cli::cli_ul(c(
    "Components (k): {ncol(V_list_internal[[1]])}",
    "Penalty type: {penalty_method_val}",
    "Lambda: {lambda_val}"
  ))
  
  cli::cli_h2("Results")
  cli::cli_ul(c(
    "Number of blocks: {n_blocks}",
    "Iterations run: {iterations_run}"
  ))
  
  if (!is.null(obj_values_val) && length(obj_values_val) > 1) {
    initial_obj <- obj_values_val[1]
    final_obj <- tail(obj_values_val, 1)
    
    cli::cli_h2("Convergence")
    cli::cli_ul(c(
      "Initial objective: {round(initial_obj, 4)}",
      "Final objective: {round(final_obj, 4)}"
    ))
  }
  
  if (!is.null(consensus_val)) {
    cli::cli_h2("Consensus Loadings")
    cli::cli_alert_info("A consensus loading matrix ({nrow(consensus_val)} x {ncol(consensus_val)}) is available in `attr(x, 'consensus')`.")
  }
  
  invisible(x)
}

#' Core engine for penalized MFA optimization
#' @keywords internal
penalized_mfa_engine <- function(Xp,
                                 ncomp,
                                 lambda,
                                 penalty_method,
                                 max_iter,
                                 nsteps_inner,
                                 learning_rate,
                                 optimizer,
                                 beta1,
                                 beta2,
                                 adam_epsilon,
                                 tol_obj,
                                 tol_inner,
                                 verbose,
                                 compute_consensus) {

  S <- length(Xp)
  n_dims <- sapply(Xp, nrow)
  p_dims <- sapply(Xp, ncol)
  consistent_p <- all(p_dims == p_dims[1])
  use_precomputed_xtx <- all(n_dims >= p_dims)
  
  penalty_method_used <- penalty_method
  if (!consistent_p && penalty_method == "projection") {
    warning("Cannot use 'projection' penalty with blocks of different dimensions; disabling penalty (lambda set to 0).")
    penalty_method_used <- "none"
    lambda <- 0
  } else if (!consistent_p && penalty_method %in% c("pairwise", "global_mean")) {
    warning("Penalty requires blocks with equal feature dimension; disabling penalty.")
    penalty_method_used <- "none"
    lambda <- 0
  }

  # Initialize loadings V_i
  V_list <- lapply(Xp, function(Xi) {
    k <- min(ncomp, ncol(Xi))
    if (k < 1) stop("Each block must have at least one column after preprocessing.")
    sv <- svd(Xi, nu = 0, nv = k)
    sv$v[, 1:k, drop = FALSE]
  })
  
  # Initialize Adam moments if needed
  if (optimizer == "adam") {
    m_list <- lapply(V_list, function(V) matrix(0, nrow(V), ncol(V)))
    v_list <- lapply(V_list, function(V) matrix(0, nrow(V), ncol(V)))
    step_counts <- integer(S)
  }
  
  # Pre-compute XtX for each block when beneficial (n >= p)
  XtX_list <- if (use_precomputed_xtx) lapply(Xp, crossprod) else vector("list", S)
  norm_X2_vec <- vapply(Xp, function(x) sum(x * x), numeric(1))

  # Objective function
  obj_fun <- function(Vs) {
    recon_cost <- 0
    for (i in seq_len(S)) {
      XV <- Xp[[i]] %*% Vs[[i]]
      recon_cost <- recon_cost + norm_X2_vec[i] - sum(XV * XV)
    }
    
    pen_cost <- if (lambda > 0 && penalty_method_used != "none") {
      if (penalty_method_used == "projection") {
        Ps <- lapply(Vs, function(V) V %*% t(V))
        P_bar <- Reduce(`+`, Ps) / S
        lambda * sum(sapply(Ps, function(P) sum((P - P_bar)^2)))
      } else { # pairwise or global_mean
        V_bar <- Reduce(`+`, Vs) / S
        lambda * S * sum(sapply(Vs, function(V) sum((V - V_bar)^2)))
      }
    } else { 0 }
    recon_cost + pen_cost
  }
  
  # Main optimization loop
  obj_values <- numeric(max_iter + 1)
  obj_values[1] <- obj_fun(V_list)
  
  if (verbose) cli::cli_alert_info("Starting optimization... Initial objective: {round(obj_values[1], 4)}")

  for (iter in 1:max_iter) {
    
    # Pre-compute means for the iteration only if dimensions are consistent
    V_bar <- if (consistent_p) Reduce(`+`, V_list) / S else NULL
    P_bar <- if (consistent_p && penalty_method_used == "projection") {
      Ps_iter <- lapply(V_list, function(V) V %*% t(V))
      Reduce(`+`, Ps_iter) / S
    } else { NULL }
    
    for (b in 1:S) { # Iterate through blocks
      
      Vi <- V_list[[b]]
      
      if (optimizer == "adam") {
        Mi <- m_list[[b]]
        V2 <- v_list[[b]]
      }
      
      for (step_inner in 1:nsteps_inner) {
        
        if (use_precomputed_xtx) {
          grad_recon <- -2 * XtX_list[[b]] %*% Vi
        } else {
          XV_block <- Xp[[b]] %*% Vi
          grad_recon <- -2 * crossprod(Xp[[b]], XV_block)
        }
        
        grad_penalty <- if (lambda > 0 && penalty_method_used != "none") {
          if (penalty_method_used == "projection") {
            Pi <- Vi %*% t(Vi)
            4 * lambda * (Pi - P_bar) %*% Vi
          } else { # pairwise or global_mean
            2 * lambda * S * (Vi - V_bar)
          }
        } else { 0 }
        
        G <- grad_recon + grad_penalty
        
        # Project gradient to tangent space of Stiefel manifold
        A <- t(Vi) %*% G
        G_tangent <- G - Vi %*% ( (A + t(A)) / 2 )
        
        # Update step
        if (optimizer == "adam") {
          step_counts[b] <- step_counts[b] + 1
          adam_res <- adam_update_block(Vi, G_tangent, Mi, V2, step_counts[b], beta1, beta2, adam_epsilon, learning_rate)
          Vi_new <- adam_res$V
          Mi <- adam_res$M
          V2 <- adam_res$V2
        } else { # gradient descent
          Vi_new <- Vi - learning_rate * G_tangent
        }
        
        # Retraction: Re-orthonormalize to stay on the manifold
        Vi_new <- qr.Q(qr(Vi_new))
        
        # Check for inner loop convergence
        if (!is.null(tol_inner) && norm(Vi_new - Vi, "F") < tol_inner) {
          Vi <- Vi_new
          break
        }
        Vi <- Vi_new
      } # End inner loop
      
      V_list[[b]] <- Vi
      if (optimizer == "adam") {
        m_list[[b]] <- Mi
        v_list[[b]] <- V2
      }
    } # End block loop
    
    # Calculate and store objective
    obj_values[iter + 1] <- obj_fun(V_list)
    
    # Check for outer loop convergence
    rel_change <- abs(obj_values[iter + 1] - obj_values[iter]) / (abs(obj_values[iter]) + 1e-12)
    
    if (verbose) {
      cli::cli_alert_info("Iter {iter}: obj={round(obj_values[iter+1], 4)}, rel_change={format(rel_change, scientific=TRUE, digits=2)}")
    }
    
    if (rel_change < tol_obj) {
      if (verbose) cli::cli_alert_success("Converged after {iter} iterations.")
      break
    }
  } # End main loop
  
  # Compute consensus only if dimensions are consistent
  consensus_v <- NULL
  if (compute_consensus && consistent_p) {
    if (verbose) cli::cli_alert_info("Computing consensus loading matrix via simple averaging.")
    Vsum <- Reduce(`+`, V_list)
    consensus_v <- qr.Q(qr(Vsum))[, 1:ncomp, drop = FALSE]
  } else if (compute_consensus && !consistent_p) {
      warning("Cannot compute consensus loading matrix because blocks have different feature dimensions.")
  }

  
  res <- list(
    V_list = V_list,
    obj_values = head(obj_values, iter + 1),
    iterations_run = iter,
    consensus = consensus_v
  )
  attr(res, "penalty_method_used") <- penalty_method_used
  res
}
