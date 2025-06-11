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
prepare_block_preprocessors <- function(data, preproc, check_consistent_ncol = TRUE) {
  proclist <- if (is.null(preproc)) {
    # If no preproc, create a pass-through preprocessor for each block
    lapply(seq_along(data), function(i) multivarious::prep(multivarious::pass()))
  } else if (inherits(preproc, "pre_processor")) {
    # If single preproc, create a fresh copy for each block
    lapply(seq_along(data), function(i) multivarious::fresh(preproc) %>% multivarious::prep(data[[i]]))
  } else if (is.list(preproc)) {
    # If list of preprocs, use them directly
    chk::chk_equal(length(preproc), length(data))
    preproc
  } else {
    stop("`preproc` must be NULL, a pre_processor object, or a list of pre_processors.")
  }
  
  # Apply preprocessing
  Xp <- multivarious::apply_transform(proclist, data)
  
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
#' This function implements a penalized MFA-like decomposition via block-coordinate
#' descent (BCD). Each block \eqn{\mathbf{X}_i \in \mathbb{R}^{n \times p_i}} has a loading
#' matrix \eqn{\mathbf{V}_i \in \mathbb{R}^{p_i \times k}} with orthonormal columns.
#' It assumes data blocks `X_i` are pre-processed (e.g., column-centered) as needed.
#'
#' We solve:
#' \deqn{
#'   \min_{\{V_i\}} \sum_{i=1}^S \|X_i - X_i V_i V_i^\top\|_F^2
#'   \;+\; \lambda \,\mathrm{Penalty}(\{V_i\}),
#' }
#' where the penalty is one of:
#' \itemize{
#'   \item \code{"pairwise"}: \deqn{\sum_{i < j} \|V_i - V_j\|_F^2}. Penalizes Euclidean distance between loading matrices. Not rotation-invariant.
#'   \item \code{"global_mean"}: \deqn{\sum_{i} \|V_i - \bar{V}\|_F^2}. Similar to "pairwise" but computationally simpler. Not rotation-invariant.
#'   \item \code{"projection"}: \deqn{\sum_{i} \|P_i - \bar{P}\|_F^2}, where \eqn{P_i = V_i V_i^\top} is the projection matrix for block `i`. This penalty is rotation-invariant.
#' }
#'
#' @section Input types:
#' The function supports several input types for `data`: `list`, `multiblock`, and `multidesign`.
#'
#' @param data A list of matrices, a `multiblock` object, or a `multidesign` object.
#' @param ncomp Integer number of latent components.
#' @param lambda Non-negative scalar controlling the penalty strength.
#' @param penalty_method One of `"pairwise"`, `"global_mean"`, or `"projection"`.
#' @param max_iter Maximum number of outer BCD iterations.
#' @param nsteps_inner Number of gradient updates per block in each iteration.
#' @param learning_rate Step size for the optimizer.
#' @param optimizer One of `"gradient"` (fixed step size) or `"adam"` (adaptive).
#' @param preproc A `pre_processor` object from `multivarious` (e.g., `center()`), a list of such objects, or `NULL` for no preprocessing.
#' @param beta1,beta2 Adam hyperparameters.
#' @param adam_epsilon Small constant for Adam's denominator.
#' @param tol_obj Numeric tolerance for the relative change in the objective function to determine convergence.
#' @param tol_inner Optional tolerance for stopping inner loop updates based on the norm of the change in `V`.
#' @param compute_consensus If `TRUE`, computes a consensus (average) loading matrix. Uses Procrustes alignment via the `vegan` package if installed.
#' @param verbose If `TRUE`, prints iteration logs using the `cli` package.
#' @param subject Required for `multidesign` method: the name of the subject variable.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `multiblock_projector` object of class `penalized_mfa`, containing the concatenated loadings (`v`), block indices, and preprocessor.
#'   Additional results are stored as attributes:
#'   \item{V_list}{The list of final, orthonormal loading matrices for each block.}
#'   \item{obj_values}{A vector of objective function values at each iteration.}
#'   \item{consensus}{The consensus loading matrix, if `compute_consensus=TRUE`.}
#'   \item{lambda, penalty_method, iterations_run}{Model and convergence information.}
#'
#' @examples
#' \dontrun{
#' # Example with a list of matrices
#' data_list <- lapply(1:3, function(i) scale(matrix(rnorm(100), 10, 10), scale=FALSE))
#' res <- penalized_mfa(data_list, ncomp=2, lambda=1, penalty_method="projection",
#'                      optimizer="adam", max_iter=50, verbose=TRUE)
#'
#' # Example with a multiblock object
#' mb <- multiblock(data_list)
#' res_mb <- penalized_mfa(mb, ncomp=2, lambda=1)
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
  
  engine_res <- penalized_mfa_engine(
    Xp, ncomp, lambda, penalty_method, max_iter, nsteps_inner,
    learning_rate, optimizer, beta1, beta2, adam_epsilon,
    tol_obj, tol_inner, verbose, compute_consensus
  )
  
  V_list <- engine_res$V_list
  v_concat <- do.call(rbind, V_list)
  
  block_indices <- list()
  current_start <- 1
  for(i in 1:length(Xp)) {
      p_i <- ncol(Xp[[i]])
      block_indices[[i]] <- current_start:(current_start + p_i - 1)
      current_start <- current_start + p_i
  }
  names(block_indices) <- names(Xp)
  
  final_preproc <- multivarious::concat_pre_processors(proclist, block_indices)
  
  result_projector <- multivarious::multiblock_projector(
      v = v_concat,
      preproc = final_preproc,
      block_indices = block_indices,
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
  
  # `split_by` creates a list of multiblock objects, one per subject.
  sdat_multiblocks <- multidesign::split_by(data, !!subject_quo)
  
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
  p_dims <- sapply(Xp, ncol)
  consistent_p <- all(p_dims == p_dims[1])
  
  # Validate penalty method based on dimension consistency
  penalty_method_used <- penalty_method
  if (penalty_method == "projection" && !consistent_p) {
    penalty_method_used <- "pairwise" # Fallback
    warning("Cannot use 'projection' penalty with blocks of different dimensions. Falling back to 'pairwise'.")
  }

  # Initialize loadings V_i
  V_list <- lapply(Xp, function(Xi) {
    if (ncol(Xi) < ncomp) {
       warning("Block has fewer columns than ncomp. Cannot initialize with SVD.", call. = FALSE)
       return(matrix(0, ncol(Xi), ncomp))
    }
    res <- if (nrow(Xi) < ncol(Xi)) svd(crossprod(Xi), nu = ncomp) else svd(Xi, nu = 0, nv = ncomp)
    res$v[, 1:ncomp, drop = FALSE]
  })
  
  # Initialize Adam moments if needed
  if (optimizer == "adam") {
    m_list <- lapply(V_list, function(V) matrix(0, nrow(V), ncol(V)))
    v_list <- lapply(V_list, function(V) matrix(0, nrow(V), ncol(V)))
    global_step <- 0
  }
  
  # Pre-compute XtX for each block
  XtX_list <- lapply(Xp, crossprod)

  # Objective function
  obj_fun <- function(Vs, Ps = NULL) {
    recon_cost <- sum(sapply(1:S, function(i) {
      resid <- Xp[[i]] - (Xp[[i]] %*% Vs[[i]]) %*% t(Vs[[i]])
      sum(resid^2)
    }))
    
    pen_cost <- if (lambda > 0 && consistent_p) { # Only apply penalty if dimensions are consistent
      if (penalty_method_used == "projection") {
        if (is.null(Ps)) Ps <- lapply(Vs, function(V) V %*% t(V))
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
      lapply(V_list, function(V) V %*% t(V)) %>% Reduce(`+`, .) / S
    } else { NULL }
    
    for (b in 1:S) { # Iterate through blocks
      
      Vi <- V_list[[b]]
      
      if (optimizer == "adam") {
        Mi <- m_list[[b]]
        V2 <- v_list[[b]]
      }
      
      for (step_inner in 1:nsteps_inner) {
        
        grad_recon <- -2 * XtX_list[[b]] %*% Vi
        
        # Only compute penalty gradient if dimensions are consistent
        grad_penalty <- if (lambda > 0 && consistent_p) {
          if (penalty_method_used == "projection") {
            Pi <- Vi %*% t(Vi)
            2 * lambda * (Pi - P_bar) %*% Vi
          } else { # pairwise or global_mean
            2 * lambda * (Vi - V_bar)
          }
        } else { 0 }
        
        G <- grad_recon + grad_penalty
        
        # Project gradient to tangent space of Stiefel manifold
        A <- t(Vi) %*% G
        G_tangent <- G - Vi %*% ( (A + t(A)) / 2 )
        
        # Update step
        if (optimizer == "adam") {
          global_step <- global_step + 1
          adam_res <- adam_update_block(Vi, G_tangent, Mi, V2, global_step, beta1, beta2, adam_epsilon, learning_rate)
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
    if (verbose) cli::cli_alert_info("Using simple averaging for consensus matrix.")
    Vsum <- Reduce(`+`, V_list)
    consensus_v <- qr.Q(qr(Vsum))[, 1:ncomp, drop = FALSE]
  } else if (compute_consensus && !consistent_p) {
      warning("Cannot compute consensus loading matrix because blocks have different numbers of features.")
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
