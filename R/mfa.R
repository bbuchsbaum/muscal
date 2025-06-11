#' Compute a similarity matrix from blocks of data
#'
#' Creates a symmetric matrix where each element [i,j] is the similarity between
#' blocks i and j, calculated using the supplied function.
#'
#' @param blocks A list of numeric matrices or data frames
#' @param FUN Function to compute similarity between two blocks
#' @param ... Additional arguments passed to FUN
#' @return A symmetric similarity matrix with dimensions length(blocks) × length(blocks)
#' @importFrom utils combn
#' @noRd
#' @keywords internal
compute_sim_mat <- function(blocks, FUN, ...) {
  pairs <- combn(length(blocks),2)
  M <- matrix(0, length(blocks), length(blocks))
  
  for (i in 1:ncol(pairs)) {
    p1 <- pairs[1,i]
    p2 <- pairs[2,i]
    sim <- FUN(blocks[[p1]], blocks[[p2]], ...)
    M[p1,p2] <- sim
    M[p2,p1] <- sim
  }
  
  M
}

#' Calculate normalization factors for blocks in MFA
#'
#' Determines weighting factors for each block depending on the selected normalization method.
#'
#' @param blocks A list of preprocessed data blocks
#' @param type The normalization method to use: "MFA", "RV", "RV2", "None", or "Frob"
#' @return A numeric vector of normalization factors, one per block
#' @noRd
#' @keywords internal
normalization_factors <- function(blocks, type=c("MFA", "RV", "RV2", "None", "Frob")) {
  type <- match.arg(type)
  message("normalization type:", type)
  alpha <- if (type == "MFA") {
    unlist(lapply(blocks, function(X) 1/(multivarious::svd_wrapper(X, ncomp=1, method="svds")$sdev[1]^2)))
  } else if (type == "RV" && length(blocks) > 2) {
    smat <- compute_sim_mat(blocks, function(x1,x2) MatrixCorrelation::RV2(x1,x2))
    diag(smat) <- 1
    abs(multivarious::svd_wrapper(smat, ncomp=1, method="svds")$u[,1])
  } else if (type == "RV2" && length(blocks) > 2) {
    smat <- compute_sim_mat(blocks, function(x1,x2) MatrixCorrelation::RV(x1,x2))
    diag(smat) <- 1
    abs(multivarious::svd_wrapper(smat, ncomp=1, method="svds")$u[,1])
  } else if (type == "Frob") {
    unlist(lapply(as.list(blocks), function(X) sum(X^2)))
  } else {
    rep(1, length(blocks))
  }
}


#' @md
#' @rdname mfa
#' @details
#' The `mfa.list` method applies the MFA algorithm to a list of data matrices or data frames.
#' This method first converts the list to a multiblock object and then calls `mfa.multiblock`.
#'
#' @examples
#' # Apply MFA to a list of matrices
#' X <- replicate(5, { matrix(rnorm(10*10), 10, 10) }, simplify=FALSE)
#' res <- mfa(X, ncomp=3, normalization="MFA")
#' 
#' @export
mfa.list <- function(data, preproc=center(), ncomp=2,
                     normalization=c("MFA", "RV", "None", "Frob", "custom"),
                     A=NULL, M=NULL, ...) {
  data <- multiblock(data)
  mfa.multiblock(data, preproc, ncomp, normalization, A, M,...)
}


#' @md
#' @rdname mfa
#' @details
#' The `mfa.multiblock` method implements Multiple Factor Analysis for a collection of 
#' data blocks. This method handles data preprocessing, block normalization, and integration 
#' of multiple data tables that share the same observations.
#' 
#' Normalization options include:
#' * `MFA`: Scales each block by its first singular value (default)
#' * `RV`: Normalizes blocks based on RV matrix correlation
#' * `None`: No scaling applied
#' * `Frob`: Uses Frobenius norm for scaling
#' * `custom`: Uses custom weight matrices provided via A and M parameters
#'
#' @examples 
#' # Create 5 random matrices of the same size
#' X <- replicate(5, { matrix(rnorm(10*10), 10, 10) }, simplify=FALSE)
#' 
#' # Apply MFA with MFA normalization
#' res <- mfa(X, ncomp=3, normalization="MFA")
#' 
#' # Project a block onto the model
#' p <- multivarious::project_block(res, X[[1]], 1)
#' 
#' # Verify number of components
#' stopifnot(ncol(multivarious::scores(res)) == 3)
#' 
#' # Create a classifier
#' labs <- letters[1:10]
#' cfier <- multivarious::classifier(res, new_data=do.call(cbind, X), labels=labs)
#' pred <- predict(cfier, X[1:2,])
#' 
#' # Create a classifier using a specific block
#' cfier2 <- multivarious::classifier(res, new_data=X[[2]], labels=labs, 
#'                                   colind=res$block_indices[[2]])
#' pred2 <- predict(cfier2, X[1:2,res$block_indices[[2]]])
#' @export
mfa.multiblock <- function(data, preproc=center(), ncomp=2,
                normalization=c("MFA", "RV", "None", "Frob", "custom"),
                A=NULL, M=NULL, ...) {
  
  
  chk::chk_true(length(data) > 1)
  for (i in 1:length(data)) {
    chk::chkor(chk::chk_matrix(data[[i]]), chk::chk_s4_class(data[[i]], "Matrix"))
  }
  
  nrs <- sapply(data, nrow)
  chk::chk_true(all(nrs == nrs[1]))
  nr <- nrs[1]
  
  normalization <- match.arg(normalization)
  
  if (normalization == "custom") {
    chk::chkor(chk::chk_not_null(A), chk::chk_not_null(M))
  }
  
  S <- length(data)
  if (is.null(names(data))) {
    names(data) <- paste0("B", 1:S)
  }
  
  # Preprocessing using utility function
  # Set check_consistent_ncol=FALSE as MFA handles potential concatenation later
  preproc_result <- prepare_block_preprocessors(data, preproc, check_consistent_ncol = FALSE)
  proclist <- preproc_result$proclist
  strata <- preproc_result$Xp # Renamed from Xp for consistency with original MFA code
  
  ## calculate block normalization factors
  if (normalization != "custom") {
    alpha <- normalization_factors(strata, type=normalization)
    A <- rep(alpha, sapply(strata, ncol))
  } else {
    alpha <- rep(1, length(strata))
  }
  
  ## compute block indicees
  block_indices <- list()
  ind <- 1
  for (i in 1:length(strata)) {
    block_indices[[i]] <- seq(ind, ind+ncol(strata[[i]])-1)
    ind <- ind + ncol(strata[[i]])
  }
  
  proc <- multivarious::concat_pre_processors(proclist, block_indices)

    ## fit genpca
  if (!requireNamespace("genpca", quietly = TRUE)) {
    stop("Package 'genpca' needed for MFA analysis. Please install it.", call. = FALSE)
  }
  Xp <- do.call(cbind, strata)
  fit <- genpca::genpca(Xp, 
            preproc=multivarious::pass(),
            A=A, 
            M=M,
            ncomp=ncomp,
            ...)
  
  fit[["block_indices"]] <- block_indices
  fit[["alpha"]] <- alpha
  fit[["normalization"]] <- normalization
  fit[["names"]] <- names(data)
  
  ## this is awkward...
  ## instead, we need a "delegation" mechanism, where a multiblock projector simply wraps a projector
  ## here, we rely on the fact that we use "pass()" pre-processing for inner genpca fit
  fit[["preproc"]] <- proc

  # Construct the final multiblock_biprojector
  mfa_result <- multivarious::multiblock_biprojector(
      v = fit$v,              # Loadings from genpca (concatenated space)
      s = fit$s,              # Scores from genpca
      sdev = fit$sdev,        # Singular values from genpca
      preproc = proc,         # Use the concatenated block-aware preprocessor
      block_indices = block_indices,
      # Pass MFA specific info via ...
      alpha = alpha,
      normalization = normalization,
      names = names(data),
      # Pass relevant genpca info via ... as well
      ou = fit$ou,
      ov = fit$ov,
      A_genpca = fit$A, # Renamed to avoid conflict if A was input
      M_genpca = fit$M,
      # Add the class
      classes = "mfa"
  )
  
  
  return(mfa_result)
}
