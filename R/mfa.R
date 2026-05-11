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
  
  first_sv <- function(mat, method_pref = "irlba") {
    method <- .muscal_svd_method(mat, ncomp = 1)
    if (!identical(method_pref, "irlba") && !identical(method, "base")) {
      method <- method_pref
    }
    tryCatch(
      multivarious::svd_wrapper(mat, ncomp = 1, method = method)$sdev[1],
      error = function(e) multivarious::svd_wrapper(mat, ncomp = 1, method = "base")$sdev[1]
    )
  }
  
  alpha <- if (type == "MFA") {
    unlist(lapply(blocks, function(X) 1 / (first_sv(X)^2)))
  } else if (type == "RV") {
    smat <- compute_sim_mat(blocks, function(x1,x2) MatrixCorrelation::RV(x1,x2))
    diag(smat) <- 1
    u1 <- multivarious::svd_wrapper(smat, ncomp=1, method = .muscal_svd_method(smat, ncomp = 1))$u[,1]
    abs(u1) / sum(abs(u1))
  } else if (type == "RV2") {
    smat <- compute_sim_mat(blocks, function(x1,x2) MatrixCorrelation::RV2(x1,x2))
    diag(smat) <- 1
    u1 <- multivarious::svd_wrapper(smat, ncomp=1, method = .muscal_svd_method(smat, ncomp = 1))$u[,1]
    abs(u1) / sum(abs(u1))
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
#' \donttest{
#' # Apply MFA to a list of matrices
#' X <- replicate(5, { matrix(rnorm(10*10), 10, 10) }, simplify=FALSE)
#' res <- mfa(X, ncomp=3, normalization="MFA")
#' }
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
#' @param missing Missing-data handling mode. `"error"` preserves the historical
#'   complete-data contract. `"regularized"` and `"em"` run iterative MFA
#'   imputation before the final complete-data MFA fit. `"nipals"` is reserved
#'   for a future direct available-data backend.
#' @param ncp_impute Integer; number of MFA components used inside the iterative
#'   imputation loop.
#' @param missing_tol Numeric convergence tolerance for iterative imputation.
#' @param missing_maxiter Integer maximum number of imputation iterations.
#' @param return_imputed Logical; if `TRUE`, attach completed blocks as
#'   `imputed_data`.
#'
#' @examples
#' \donttest{
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
#' pred <- predict(cfier, do.call(cbind, X)[1:2,])
#'
#' # Create a classifier using a specific block
#' cfier2 <- multivarious::classifier(res, new_data=X[[2]], labels=labs,
#'                                   colind=res$block_indices[[2]])
#' pred2 <- predict(cfier2, X[[2]][1:2,])
#' }
#' @export
mfa.multiblock <- function(data, preproc=center(), ncomp=2,
                normalization=c("MFA", "RV", "None", "Frob", "custom"),
                A=NULL, M=NULL,
                missing = c("error", "regularized", "em", "nipals"),
                ncp_impute = ncomp,
                missing_tol = 1e-6,
                missing_maxiter = 100,
                return_imputed = FALSE,
                ...) {
  fit_call <- match.call(expand.dots = FALSE)
  fit_dots <- list(...)
  missing <- match.arg(missing)

  has_na <- any(vapply(data, anyNA, logical(1L)))
  if (!has_na) {
    return(.mfa_complete(
      data = data,
      preproc = preproc,
      ncomp = ncomp,
      normalization = normalization,
      A = A,
      M = M,
      fit_call = fit_call,
      fit_dots = fit_dots,
      ...
    ))
  }

  if (identical(missing, "error")) {
    stop(
      "mfa() does not accept NA values unless `missing` is set to ",
      "'regularized', 'em', or 'nipals'.",
      call. = FALSE
    )
  }

  if (identical(missing, "nipals")) {
    stop(
      "missing = 'nipals' is not implemented yet; use missing = 'regularized' or 'em'.",
      call. = FALSE
    )
  }

  missing_controls <- .mfa_validate_missing_controls(
    ncp_impute = ncp_impute,
    missing_tol = missing_tol,
    missing_maxiter = missing_maxiter,
    return_imputed = return_imputed
  )
  ncp_impute <- missing_controls$ncp_impute
  missing_tol <- missing_controls$missing_tol
  missing_maxiter <- missing_controls$missing_maxiter
  return_imputed <- missing_controls$return_imputed

  imp <- .mfa_impute_iterative(
    data = data,
    preproc = preproc,
    ncp = ncp_impute,
    normalization = normalization,
    A = A,
    M = M,
    method = missing,
    tol = missing_tol,
    maxiter = missing_maxiter,
    ...
  )

  fit <- .mfa_complete(
    data = imp$data,
    preproc = preproc,
    ncomp = ncomp,
    normalization = normalization,
    A = A,
    M = M,
    fit_call = fit_call,
    fit_dots = fit_dots,
    ...
  )
  fit$missing <- list(
    method = missing,
    ncp_impute = ncp_impute,
    iterations = imp$iterations,
    converged = imp$converged,
    delta = imp$delta,
    mask = imp$mask
  )
  fit$fit_spec$refit <- .mfa_missing_refit_spec(
    data = data,
    preproc = preproc,
    ncomp = ncomp,
    normalization = normalization,
    A = A,
    M = M,
    missing = missing,
    ncp_impute = ncp_impute,
    missing_tol = missing_tol,
    missing_maxiter = missing_maxiter,
    return_imputed = return_imputed,
    fit_dots = fit_dots
  )
  fit$fit_spec$refit_supported <- TRUE
  if (isTRUE(return_imputed)) {
    fit$imputed_data <- imp$data
  }
  fit
}

.mfa_complete <- function(data, preproc=center(), ncomp=2,
                normalization=c("MFA", "RV", "None", "Frob", "custom"),
                A=NULL, M=NULL,
                build_diagnostics = TRUE,
                fit_call = NULL,
                fit_dots = list(),
                ...) {
  if (is.null(fit_call)) {
    fit_call <- match.call(expand.dots = FALSE)
  }
  if (is.null(fit_dots)) {
    fit_dots <- list(...)
  }
  A_input <- A
  M_input <- M
  
  
  chk::chk_true(length(data) > 1)
  for (i in 1:length(data)) {
    chk::chkor_vld(chk::vld_matrix(data[[i]]), chk::vld_s4_class(data[[i]], "Matrix"))
  }
  
  nrs <- sapply(data, nrow)
  chk::chk_true(all(nrs == nrs[1]))
  nr <- nrs[1]
  
  normalization <- match.arg(normalization)
  
  if (normalization == "custom" && is.null(A) && is.null(M)) {
    stop('At least one of A or M must be provided when normalization = "custom"')
  }
  
  S <- length(data)
  if (is.null(names(data))) {
    names(data) <- paste0("B", 1:S)
  }
  data_refit <- lapply(data, function(x) as.matrix(x))
  
  # Preprocessing using utility function
  # Set check_consistent_ncol=FALSE as MFA handles potential concatenation later
  preproc_result <- prepare_block_preprocessors(data, preproc, check_consistent_ncol = FALSE)
  proclist <- preproc_result$proclist
  fitted_proclist <- .muscal_materialize_block_preprocessors(data, proclist)
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
  
  proc <- multivarious::concat_pre_processors(fitted_proclist, block_indices)

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

  # -------------------------------------------------------------------------
  # Derived quantities for plotting / diagnostics
  # -------------------------------------------------------------------------
  if (isTRUE(build_diagnostics)) {
    partial_scores <- lapply(seq_along(strata), function(i) {
      idx <- block_indices[[i]]
      strata[[i]] %*% fit$v[idx, , drop = FALSE]
    })
    names(partial_scores) <- names(data)

    cor_loadings <- tryCatch(
      {
        C <- stats::cor(Xp, fit$s)
        C[!is.finite(C)] <- 0
        C
      },
      error = function(e) NULL
    )

    rv <- tryCatch(
      {
        M <- compute_sim_mat(strata, function(x1, x2) MatrixCorrelation::RV(x1, x2))
        diag(M) <- 1
        rownames(M) <- colnames(M) <- names(data)
        M
      },
      error = function(e) NULL
    )

    rv2 <- tryCatch(
      {
        M <- compute_sim_mat(strata, function(x1, x2) MatrixCorrelation::RV2(x1, x2))
        diag(M) <- 1
        rownames(M) <- colnames(M) <- names(data)
        M
      },
      error = function(e) NULL
    )
  } else {
    partial_scores <- NULL
    cor_loadings <- NULL
    rv <- NULL
    rv2 <- NULL
  }

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
      partial_scores = partial_scores,
      cor_loadings = cor_loadings,
      rv = rv,
      rv2 = rv2,
      # Pass relevant genpca info via ... as well
      ou = fit$ou,
      ov = fit$ov,
      A_genpca = fit$A, # Renamed to avoid conflict if A was input
      M_genpca = fit$M,
      # Add the class
      classes = "mfa"
  )
  
  mfa_result <- .muscal_attach_fit_contract(
    mfa_result,
    method = "mfa",
    task = "reconstruction",
    oos_types = c("scores", "reconstruction"),
    fit_call = fit_call,
    refit_supported = TRUE,
    prediction_target = "blocks",
    refit = .muscal_make_refit_spec(
      data = data_refit,
      fit_fn = function(data) {
        do.call(
          mfa,
          c(
            list(
              data = data,
              preproc = preproc,
              ncomp = ncomp,
              normalization = normalization,
              A = A_input,
              M = M_input
            ),
            fit_dots
          )
        )
      },
      bootstrap_fn = function(data) {
        n <- nrow(data[[1L]])
        idx <- sample.int(n, size = n, replace = TRUE)
        lapply(data, function(block) block[idx, , drop = FALSE])
      },
      permutation_fn = function(data) {
        lapply(data, function(block) block[sample.int(nrow(block)), , drop = FALSE])
      },
      resample_unit = "rows"
    )
  )

  return(mfa_result)
}

.mfa_impute_iterative <- function(data, preproc, ncp, normalization, A, M,
                                  method = c("regularized", "em"),
                                  tol = 1e-6, maxiter = 100, ...) {
  method <- match.arg(method)
  X <- lapply(data, function(x) as.matrix(x))
  if (is.null(names(X))) {
    names(X) <- paste0("B", seq_along(X))
  }
  mask <- lapply(X, is.na)

  for (b in seq_along(X)) {
    bad_cols <- which(colSums(!mask[[b]]) == 0L)
    if (length(bad_cols) > 0L) {
      stop(
        sprintf(
          "Cannot impute variable/column %d with all values missing in block '%s'.",
          bad_cols[[1L]], names(X)[[b]]
        ),
        call. = FALSE
      )
    }
  }

  Ximp <- lapply(X, function(block) {
    out <- block
    cm <- colMeans(out, na.rm = TRUE)
    miss <- is.na(out)
    if (any(miss)) {
      out[miss] <- cm[col(out)[miss]]
    }
    out
  })

  converged <- FALSE
  delta <- Inf
  iter <- 0L
  for (iter in seq_len(maxiter)) {
    fit <- .mfa_complete(
      data = Ximp,
      preproc = preproc,
      ncomp = ncp,
      normalization = normalization,
      A = A,
      M = M,
      build_diagnostics = FALSE,
      ...
    )
    Xhat <- .mfa_reconstruct_blocks(fit, Ximp, ncomp = ncp, method = method)

    old_missing <- unlist(Map(function(x, m) x[m], Ximp, mask), use.names = FALSE)
    for (b in seq_along(Ximp)) {
      Ximp[[b]][mask[[b]]] <- Xhat[[b]][mask[[b]]]
    }
    new_missing <- unlist(Map(function(x, m) x[m], Ximp, mask), use.names = FALSE)
    denom <- max(1, sqrt(sum(old_missing^2, na.rm = TRUE)))
    delta <- sqrt(sum((new_missing - old_missing)^2, na.rm = TRUE)) / denom
    if (is.finite(delta) && delta < tol) {
      converged <- TRUE
      break
    }
  }

  list(
    data = Ximp,
    mask = mask,
    iterations = iter,
    converged = converged,
    delta = delta
  )
}

.mfa_validate_missing_controls <- function(ncp_impute, missing_tol,
                                           missing_maxiter, return_imputed) {
  if (!is.numeric(ncp_impute) ||
      length(ncp_impute) != 1L ||
      is.na(ncp_impute) ||
      !is.finite(ncp_impute) ||
      ncp_impute < 1 ||
      abs(ncp_impute - round(ncp_impute)) > 1e-8) {
    stop("`ncp_impute` must be a single integer-like value >= 1.", call. = FALSE)
  }

  if (!is.numeric(missing_tol) ||
      length(missing_tol) != 1L ||
      is.na(missing_tol) ||
      !is.finite(missing_tol) ||
      missing_tol <= 0) {
    stop("`missing_tol` must be a single finite value > 0.", call. = FALSE)
  }

  if (!is.numeric(missing_maxiter) ||
      length(missing_maxiter) != 1L ||
      is.na(missing_maxiter) ||
      !is.finite(missing_maxiter) ||
      missing_maxiter < 1 ||
      abs(missing_maxiter - round(missing_maxiter)) > 1e-8) {
    stop("`missing_maxiter` must be a single integer-like value >= 1.", call. = FALSE)
  }

  if (!is.logical(return_imputed) || length(return_imputed) != 1L || is.na(return_imputed)) {
    stop("`return_imputed` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  list(
    ncp_impute = as.integer(round(ncp_impute)),
    missing_tol = as.numeric(missing_tol),
    missing_maxiter = as.integer(round(missing_maxiter)),
    return_imputed = isTRUE(return_imputed)
  )
}

.mfa_reconstruct_blocks <- function(fit, data, ncomp, method = c("regularized", "em")) {
  method <- match.arg(method)
  comp <- seq_len(min(as.integer(ncomp), ncol(fit$v), ncol(fit$s)))
  Xconcat <- do.call(cbind, data)

  Xhat <- if (identical(method, "em")) {
    multivarious::reconstruct_new(fit, Xconcat, comp = comp)
  } else {
    sdev <- fit$sdev[comp]
    lambda <- sdev^2
    sigma2 <- if (length(lambda) > 1L) min(lambda, na.rm = TRUE) else 0
    shrink <- if (is.finite(sigma2) && sigma2 > 0) {
      lambda / (lambda + sigma2)
    } else {
      rep(1, length(comp))
    }
    scores_new <- multivarious::project(fit, Xconcat)[, comp, drop = FALSE]
    Xhat_p <- scores_new %*%
      diag(shrink, nrow = length(comp), ncol = length(comp)) %*%
      t(fit$v[, comp, drop = FALSE])
    multivarious::inverse_transform(fit$preproc, Xhat_p)
  }

  lapply(fit$block_indices, function(idx) Xhat[, idx, drop = FALSE])
}

.mfa_missing_refit_spec <- function(data, preproc, ncomp, normalization, A, M,
                                    missing, ncp_impute, missing_tol,
                                    missing_maxiter, return_imputed, fit_dots) {
  data_refit <- lapply(data, function(x) as.matrix(x))
  .muscal_make_refit_spec(
    data = data_refit,
    fit_fn = function(data) {
      do.call(
        mfa,
        c(
          list(
            data = data,
            preproc = preproc,
            ncomp = ncomp,
            normalization = normalization,
            A = A,
            M = M,
            missing = missing,
            ncp_impute = ncp_impute,
            missing_tol = missing_tol,
            missing_maxiter = missing_maxiter,
            return_imputed = return_imputed
          ),
          fit_dots
        )
      )
    },
    bootstrap_fn = function(data) {
      n <- nrow(data[[1L]])
      idx <- sample.int(n, size = n, replace = TRUE)
      lapply(data, function(block) block[idx, , drop = FALSE])
    },
    permutation_fn = function(data) {
      lapply(data, function(block) block[sample.int(nrow(block)), , drop = FALSE])
    },
    resample_unit = "rows"
  )
}
