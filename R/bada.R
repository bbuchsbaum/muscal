#' @import abind
#' @importFrom stats quantile sd
#' @noRd
# A helper function to summarize a matrix
summarize_matrix <- function(matrices, .f) {
  # Stack matrices
  stacked_mat <- abind::abind(matrices, along = 3)
  # Apply function along the third dimension
  apply(stacked_mat, c(1, 2), .f)
}


#' Summarize Bootstrap Resampling
#'
#' This function summarizes the results of bootstrap resampling by computing the mean, standard deviation,
#' and percentiles for each list of matrices in the input tibble.
#'
#' @param boot_ret A tibble containing the results of bootstrap resampling. It should contain two columns:
#'        - boot_i: A list column where each list element is a matrix.
#'        - boot_j: A list column where each list element is a matrix.
#' @param alpha The percentile level for the computation of the lower and upper percentiles (default is 0.05).
#' @return A list containing the mean, standard deviation, upper and lower percentiles for each set of matrices.
#' @importFrom stats quantile sd
#' @noRd
summarize_boot <- function(boot_ret, alpha=.05){
  
  # Apply summarize_matrix for mean, sd, and percentiles
  boot_i_mean <- summarize_matrix(boot_ret$boot_i, mean)
  boot_i_sd <- summarize_matrix(boot_ret$boot_i, sd)
  boot_i_upper <- summarize_matrix(boot_ret$boot_i, function(x) quantile(x, 1-alpha))
  boot_i_lower <- summarize_matrix(boot_ret$boot_i, function(x) quantile(x, alpha))
  
  boot_j_mean <- summarize_matrix(boot_ret$boot_j, mean)
  boot_j_sd <- summarize_matrix(boot_ret$boot_j, sd)
  boot_j_upper <- summarize_matrix(boot_ret$boot_j, function(x) quantile(x, 1-alpha))
  boot_j_lower <- summarize_matrix(boot_ret$boot_j, function(x) quantile(x, alpha))
  
  # Create a list of results
  list(
    boot_scores_mean = boot_i_mean,
    boot_scores_sd = boot_i_sd,
    boot_scores_upper = boot_i_upper,
    boot_scores_lower = boot_i_lower,
    boot_lds_mean = boot_j_mean,
    boot_lds_sd = boot_j_sd,
    boot_lds_upper = boot_j_upper,
    boot_lds_lower = boot_j_lower
  )
}


#' Bootstrap Resampling for bada Multivariate Models
#'
#' Perform bootstrap resampling on a multivariate bada model to estimate the variability 
#' of components and scores.
#'
#' @param x A fitted bada model object that has been fit to a training dataset.
#' @param data The dataset on which the bootstrap resampling should be performed.
#' @param nboot An integer specifying the number of bootstrap resamples to perform (default is 500).
#' @param alpha The percentile level for the computation of the lower and upper percentiles (default is 0.05).
#' @param verbose Logical; if `TRUE`, progress information is printed during
#'   resampling. Defaults to `FALSE`.
#' @param ... Additional arguments to be passed to the specific model implementation of `bootstrap`.
#' @details The function returns a list containing the summarized bootstrap resampled components and scores for the model. 
#' The returned list contains eight elements:
#'    * `boot_scores_mean`: A matrix which is the mean of all bootstrapped scores matrices.
#'    * `boot_scores_sd`: A matrix which is the standard deviation of all bootstrapped scores matrices.
#'    * `boot_scores_upper`: A matrix which is the upper alpha percentile of all bootstrapped scores matrices.
#'    * `boot_scores_lower`: A matrix which is the lower alpha percentile of all bootstrapped scores matrices.
#'    * `boot_lds_mean`: A matrix which is the mean of all bootstrapped loadings matrices.
#'    * `boot_lds_sd`: A matrix which is the standard deviation of all bootstrapped loadings matrices.
#'    * `boot_lds_upper`: A matrix which is the upper alpha percentile of all bootstrapped loadings matrices.
#'    * `boot_lds_lower`: A matrix which is the lower alpha percentile of all bootstrapped loadings matrices.
#' The dimensions of each matrix in the list correspond to the dimensions of the respective matrices in the input data.
#' @importFrom furrr future_map furrr_options
#' @importFrom dplyr tibble bind_rows
#' @importFrom purrr map
#' @importFrom multivarious bootstrap
#' @importFrom stats sd
#' @rdname bada
#' @export
#' @method bootstrap bada
bootstrap.bada <- function(x, data, nboot=500, alpha=.05, verbose = FALSE, ...) {
  sdat <- split(data, x$subjects)
  ## subject-split and preprocessed data
  strata <- seq_along(sdat) %>% purrr::map(function(i) {
    p <- x$proclist[[i]]
    Xi <- sdat[[i]]$x
    Xout <- init_transform(p, Xi)
    multidesign(Xout, sdat[[i]]$design)
  })


  boot_ret <- furrr::future_map(1:nboot, function(i) {
    if (isTRUE(verbose)) message(i)
    boot_indices <- sample(1:length(strata), replace=TRUE)
    boot_strata <- strata[boot_indices]
    ## group designs
    Dc <- boot_strata %>% purrr::map(function(s) {
      summarize_by(s, !!x$y_var)
    })
    
    ## group barycenters
    Xc <-Reduce("+", lapply(Dc, "[[", "x"))/length(Dc)
    
    ## requires consistent ordering. Better to extract from design
    row.names(Xc) <- x$label_set
    ncomp <- min(x$ncomp, nrow(Xc))
    
    boot_i <- Xc %*% x$v
    
    variance <- sdev(x)^2
    #sc <- x$fscores
    #sc <- boot_j
    boot_j <- t(Xc) %*% (x$fscores) %*% diag(1/variance, nrow=length(variance), ncol=length(variance))
    
    dplyr::tibble(i=i, boot_i=list(boot_i), boot_j=list(boot_j))
    
    ## group pca
    ##pca_group <- pca(Xc, ncomp=ncomp, preproc=pass())
  },.options = furrr::furrr_options(seed=TRUE)) %>% bind_rows()
    
  ret <- summarize_boot(boot_ret,alpha)
  class(ret) <- c("bootstrap_bada_result", "list")
  ret
}


#' @import chk
#' @param resdim pca dimensionality for residual analysis (only relevant if `rescomp` > 0)
#' @param rescomp number of final residual components (default = 0, no residual aanalysis)
#' @param ... Additional arguments passed to methods (currently unused).
#' @inheritParams bada
#' @rdname bada
#' @importFrom multidesign summarize_by
#' @importFrom multivarious pca sdev pass init_transform prep discriminant_projector within_class_scatter between_class_scatter scores
#' @importFrom stats prcomp qr qr.Q eigen solve interaction colMeans sd
#' @importFrom rlang enquo quo_get_expr
#' @importFrom dplyr select pull
#' @importFrom chk chk_true
#' @export
bada.multidesign <- function(data, y, subject, preproc=center(), ncomp=2,
                             resdim=20, rescomp=0, ...) {
  y_quo <- rlang::enquo(y)
  subject_quo <- rlang::enquo(subject)
  
  labels <- factor(data$design %>% dplyr::select(!!y_quo) %>% dplyr::pull(!!y_quo))
  label_set <- levels(labels)
  
  subjects <- factor(data$design %>% dplyr::select(!!subject_quo) %>% dplyr::pull(!!subject_quo))
  subject_set <- levels(subjects)
  

  
  ## data split by subject
  sdat <- split(data, subject)
  
  
  ## pre-processors, one per subject
  proclist <- lapply(seq_along(sdat), function(sd) {
    multivarious:::fresh(preproc) %>% prep()
  })
  
  names(proclist) <- as.character(subject_set)
  
  ## subject-split and preprocessed data
  strata <- seq_along(sdat) %>% purrr::map(function(i) {
    p <- multivarious::prep(proclist[[i]], sdat[[i]]$x) 
    Xi <- sdat[[i]]$x
    Xout <- multivarious::init_transform(p, Xi)
    multidesign(Xout, sdat[[i]]$design)
  })
  
  block_indices <- list()
  ind <- 1
  for (i in 1:length(strata)) {
    block_indices[[i]] <- seq(ind, ind+ncol(strata[[i]]$x)-1)
    ind <- ind + ncol(strata[[i]]$x)
  }
  
  names(block_indices) <- as.character(subject_set)
  
  
  ## group designs
  Dc <- strata %>% purrr::map(function(s) {
    summarize_by(s, !!y_quo)
  })
  
  ## group barycenters
  Xc <-Reduce("+", lapply(Dc, "[[", "x"))/length(Dc)
  
  ## requires consistent ordering. Better to extract from design
  row.names(Xc) <- label_set
  ncomp <- min(ncomp, nrow(Xc))

  ## group pca
  pca_group <- multivarious::pca(Xc, ncomp=ncomp, preproc=multivarious::pass())
  
  ## residual analysis
  if (rescomp > 0) {
    chk::chk_true(resdim > 0)
    residual_strata <- strata %>% purrr::map(function(s) {
      levs <- s$design %>% dplyr::pull(!!y_quo)
      s$x <- s$x - Xc[levs,,drop=FALSE]
      s
    })
  
    Xresid <- do.call(rbind, residual_strata %>% purrr::map( ~ .x$x))
  
    pca_resid <- pca(Xresid, ncomp=resdim, method="irlba")
    Xpca_resid <- multivarious::scores(pca_resid)
  
    Sw <- multivarious::within_class_scatter(Xpca_resid, interaction(subjects, labels))
    Sb <- multivarious::between_class_scatter(Xpca_resid, interaction(subjects, labels), 
                              colMeans(Xpca_resid))
  
    eigout <- eigen(solve(Sw, Sb))
    #scores <- Xpca_resid %*% eigout$vectors
    resid_v <- pca_resid$v %*% eigout$vectors
    v <- cbind(pca_group$v, resid_v)
    vq <- qr(v)
    v <- qr.Q(vq)
  } else {
    resdim <- 0
    rescomp <- 0
    v <- pca_group$v
  }
  
  ## compute projections, one subject at a time.
  s <- do.call(rbind, strata %>% purrr::map(function(s) {
    s$x %*% v
  }))
  
  proc <- multivarious:::concat_pre_processors(proclist, block_indices)
  
  multivarious:::discriminant_projector(v=v, 
                                        s=s, 
                                        fscores=pca_group$s,
                                        sdev=apply(s, 2, stats::sd),
                                        preproc = proc,
                                        proclist = proclist,
                                        labels=labels, 
                                        label_set=label_set,
                                        resdim=resdim,
                                        rescomp=rescomp,
                                        subjects=subjects,
                                        barycenters=Xc,
                                        block_indices=block_indices,
                                        subject_var=subject,
                                        y_var=y,
                                        classes="bada")
 
}


#' Reprocess New Data for a Barycentric Discriminant Analysis (BaDA) Model
#'
#' This function transforms new data using the preprocessing pipeline from a fitted BaDA model.
#' It handles different preprocessing scenarios based on the provided parameters.
#'
#' @param x A fitted BaDA model object.
#' @param new_data A numeric matrix of new data to be processed.
#' @param colind An optional integer vector specifying column indices to use within blocks.
#'   If NULL and block is also NULL, all blocks are used. If NULL but block is provided,
#'   all columns in the specified block are used.
#' @param block An optional character string specifying which block's preprocessing to apply.
#'   If NULL and colind is also NULL, preprocessing is averaged across all blocks.
#'   If provided, only the specified block's preprocessing is applied.
#'
#' @details
#' The function handles three scenarios:
#' 1. When both colind and block are NULL: The function applies preprocessing for each block
#'    and averages the results.
#' 2. When block is provided: The function applies preprocessing specific to the named block.
#'    If colind is also provided, it acts as a relative subset within the block.
#' 3. When only colind is provided: The function applies preprocessing for each block
#'    using the specified column indices and averages the results.
#'
#' @return A preprocessed numeric matrix with the same number of rows as the input data.
#'
#' @rdname bada
#' @export
reprocess.bada <- function(x, new_data, colind=NULL, block=NULL) {
  if (is.null(colind) && is.null(block)) {
    ## how to pre-process when you don't know the subject?
    ## we pre-process every way and average.
    chk::chk_equal(ncol(new_data), shape(x)[1])
    
    ## pre-process every which way...
    Reduce("+", lapply(seq_along(x$block_indices), function(i) {
      multivarious::apply_transform(x$preproc, new_data, colind=x$block_indices[[i]])
    }))/length(x$block_indices)
    
  } else if (!is.null(block)) {
    ## pre-process along one block
    chk::chk_character(block)
    chk::chk_equal(ncol(new_data), shape(x)[1])
    sind <- x$block_indices[[block]]
    if (!is.null(colind)) {
      ## relative subset using colind
      sind <- sind[colind]
    }
    multivarious::apply_transform(x$preproc, new_data, colind=sind)
  } else {
    ## colind not null. pre-process every which way using colind per block
    Reduce("+", lapply(seq_along(x$block_indices), function(i) {
      multivarious::apply_transform(x$preproc, new_data, colind=colind)
    }))/length(x$block_indices)
  }
  
}

#' Project New Data onto a Barycentric Discriminant Analysis (BaDA) Model
#'
#' This function projects new data onto a previously fitted BaDA model,
#' returning the scores of the new data in the space of the original model.
#'
#' @param x A fitted BaDA model object.
#' @param new_data A numeric matrix of new data to be projected.
#' @param block An optional character string specifying which block's preprocessing
#'   to apply before projection. If missing, the generic method is used.
#'
#' @details
#' When a specific block is provided, the function first reprocesses the data using
#' that block's preprocessing pipeline, then projects it onto the model space.
#' If no block is specified, it delegates to the default project method.
#'
#' @return A numeric matrix of projected scores.
#'
#' @rdname bada
#' @export
project.bada <- function(x, new_data, block) {
  if (missing(block)) {
    NextMethod(x,new_data)
  } else {
    #Xp <- multivarious::apply_transform(preproc, new_data)
    reprocess(x, new_data, block=block) %*% multivarious::coef.projector(x)
  }
}

