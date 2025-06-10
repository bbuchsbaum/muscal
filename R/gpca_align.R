#' Generalized PCA Alignment for Multi-Block Data
#'
#' @title Generalized PCA Alignment
#'
#' @description
#' Performs Generalized Principal Component Analysis (GPCA) alignment on multi-block data
#' with a focus on both within-block and between-block relationships. This method extends
#' traditional PCA by incorporating a custom similarity metric and balancing within-group
#' and between-group structure.
#'
#' @details
#' The GPCA alignment method proceeds through several steps:
#'
#' 1. Preprocessing:
#'    * Applies specified preprocessing to each data block
#'    * Concatenates preprocessed blocks
#'
#' 2. Similarity Matrix Construction:
#'    * Computes within-block similarity matrices
#'    * Computes between-block similarity matrix
#'    * Combines matrices with mixing parameter u
#'
#' 3. GPCA Computation:
#'    * Normalizes similarity matrix
#'    * Performs generalized matrix decomposition
#'    * Projects data onto principal components
#'
#' The method is particularly useful for:
#' * Aligning multiple data blocks with shared structure
#' * Balancing within and between-block variation
#' * Supervised dimensionality reduction
#'
#' @param data A hyperdesign object containing multiple data blocks
#' @param y The response variable (unquoted column name)
#' @param preproc Preprocessing function (default: center())
#' @param ncomp Number of components to extract (default: 2)
#' @param simfun Function to compute similarity matrix
#' @param csimfun Optional custom similarity function for cross-block comparisons
#' @param u Mixing parameter between within and between similarity (default: 0.5)
#' @param lambda Ridge regularization parameter (default: 0.1)
#'
#' @return A `gpca_align` object containing:
#' \itemize{
#'   \item \code{v} - Loading vectors
#'   \item \code{s} - Component scores
#'   \item \code{sdev} - Standard deviations of components
#'   \item \code{preproc} - Preprocessing functions used
#'   \item \code{block_indices} - Indices for each data block
#'   \item \code{labels} - Response variable labels
#' }
#'
#' @examples
#' \donttest {
#' # Create example data blocks
#' data1 <- matrix(rnorm(100*10), 100, 10)
#' data2 <- matrix(rnorm(100*15), 100, 15)
#' 
#' # Create labels
#' labels <- factor(rep(1:2, each=50))
#' 
#' # Define similarity function
#' simfun <- function(y) {
#'   outer(y, y, "==") * 1
#' }
#' 
#' # Create hyperdesign
#' design <- data.frame(y=labels)
#' hd <- hyperdesign(list(block1=data1, block2=data2), design)
#' 
#' # Perform GPCA alignment
#' result <- gpca_align(hd, y, simfun=simfun)
#' print(result)
#' }
#'
#'
#' @seealso
#' \code{\link{genpca}} for the underlying GPCA implementation
#'
#' @import PRIMME
#' @importFrom Matrix bdiag Diagonal
#' @importFrom purrr map
#' @importFrom dplyr select pull
#' @import genpca
#' @export
gpca_align <- function(data, y, ...) {
  UseMethod("gpca_align")
}

#' @rdname gpca_align
#' @import PRIMME
#' @importFrom Matrix bdiag Diagonal
#' @importFrom purrr map
#' @importFrom dplyr select pull
#' @import genpca
#' @import multivarious
#' @export
gpca_align.hyperdesign <- function(data, y, 
                    preproc=multivarious::center(), 
                    ncomp=2,
                    simfun,
                    csimfun=NULL,
                    u=.5,
                    lambda=.1,
                    ...) {
  
  # TODO: Add check/warning if csimfun is provided, as it's not used.
  
  y <- rlang::enquo(y)
  
  label_list <- purrr::map(data, function(x) x$design %>% select(!!y) %>% pull(!!y))
  labels <- factor(unlist(label_list))
  label_set <- levels(labels)
  
  #subjects <- purrr::map(data, function(x) x$design %>% select(subject) %>% pull())
  M_within <- Matrix::bdiag(lapply(label_list, function(l) Matrix(simfun(l), sparse=TRUE)))
  M_between <- simfun(labels) - M_within
  M_between <- M_between/length(label_list)
  M <- u*M_within + (1-u)*M_between
  
  ninstances <- length(labels)
  nsets <- length(data)
  
  pdata <- multivarious::init_transform(data, preproc) 
  proclist <- attr(pdata, "preproc")
  
  names(proclist) <- names(pdata)
  
  block_indices <- block_indices(pdata)
  proc <- multivarious::concat_pre_processors(proclist, block_indices)
  names(block_indices) <- names(pdata)
  #browser()
  #M <- simfun(labels)
  M <- M + Matrix::Diagonal(x=rep(lambda, nrow(M)))
  evm <- PRIMME::eigs_sym(M, NEig=1,  which="LA",method='PRIMME_DEFAULT_MIN_MATVECS')
  M <- M/evm$values[1]
  
  X_block <- Matrix::bdiag(lapply(pdata, function(x) x$x))
  ret <- genpca::genpca(X_block, M=M, ncomp=ncomp, preproc=multivarious::pass())

  multiblock_biprojector(
    v=ret$v,
    s=ret$s,
    sdev=ret$sdev,
    preproc=proc,
    block_indices=block_indices,
    labels=labels,
    classes="gpca_align"
  )
}

#' @export
print.gpca_align <- function(x, ...) {
  cat("\n")
  cat(crayon::bold(crayon::blue("✧ GPCA Alignment Results ✧")), "\n\n")
  
  # Model dimensions
  cat(crayon::bold("Model Information:"), "\n")
  cat("  • Number of components:", crayon::green(length(x$sdev)), "\n")
  cat("  • Number of blocks:", crayon::green(length(x$block_indices)), "\n")
  
  # Block information
  cat("\n", crayon::bold("Block Information:"), "\n")
  for (i in seq_along(x$block_indices)) {
    block_name <- names(x$block_indices)[i]
    if (is.null(block_name) || block_name == "") block_name <- paste("Block", i)
    block_size <- length(x$block_indices[[i]]) # Use length of index vector
    cat(sprintf("  • %s: %s variables\n", 
                crayon::blue(block_name), 
                crayon::green(block_size)))
  }
  
  # Explained variance
  var_explained <- (x$sdev^2 / sum(x$sdev^2)) * 100
  cumvar <- cumsum(var_explained)
  
  cat("\n", crayon::bold("Variance Explained:"), "\n")
  for (i in seq_along(var_explained)) {
    cat(sprintf("  • Component %d: %s%% (Cumulative: %s%%)\n", 
                i, 
                crayon::yellow(format(var_explained[i], digits=2)), 
                crayon::yellow(format(cumvar[i], digits=2))))
  }
  
  # Label information
  if (!is.null(x$labels)) {
    n_classes <- length(unique(x$labels))
    cat("\n", crayon::bold("Response Variable used for Alignment:"), "\n")
    cat("  • Number of classes/levels:", crayon::green(n_classes), "\n")
    class_table <- table(x$labels)
    # Show only a few levels if many
    max_levels_show <- 10
    levels_to_show <- names(class_table)
    if (length(levels_to_show) > max_levels_show) {
        levels_to_show <- c(head(levels_to_show, max_levels_show %/% 2), "...", tail(levels_to_show, max_levels_show %/% 2))
    }
    for (level_name in levels_to_show) {
        if (level_name == "...") {
            cat(sprintf("  • %s\n", crayon::blue(level_name)))
        } else {
            cat(sprintf("  • %s: %s samples\n", 
                      crayon::blue(level_name), 
                      crayon::green(class_table[level_name])))
        }
    }
    if (length(class_table) > max_levels_show) {
         cat(sprintf("  • (Total %d levels)\n", length(class_table)))
    }
  }
  
  cat("\n")
  invisible(x)
}

#' @rdname gpca_align
#' @export
gpca_align.default <- function(data, y, ...) {
  stop("`gpca_align` requires a `hyperdesign` object as input. The provided object is of class: ", 
       paste(class(data), collapse=", "), ". \nPlease construct a hyperdesign using `hyperdesign()` or `df_to_hyperdesign()`.",
       call. = FALSE)
}