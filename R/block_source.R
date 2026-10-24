#' Out-of-Core Block Sources
#'
#' A block source describes a list of data blocks without holding them in
#' memory. Fitting functions that accept a block source (currently
#' [aligned_mfa()]) read it one column chunk of one block at a time, so peak
#' memory is one chunk plus the method's small summary statistics.
#'
#' `block_source()` wraps a loader function, which is the general mechanism:
#' any storage (HDF5, file-backed matrices, NIfTI volumes, a database) can be
#' exposed by writing `loader(k, cols)`. `as_block_source()` converts an
#' in-memory list of matrices, which is mainly useful for testing the
#' out-of-core code paths against the in-memory ones.
#'
#' @param loader Function `function(k, cols)` returning block `k` restricted to
#'   the integer column indices `cols`, as a numeric matrix with `nrow[k]` rows
#'   and `length(cols)` columns. Every call for block `k` must return the rows
#'   in the same order. Missing and infinite values are not supported.
#' @param nrow,ncol Integer vectors giving the number of rows and columns of
#'   each block.
#' @param names Optional block names; defaults to `names(nrow)` or `X1, X2, ...`.
#' @param chunk_size Maximum number of columns read per call to `loader`. Peak
#'   memory per read is about `8 * nrow[k] * chunk_size` bytes.
#' @param x An object to convert: a list of matrices/data frames, or an
#'   existing block source.
#' @param ... Passed to `block_source()` (e.g. `chunk_size`).
#'
#' @return An object of class `muscal_block_source`.
#'
#' @examples
#' X <- list(A = matrix(rnorm(40), 10), B = matrix(rnorm(60), 10))
#' src <- block_source(
#'   loader = function(k, cols) X[[k]][, cols, drop = FALSE],
#'   nrow = c(10, 10), ncol = c(4, 6), names = c("A", "B"), chunk_size = 3
#' )
#' src
#' @export
block_source <- function(loader, nrow, ncol, names = NULL, chunk_size = 5000L) {
  chk::chk_function(loader)
  nrow <- as.integer(nrow)
  ncol <- as.integer(ncol)
  if (length(nrow) < 1L || length(nrow) != length(ncol)) {
    stop("`nrow` and `ncol` must be non-empty and of equal length.", call. = FALSE)
  }
  if (anyNA(nrow) || anyNA(ncol) || any(nrow < 1L) || any(ncol < 1L)) {
    stop("`nrow` and `ncol` must be positive integers.", call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)
  chk::chk_scalar(chunk_size)
  chk::chk_gte(chunk_size, 1L)

  if (is.null(names)) {
    names <- names(nrow) %||% paste0("X", seq_along(nrow))
  }
  names <- as.character(names)
  if (length(names) != length(nrow) || anyDuplicated(names)) {
    stop("`names` must be unique and have one entry per block.", call. = FALSE)
  }

  structure(
    list(
      loader = loader,
      nrow = stats::setNames(nrow, names),
      ncol = stats::setNames(ncol, names),
      names = names,
      chunk_size = chunk_size
    ),
    class = "muscal_block_source"
  )
}

#' @rdname block_source
#' @export
as_block_source <- function(x, ...) {
  UseMethod("as_block_source")
}

#' @rdname block_source
#' @export
as_block_source.muscal_block_source <- function(x, ...) {
  x
}

#' @rdname block_source
#' @export
as_block_source.list <- function(x, ...) {
  x <- lapply(x, as.matrix)
  nms <- names(x) %||% paste0("X", seq_along(x))
  block_source(
    loader = function(k, cols) x[[k]][, cols, drop = FALSE],
    nrow = vapply(x, nrow, integer(1L)),
    ncol = vapply(x, ncol, integer(1L)),
    names = nms,
    ...
  )
}

#' @rdname block_source
#' @export
as_block_source.default <- function(x, ...) {
  stop("Cannot convert an object of class '", class(x)[[1L]], "' to a block source.",
       call. = FALSE)
}

#' @export
print.muscal_block_source <- function(x, ...) {
  cat(sprintf("<muscal_block_source> %d blocks, chunk_size = %d\n",
              length(x$names), x$chunk_size))
  dims <- data.frame(block = x$names, nrow = unname(x$nrow), ncol = unname(x$ncol))
  print(dims, row.names = FALSE)
  invisible(x)
}

#' @rdname block_source
#' @param k Block index or name.
#' @param cols Integer column indices to read; defaults to all columns.
#' @return For `read_block()`, the requested numeric matrix.
#' @export
read_block <- function(x, k, cols = NULL) {
  k <- .muscal_source_index(x, k)
  if (is.null(cols)) cols <- seq_len(x$ncol[[k]])
  .muscal_source_read(x, k, cols)
}

#' @export
names.muscal_block_source <- function(x) {
  x$names
}

#' @export
length.muscal_block_source <- function(x) {
  length(unclass(x)$names)
}

# ---------------------------------------------------------------------------
# Internal one-pass operations shared by out-of-core fitting paths
# ---------------------------------------------------------------------------

.muscal_source_index <- function(src, k) {
  if (is.character(k)) {
    i <- match(k, src$names)
    if (is.na(i)) stop("Unknown block '", k, "'.", call. = FALSE)
    return(i)
  }
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 1L || k > length(src$names)) {
    stop("Block index out of range.", call. = FALSE)
  }
  k
}

.muscal_source_chunks <- function(src, k) {
  p <- src$ncol[[k]]
  starts <- seq(1L, p, by = src$chunk_size)
  lapply(starts, function(s) s:min(p, s + src$chunk_size - 1L))
}

.muscal_source_read <- function(src, k, cols) {
  X <- src$loader(k, cols)
  if (is.data.frame(X)) X <- as.matrix(X)
  if (inherits(X, "Matrix")) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X)) {
    stop(sprintf("Loader for block '%s' must return a numeric matrix.", src$names[[k]]),
         call. = FALSE)
  }
  if (nrow(X) != src$nrow[[k]] || ncol(X) != length(cols)) {
    stop(
      sprintf(
        "Loader for block '%s' returned a %d x %d matrix; expected %d x %d.",
        src$names[[k]], nrow(X), ncol(X), src$nrow[[k]], length(cols)
      ),
      call. = FALSE
    )
  }
  if (!all(is.finite(X))) {
    stop(
      sprintf("Block '%s' contains missing or infinite values, which block sources do not support.",
              src$names[[k]]),
      call. = FALSE
    )
  }
  X
}

# Identify a column-local preprocessor ("none", "center", "standardize") by
# probing it, since only those can be fitted chunk by chunk.
.muscal_stream_preproc_kinds <- function(preproc, n_blocks) {
  defs <- if (is.null(preproc) ||
              inherits(preproc, "prepper") ||
              inherits(preproc, "pre_processor") ||
              !is.list(preproc)) {
    rep(list(preproc), n_blocks)
  } else {
    if (length(preproc) != n_blocks) {
      stop(sprintf("If 'preproc' is a list, its length (%d) must match the number of blocks (%d).",
                   length(preproc), n_blocks), call. = FALSE)
    }
    preproc
  }
  probe <- matrix(c(1, 4, 2, 9, 7, 3, 5, 8, 2, 6, 1, 4), nrow = 4, ncol = 3)
  centered <- sweep(probe, 2L, colMeans(probe))
  standardized <- sweep(centered, 2L, apply(probe, 2L, stats::sd), `/`)
  vapply(defs, function(def) {
    if (is.null(def)) return("none")
    out <- tryCatch(
      {
        fitted <- multivarious::fit(multivarious::fresh(def), probe)
        as.matrix(multivarious::transform(fitted, probe))
      },
      error = function(e) NULL
    )
    same <- function(a) !is.null(out) && all(dim(out) == dim(a)) && max(abs(out - a)) < 1e-10
    if (same(probe)) return("none")
    if (same(centered)) return("center")
    if (same(standardized)) return("standardize")
    stop(
      "Block sources support only column-local preprocessing: ",
      "multivarious::pass(), center(), or standardize() (or NULL).",
      call. = FALSE
    )
  }, character(1L))
}

.muscal_stream_apply_prep <- function(X, kind, center, scale) {
  if (identical(kind, "none")) return(X)
  X <- sweep(X, 2L, center)
  if (identical(kind, "standardize")) X <- sweep(X, 2L, scale, `/`)
  X
}

# Visit every column chunk of every block once. `callback(k, cols, Xp)` receives
# the preprocessed chunk; `prep` holds per-block `center`/`scale` vectors.
.muscal_stream_pass <- function(src, kinds, prep, callback, blocks = seq_along(src$names)) {
  for (k in blocks) {
    for (cols in .muscal_source_chunks(src, k)) {
      X <- .muscal_source_read(src, k, cols)
      Xp <- .muscal_stream_apply_prep(X, kinds[[k]], prep[[k]]$center[cols], prep[[k]]$scale[cols])
      callback(k, cols, Xp)
    }
  }
  invisible(NULL)
}

# One pass: column-local preprocessing statistics and the row Gram matrix
# G_k = Xp_k Xp_k' of every block. Each chunk holds all rows of its columns, so
# its means and SDs are exact.
.muscal_stream_gram <- function(src, kinds) {
  n_blocks <- length(src$names)
  G <- lapply(seq_len(n_blocks), function(k) matrix(0, src$nrow[[k]], src$nrow[[k]]))
  prep <- lapply(seq_len(n_blocks), function(k) {
    list(center = numeric(src$ncol[[k]]), scale = rep(1, src$ncol[[k]]),
         col_ss = numeric(src$ncol[[k]]))
  })
  n_bad_sd <- 0L

  for (k in seq_len(n_blocks)) {
    for (cols in .muscal_source_chunks(src, k)) {
      X <- .muscal_source_read(src, k, cols)
      if (!identical(kinds[[k]], "none")) {
        cm <- colMeans(X)
        prep[[k]]$center[cols] <- cm
        if (identical(kinds[[k]], "standardize")) {
          sds <- apply(X, 2L, stats::sd)
          bad <- !is.finite(sds) | sds < .Machine$double.eps
          n_bad_sd <- n_bad_sd + sum(bad)
          sds[bad] <- 1
          prep[[k]]$scale[cols] <- sds
        }
      }
      Xp <- .muscal_stream_apply_prep(X, kinds[[k]], prep[[k]]$center[cols], prep[[k]]$scale[cols])
      prep[[k]]$col_ss[cols] <- colSums(sweep(Xp, 2L, colMeans(Xp))^2)
      G[[k]] <- G[[k]] + tcrossprod(Xp)
    }
    G[[k]] <- (G[[k]] + t(G[[k]])) / 2
  }
  if (n_bad_sd > 0L) {
    warning(sprintf("%d column(s) have zero or non-finite standard deviation; scale set to 1.",
                    n_bad_sd), call. = FALSE)
  }
  names(G) <- names(prep) <- src$names
  list(G = G, prep = prep)
}

# One pass: Xp_k' M_k for a list of n_k x m matrices M_k (p_k x m results).
.muscal_stream_crossprod <- function(src, kinds, prep, M_list) {
  out <- lapply(seq_along(src$names), function(k) {
    matrix(0, src$ncol[[k]], ncol(M_list[[k]]))
  })
  .muscal_stream_pass(src, kinds, prep, function(k, cols, Xp) {
    out[[k]][cols, ] <<- crossprod(Xp, M_list[[k]])
  })
  names(out) <- src$names
  out
}

# Square root R of a PSD Gram matrix (R R' = G) with numerically null
# directions dropped; R keeps the row order of the block.
.muscal_gram_root <- function(G, tol = 1e-10) {
  e <- eigen(G, symmetric = TRUE)
  top <- max(e$values[[1L]], 0)
  keep <- e$values > tol * top
  if (!any(keep)) {
    return(matrix(0, nrow(G), 1L))
  }
  e$vectors[, keep, drop = FALSE] %*%
    diag(sqrt(e$values[keep]), nrow = sum(keep), ncol = sum(keep))
}

.muscal_stream_fitted_preproc <- function(kind, prep, p) {
  dummy <- matrix(0, nrow = 1L, ncol = p)
  def <- switch(
    kind,
    none = multivarious::pass(),
    center = multivarious::center(cmeans = prep$center),
    standardize = multivarious::standardize(cmeans = prep$center, sds = prep$scale)
  )
  multivarious::fit(def, dummy)
}
