.dmfa_coordinate_matrix <- function(x, what) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop(sprintf("`%s` must be a numeric coordinate matrix or data frame.", what), call. = FALSE)
  }
  out <- as.matrix(x)
  if (!is.numeric(out) || nrow(out) < 1L || ncol(out) < 1L || any(!is.finite(out))) {
    stop(sprintf("`%s` must be a nonempty finite numeric coordinate matrix.", what), call. = FALSE)
  }
  out
}

.dmfa_coordinate_names <- function(x, supplied, prefix, what) {
  out <- if (is.null(supplied)) rownames(x) else as.character(supplied)
  if (is.null(out)) out <- paste0(prefix, seq_len(nrow(x)))
  if (length(out) != nrow(x) || anyNA(out) || any(out == "") || anyDuplicated(out)) {
    stop(sprintf("`%s` names must be unique, non-missing, and non-empty.", what), call. = FALSE)
  }
  out
}

#' Build a Sparse Coordinate-to-Feature Map
#'
#' @description
#' Constructs a sparse operator that maps a field on named reference
#' coordinates into a subject's observed coordinate space. Each subject row is
#' connected only to its `k` nearest reference coordinates, found with a
#' k-d tree rather than an all-pairs distance matrix.
#'
#' @details
#' The returned matrix has subject features in rows and reference features in
#' columns, matching the `feature_map` contract of [design_mfa()]. Gaussian and
#' inverse-distance kernels are row-normalized by default. The construction is
#' invariant to applying the same translation or orthogonal transformation to
#' both coordinate sets.
#'
#' Coordinates may be MNI voxels, surface coordinates, or another common
#' Euclidean system. This function does not estimate anatomical registration;
#' any transform into a common coordinate system must be applied beforehand.
#'
#' @param reference_coords Numeric reference-feature-by-coordinate matrix.
#' @param subject_coords Numeric subject-feature-by-coordinate matrix.
#' @param k Number of nearest reference coordinates retained per subject
#'   feature. Values larger than the number of reference features are truncated.
#' @param bandwidth Positive Gaussian bandwidth. When `NULL`, the median
#'   positive retained-neighbor distance is used.
#' @param max_distance Maximum permitted coordinate distance. Subject features
#'   with no retained neighbor inside this radius are rejected.
#' @param kernel `"gaussian"`, `"inverse_distance"`, or `"nearest"`.
#' @param normalize_rows Normalize each subject row to sum to one.
#' @param reference_names Optional reference feature names; otherwise row names
#'   or deterministic generated names are used.
#' @param subject_names Optional subject feature names; otherwise row names or
#'   deterministic generated names are used.
#'
#' @return A sparse subject-feature-by-reference-feature matrix.
#'
#' @examples
#' reference <- cbind(x = 0:3, y = c(0, 0, 1, 1), z = 0)
#' rownames(reference) <- paste0("r", 1:4)
#' subject <- reference + 0.05
#' rownames(subject) <- paste0("s", 1:4)
#' coordinate_feature_map(reference, subject, k = 2, bandwidth = 0.2)
#' @export
coordinate_feature_map <- function(reference_coords,
                                   subject_coords,
                                   k = 4L,
                                   bandwidth = NULL,
                                   max_distance = Inf,
                                   kernel = c("gaussian", "inverse_distance", "nearest"),
                                   normalize_rows = TRUE,
                                   reference_names = NULL,
                                   subject_names = NULL) {
  reference_coords <- .dmfa_coordinate_matrix(reference_coords, "reference_coords")
  subject_coords <- .dmfa_coordinate_matrix(subject_coords, "subject_coords")
  if (ncol(reference_coords) != ncol(subject_coords)) {
    stop("Reference and subject coordinates must have the same coordinate dimensions.", call. = FALSE)
  }
  if (!is.null(colnames(reference_coords)) && !is.null(colnames(subject_coords))) {
    if (!setequal(colnames(reference_coords), colnames(subject_coords))) {
      stop("Reference and subject coordinates must have the same named coordinate dimensions.", call. = FALSE)
    }
    subject_coords <- subject_coords[, colnames(reference_coords), drop = FALSE]
  }
  if (anyDuplicated(as.data.frame(reference_coords))) {
    stop("`reference_coords` contains duplicate coordinates with ambiguous reference identity.", call. = FALSE)
  }

  kernel <- match.arg(kernel)
  if (!is.numeric(k)) stop("`k` must be one positive integer.", call. = FALSE)
  raw_k <- as.numeric(k)
  if (length(raw_k) != 1L || !is.finite(raw_k) || raw_k < 1L ||
      abs(raw_k - round(raw_k)) > sqrt(.Machine$double.eps)) {
    stop("`k` must be one positive integer.", call. = FALSE)
  }
  k <- as.integer(round(raw_k))
  k_effective <- min(k, nrow(reference_coords))
  if (!is.numeric(max_distance) || length(max_distance) != 1L ||
      is.na(max_distance) || max_distance < 0) {
    stop("`max_distance` must be one nonnegative number or `Inf`.", call. = FALSE)
  }
  if (!is.logical(normalize_rows) || length(normalize_rows) != 1L || is.na(normalize_rows)) {
    stop("`normalize_rows` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(bandwidth) &&
      (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
       !is.finite(bandwidth) || bandwidth <= 0)) {
    stop("`bandwidth` must be NULL or one positive finite number.", call. = FALSE)
  }

  reference_names <- .dmfa_coordinate_names(
    reference_coords,
    reference_names,
    "reference_",
    "reference feature"
  )
  subject_names <- .dmfa_coordinate_names(
    subject_coords,
    subject_names,
    "subject_feature_",
    "subject feature"
  )

  neighbors <- RANN::nn2(
    data = reference_coords,
    query = subject_coords,
    k = k_effective,
    searchtype = "standard"
  )
  index <- matrix(neighbors$nn.idx, nrow = nrow(subject_coords), ncol = k_effective)
  distance <- matrix(neighbors$nn.dists, nrow = nrow(subject_coords), ncol = k_effective)
  allowed <- is.finite(distance) & distance <= max_distance
  if (identical(kernel, "nearest") && k_effective > 1L) {
    allowed[, -1L] <- FALSE
  }
  unsupported <- which(rowSums(allowed) == 0L)
  if (length(unsupported) > 0L) {
    labels <- paste(utils::head(subject_names[unsupported], 5L), collapse = ", ")
    stop(
      sprintf(
        "%d subject coordinate(s) have no reference neighbor within `max_distance` (first: %s).",
        length(unsupported),
        labels
      ),
      call. = FALSE
    )
  }

  if (identical(kernel, "gaussian")) {
    if (is.null(bandwidth)) {
      positive <- distance[allowed & distance > 0]
      bandwidth <- if (length(positive) == 0L) 1 else stats::median(positive)
    }
    weight <- exp(-0.5 * (distance / bandwidth)^2)
  } else if (identical(kernel, "inverse_distance")) {
    weight <- matrix(0, nrow(distance), ncol(distance))
    for (i in seq_len(nrow(distance))) {
      exact <- allowed[i, ] & distance[i, ] == 0
      if (any(exact)) {
        weight[i, exact] <- 1
      } else {
        weight[i, allowed[i, ]] <- 1 / distance[i, allowed[i, ]]
      }
    }
  } else {
    weight <- matrix(1, nrow(distance), ncol(distance))
  }
  weight[!allowed | !is.finite(weight)] <- 0
  zero_weight <- which(rowSums(weight) == 0L)
  if (length(zero_weight) > 0L) {
    nearest_allowed <- max.col(allowed[zero_weight, , drop = FALSE], ties.method = "first")
    weight[cbind(zero_weight, nearest_allowed)] <- 1
  }
  if (isTRUE(normalize_rows)) {
    weight <- weight / rowSums(weight)
  }

  row_index <- rep(seq_len(nrow(subject_coords)), each = k_effective)
  column_index <- as.vector(t(index))
  values <- as.vector(t(weight))
  keep <- values > 0
  out <- Matrix::sparseMatrix(
    i = row_index[keep],
    j = column_index[keep],
    x = values[keep],
    dims = c(nrow(subject_coords), nrow(reference_coords)),
    dimnames = list(subject_names, reference_names)
  )
  attr(out, "kernel") <- kernel
  attr(out, "bandwidth") <- if (identical(kernel, "gaussian")) bandwidth else NULL
  attr(out, "max_distance") <- max_distance
  attr(out, "k") <- k_effective
  attr(out, "nearest_distance") <- distance[, 1L]
  out
}

#' Build Coordinate Maps for Several Subjects
#'
#' @param reference_coords Numeric reference-feature-by-coordinate matrix.
#' @param subject_coords A list of subject coordinate matrices.
#' @param require_full_reference Require every reference feature to receive
#'   support from at least one subject map.
#' @param ... Arguments passed to [coordinate_feature_map()].
#'
#' @return A named list of sparse feature maps. Its `coverage` attribute reports
#'   subject support and aggregate map weight for every reference feature.
#'
#' @examples
#' reference <- cbind(x = 0:3, y = c(0, 0, 1, 1), z = 0)
#' rownames(reference) <- paste0("r", 1:4)
#' subjects <- list(
#'   S1 = structure(reference + 0.03, dimnames = list(paste0("S1_", 1:4), colnames(reference))),
#'   S2 = structure(reference - 0.03, dimnames = list(paste0("S2_", 1:4), colnames(reference)))
#' )
#' maps <- coordinate_feature_maps(subject_coords = subjects,
#'                                 reference_coords = reference,
#'                                 k = 2, bandwidth = 0.2)
#' attr(maps, "coverage")
#' @export
coordinate_feature_maps <- function(reference_coords,
                                    subject_coords,
                                    require_full_reference = TRUE,
                                    ...) {
  if (!is.list(subject_coords) || length(subject_coords) < 1L) {
    stop("`subject_coords` must be a nonempty list of coordinate matrices.", call. = FALSE)
  }
  if (is.null(names(subject_coords))) names(subject_coords) <- paste0("Subject_", seq_along(subject_coords))
  if (anyNA(names(subject_coords)) || any(names(subject_coords) == "") || anyDuplicated(names(subject_coords))) {
    stop("`subject_coords` must have unique, non-empty subject names.", call. = FALSE)
  }
  if (!is.logical(require_full_reference) || length(require_full_reference) != 1L ||
      is.na(require_full_reference)) {
    stop("`require_full_reference` must be TRUE or FALSE.", call. = FALSE)
  }

  maps <- lapply(subject_coords, function(coords) {
    coordinate_feature_map(reference_coords, coords, ...)
  })
  names(maps) <- names(subject_coords)
  reference_names <- colnames(maps[[1L]])
  compatible <- vapply(maps, function(P) identical(colnames(P), reference_names), logical(1))
  if (!all(compatible)) {
    stop("All coordinate maps must resolve to the same reference feature schema.", call. = FALSE)
  }
  subject_support <- Reduce(`+`, lapply(maps, function(P) {
    as.numeric(Matrix::colSums(P != 0) > 0)
  }))
  map_weight <- Reduce(`+`, lapply(maps, function(P) as.numeric(Matrix::colSums(P))))
  coverage <- data.frame(
    feature = reference_names,
    subject_support = subject_support,
    coverage_fraction = subject_support / length(maps),
    map_weight = map_weight,
    stringsAsFactors = FALSE
  )
  if (isTRUE(require_full_reference) && any(subject_support == 0L)) {
    uncovered <- coverage$feature[subject_support == 0L]
    stop(
      sprintf(
        "%d reference feature(s) are not covered by any subject map (first: %s).",
        length(uncovered),
        paste(utils::head(uncovered, 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  attr(maps, "coverage") <- coverage
  maps
}
