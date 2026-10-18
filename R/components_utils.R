#' Estimate the number of components and loading reliability (shared helpers)
#'
#' These helpers provide a common interface for component retention and loading
#' reliability across \code{bada}, \code{mfa}, \code{penalized_mfa}, and related
#' models. They are intentionally modular so alternative criteria (e.g., AIC,
#' permutations) can be added later without changing call sites.
#'
#' @param fit A fitted object (e.g., \code{bada}, \code{mfa},
#'   \code{penalized_mfa}, \code{penalized_mfa_clusterwise}).
#' @param method Component or reliability method. For \code{estimate_components},
#'   one of \code{"rmt"} (Marchenko–Pastur edge) or \code{"variance"}
#'   (keep all non-zero). For \code{loading_reliability}, currently
#'   \code{"bootstrap"}.
#' @param sdev Optional singular values; if \code{NULL}, taken from
#'   \code{fit$sdev} when available.
#' @param V_list Optional list of loading matrices; if \code{NULL}, the helper
#'   attempts to extract from \code{fit} (using \code{attr(fit, "V_list")} or
#'   splitting \code{fit$v} by \code{block_indices} when present).
#' @param n Optional number of observations (rows); defaults to \code{nrow}
#'   of scores if available.
#' @param tail_q Numeric; fraction of smallest eigenvalues used to estimate
#'   noise variance in the historical automatic RMT heuristic. Retained for
#'   compatibility; automatic noise estimation is now unavailable.
#' @param noise_variance,spectrum,noise_model,n_eff Supplied-noise RMT model and
#'   spectral units; see [significant_components()]. No arbitrary weighted
#'   compromise spectrum is assumed to satisfy the isotropic covariance model.
#' @param alpha Tail probability for bootstrap intervals in `loading_reliability`.
#'   In `estimate_components`, a significance threshold only for explicitly
#'   enabled experimental ICC. Explicit alpha on an edge-only or variance-only
#'   call warns that it is inapplicable; it does not change an MP edge.
#' @param check_icc,common_features Optional experimental common-feature ICC
#'   criterion and required feature identity contract; see [significant_components()].
#' @param boot_loadings For \code{loading_reliability}: a list of loading
#'   matrices from bootstrap refits (all same dimensions as the reference
#'   loadings).
#' @param V_ref For \code{loading_reliability}: optional reference loading
#'   matrix. If \code{NULL}, extracted from \code{fit} using the same logic
#'   as \code{V_list}.
#' @param gap_tol For loading reliability, a nonnegative tolerance for gaps
#'   between squared reference component scales, normalized by the largest scale
#'   squared. Reference `fit$sdev` must follow the supplied reference column order.
#' @return
#'   \item{estimate_components}{A list with \code{keep}, \code{criterion},
#'   and \code{method_used}. RMT additionally returns \code{diagnostics} with
#'   feature dimensions, applicability and criterion details. Missing model or
#'   observation information yields \code{method_used="rmt_unavailable"}, NA
#'   decisions and no retained components, never a silent variance fallback.}
#'   \item{loading_reliability}{A list with matrices \code{mean}, \code{sd},
#'   \code{lower}, \code{upper}, and \code{sign_consistency}, plus \code{alignment}
#'   diagnostics. Exact one-to-one matching uses absolute normalized loading inner
#'   products, followed by sign alignment. Known tied/zero reference roots have
#'   unavailable (NA) coordinate summaries and geometric subspace distances.
#'   Without reference component scales, intervals are conditional on the supplied
#'   reference alignment and identifiability is reported as unknown. Only the
#'   retained reference gaps are checked; omitted-axis gaps and bootstrap coverage
#'   are not established by these summaries.}
#' @export
estimate_components <- function(fit,
                                method = c("rmt", "variance"),
                                sdev = NULL,
                                V_list = NULL,
                                n = NULL,
                                tail_q = 0.2,
                                alpha = 0.05,
                                noise_variance = NULL,
                                spectrum = c("unknown", "covariance", "singular"),
                                noise_model = c("unknown", "isotropic"),
                                n_eff = n,
                                check_icc = FALSE,
                                common_features = NULL) {
  method <- match.arg(method)
  chk::chk_flag(check_icc)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a finite scalar strictly between 0 and 1.", call. = FALSE)
  }
  if (!missing(alpha) && (method == "variance" || !check_icc)) {
    warning("`alpha` applies only to explicitly enabled ICC; a support edge is not an alpha-level test.",
            call. = FALSE)
  }
  sdev <- sdev %||% fit$sdev %||% stop("No singular values available.")

  if (method == "variance") {
    return(list(keep = which(sdev > 0), criterion = sdev, method_used = "variance"))
  }

  V_list <- V_list %||% .extract_V_list(fit)
  n <- n %||% .infer_n(fit)
  if (!is.list(V_list) || !length(V_list) ||
      any(!vapply(V_list, function(V) is.matrix(V) && is.numeric(V) &&
                    nrow(V) > 0L && ncol(V) > 0L && all(is.finite(V)), logical(1)))) {
    stop("V_list must contain nonempty finite numeric loading matrices.", call. = FALSE)
  }
  k_vec <- vapply(V_list, nrow, integer(1))
  if (!is.null(fit$block_indices) &&
      (!is.list(fit$block_indices) ||
       !identical(unname(lengths(fit$block_indices)), unname(k_vec)))) {
    stop("Stored block indices disagree with loading feature counts.", call. = FALSE)
  }
  if (!is.null(n)) {
    args <- list(
      fit = list(V_list = V_list, sdev = sdev),
      n = n,
      k_vec = k_vec,
      check_rmt = TRUE,
      tail_frac = tail_q,
      noise_variance = noise_variance,
      spectrum = match.arg(spectrum),
      noise_model = match.arg(noise_model),
      n_eff = n_eff,
      check_icc = check_icc,
      common_features = common_features
    )
    if (check_icc) args$alpha <- alpha
    res <- do.call(significant_components, args)
    criterion <- if (check_icc) res$rmt_pass & res$icc_pass else res$rmt_pass
    return(list(keep = res$keep, criterion = criterion,
                method_used = if (res$rmt_status == "evaluated") {
                  if (check_icc) "rmt+icc" else "rmt"
                } else "rmt_unavailable",
                diagnostics = res))
  }

  list(keep = integer(0), criterion = rep(NA, length(sdev)),
       method_used = "rmt_unavailable",
       diagnostics = list(k_vec = k_vec, rmt_status = "unavailable",
                          rmt_reason = "Observation count is unavailable."))
}

#' @rdname estimate_components
#' @export
loading_reliability <- function(fit,
                                method = c("bootstrap"),
                                boot_loadings = NULL,
                                alpha = 0.05,
                                V_ref = NULL,
                                gap_tol = 1e-6) {
  method <- match.arg(method)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be strictly between zero and one.", call. = FALSE)
  }
  if (!is.numeric(gap_tol) || length(gap_tol) != 1L || !is.finite(gap_tol) || gap_tol < 0) {
    stop("`gap_tol` must be finite and nonnegative.", call. = FALSE)
  }
  if (method == "bootstrap") {
    if (is.null(boot_loadings) || length(boot_loadings) == 0) {
      stop("boot_loadings must be provided for bootstrap reliability.")
    }
    V_ref <- V_ref %||% .extract_V_list(fit)
    V_ref_mat <- .muscal_as_component_basis_matrix(V_ref, what = "V_ref")
    if (length(boot_loadings) < 2L) stop("At least two bootstrap loading matrices are required.", call. = FALSE)
    .muscal_unit_loading_basis(V_ref_mat)
    aligned <- lapply(seq_along(boot_loadings), function(i) {
      boot_i <- .muscal_as_component_basis_matrix(boot_loadings[[i]],
        ncomp = ncol(V_ref_mat), what = sprintf("boot_loadings[[%d]]", i))
      .muscal_align_loading_basis(V_ref_mat, boot_i)
    })
    boot_mats <- lapply(aligned, `[[`, "basis")
    K <- ncol(V_ref_mat)
    ambiguous <- rep(FALSE, K)
    groups <- list(seq_len(K))
    identifiability <- "unknown_reference_spectrum"
    if (!is.null(fit$sdev)) {
      roots <- fit$sdev
      if (!is.numeric(roots) || !is.null(dim(roots)) || length(roots) < K ||
          any(!is.finite(roots)) || any(roots < 0)) {
        stop("Reference spectrum must be a vector of finite nonnegative component scales.", call. = FALSE)
      }
      roots <- roots[seq_len(K)]
      roots <- if (max(roots) > 0) (roots / max(roots))^2 else roots
      ordering <- order(roots, decreasing = TRUE)
      groups <- unname(split(ordering, cumsum(c(1L, abs(diff(roots[ordering])) > gap_tol))))
      ambiguous <- roots == 0
      for (group in groups) if (length(group) > 1L) ambiguous[group] <- TRUE
      groups <- Filter(function(group) any(ambiguous[group]), groups)
      identifiability <- if (any(ambiguous)) "ambiguous_reference_roots" else "retained_reference_roots_separated"
    }
    distances <- matrix(0, length(boot_mats), length(groups))
    for (j in seq_along(groups)) {
      group <- groups[[j]]
      distances[, j] <- vapply(boot_mats, function(B) .muscal_basis_subspace_distance(
        V_ref_mat[, group, drop = FALSE], B[, group, drop = FALSE]), numeric(1))
    }

    boot_mean <- Reduce("+", boot_mats) / length(boot_mats)
    arr <- abind::abind(boot_mats, along = 3)
    sd_mat <- apply(arr, c(1, 2), stats::sd)
    lower <- apply(arr, c(1, 2), stats::quantile, probs = alpha / 2, type = 1)
    upper <- apply(arr, c(1, 2), stats::quantile, probs = 1 - alpha / 2, type = 1)
    sign_cons <- apply(sign(arr), c(1, 2), function(v) {
      max(mean(v > 0), mean(v < 0))
    })
    summaries <- list(mean = unname(boot_mean), sd = unname(sd_mat),
      lower = unname(lower), upper = unname(upper), sign_consistency = unname(sign_cons))
    if (any(ambiguous)) summaries <- lapply(summaries, function(X) {
      X[, ambiguous] <- NA_real_
      X
    })
    summaries$alignment <- list(
      permutations = lapply(aligned, `[[`, "permutation"),
      signs = lapply(aligned, `[[`, "signs"),
      ambiguous_components = ambiguous, identifiability = identifiability,
      subspace_groups = groups, subspace_distance = distances)
    return(summaries)
  }
  stop("Unknown method for loading_reliability.")
}

# --- internal helpers ---------------------------------------------------------
.muscal_unit_loading_basis <- function(V) {
  if (!is.matrix(V) || !is.numeric(V) || !nrow(V) || !ncol(V) || any(!is.finite(V))) {
    stop("Loading bases must be nonempty finite numeric matrices.", call. = FALSE)
  }
  amplitude <- apply(abs(V), 2L, max)
  if (any(amplitude == 0)) stop("Cannot align zero-norm loading directions.", call. = FALSE)
  V <- sweep(V, 2L, amplitude, "/")
  sweep(V, 2L, sqrt(colSums(V^2)), "/")
}

# Exact square assignment by augmenting paths and dual potentials, O(K^3).
.muscal_min_assignment <- function(cost) {
  if (!is.matrix(cost) || !is.numeric(cost) || !nrow(cost) ||
      nrow(cost) != ncol(cost) || any(!is.finite(cost))) {
    stop("Assignment costs must be a nonempty finite square matrix.", call. = FALSE)
  }
  n <- nrow(cost)
  if (n == 1L) return(1L)
  u <- v <- numeric(n + 1L)
  matched <- integer(n + 1L)
  for (i in seq_len(n)) {
    matched[1L] <- i
    j0 <- 1L
    minimum <- rep(Inf, n + 1L)
    used <- rep(FALSE, n + 1L)
    previous <- integer(n + 1L)
    repeat {
      used[j0] <- TRUE
      i0 <- matched[j0]
      delta <- Inf
      j1 <- 0L
      for (j in seq.int(2L, n + 1L)) if (!used[j]) {
        reduced <- cost[i0, j - 1L] - u[i0 + 1L] - v[j]
        if (reduced < minimum[j]) {
          minimum[j] <- reduced
          previous[j] <- j0
        }
        if (minimum[j] < delta) {
          delta <- minimum[j]
          j1 <- j
        }
      }
      for (j in seq_len(n + 1L)) {
        if (used[j]) {
          u[matched[j] + 1L] <- u[matched[j] + 1L] + delta
          v[j] <- v[j] - delta
        } else minimum[j] <- minimum[j] - delta
      }
      j0 <- j1
      if (matched[j0] == 0L) break
    }
    repeat {
      j1 <- previous[j0]
      matched[j0] <- matched[j1]
      j0 <- j1
      if (j0 == 1L) break
    }
  }
  assignment <- integer(n)
  assignment[matched[-1L]] <- seq_len(n)
  assignment
}

.muscal_align_loading_basis <- function(reference, candidate) {
  if (!identical(dim(reference), dim(candidate))) {
    stop("Bootstrap loading dimensions must match the reference basis.", call. = FALSE)
  }
  if (!is.null(rownames(reference)) && !identical(rownames(reference), rownames(candidate))) {
    if (is.null(rownames(candidate)) || anyDuplicated(rownames(reference)) ||
        anyDuplicated(rownames(candidate)) || !setequal(rownames(reference), rownames(candidate))) {
      stop("Bootstrap feature identities must match the reference basis.", call. = FALSE)
    }
    candidate <- candidate[match(rownames(reference), rownames(candidate)), , drop = FALSE]
  }
  R <- .muscal_unit_loading_basis(reference)
  B <- .muscal_unit_loading_basis(candidate)
  similarity <- crossprod(R, B)
  similarity[] <- pmax(-1, pmin(1, similarity))
  permutation <- .muscal_min_assignment(1 - abs(similarity))
  signs <- .muscal_component_signs_from_basis(B[, permutation, drop = FALSE], reference = R)
  list(basis = sweep(candidate[, permutation, drop = FALSE], 2L, signs, "*"),
       permutation = permutation, signs = signs)
}

.muscal_basis_subspace_distance <- function(A, B) {
  thin <- function(X) {
    q <- qr(.muscal_unit_loading_basis(X))
    qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  }
  A <- thin(A)
  B <- thin(B)
  # Projector Frobenius distance via residuals, avoiding P x P matrices and
  # cancellation when the subspaces nearly agree.
  sqrt(sum((B - A %*% crossprod(A, B))^2) +
       sum((A - B %*% crossprod(B, A))^2))
}

.extract_V_list <- function(fit) {
  if (!is.null(attr(fit, "V_list"))) return(attr(fit, "V_list"))
  if (!is.null(fit$V_list)) return(fit$V_list)
  if (!is.null(fit$v) && !is.null(fit$block_indices)) {
    v <- fit$v
    idx <- fit$block_indices
    return(lapply(idx, function(ii) v[ii, , drop = FALSE]))
  }
  stop("Unable to extract V_list from fit; please supply V_list explicitly.")
}

.infer_n <- function(fit) {
  if (!is.null(fit$s)) return(nrow(fit$s))
  if (!is.null(fit$roi_scores)) return(nrow(fit$roi_scores))
  if (!is.null(fit$subject_scores)) return(nrow(fit$subject_scores))
  NULL
}
