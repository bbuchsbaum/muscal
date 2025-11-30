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
#' @param alpha Confidence level for bootstrap intervals.
#' @param boot_loadings For \code{loading_reliability}: a list of loading
#'   matrices from bootstrap refits (all same dimensions as the reference
#'   loadings).
#' @return
#'   \item{estimate_components}{A list with \code{keep}, \code{criterion},
#'   and \code{method_used}.}
#'   \item{loading_reliability}{A list with matrices \code{mean}, \code{sd},
#'   \code{lower}, \code{upper}, and \code{sign_consistency}.}
#' @export
estimate_components <- function(fit,
                                method = c("rmt", "variance"),
                                sdev = NULL,
                                V_list = NULL,
                                n = NULL,
                                tail_q = 0.2,
                                alpha = 0.05) {
  method <- match.arg(method)
  sdev <- sdev %||% fit$sdev %||% stop("No singular values available.")

  V_list <- V_list %||% .extract_V_list(fit)
  n <- n %||% .infer_n(fit)

  if (method == "rmt" && !is.null(n) && !is.null(V_list)) {
    k_vec <- rep(ncol(V_list[[1]]), length(V_list))
    res <- significant_components(
      list(V_list = V_list, sdev = sdev),
      n = n,
      k_vec = k_vec,
      check_rmt = TRUE,
      tail_quantile = tail_q
    )
    return(list(keep = res$keep, criterion = res$rmt_pass, method_used = "rmt"))
  }

  # Fallback: keep all positive singular values
  keep <- which(sdev > 0)
  list(keep = keep, criterion = sdev, method_used = "variance")
}

#' @rdname estimate_components
#' @export
loading_reliability <- function(fit,
                                method = c("bootstrap"),
                                boot_loadings = NULL,
                                alpha = 0.05,
                                V_ref = NULL) {
  method <- match.arg(method)
  if (method == "bootstrap") {
    if (is.null(boot_loadings) || length(boot_loadings) == 0) {
      stop("boot_loadings must be provided for bootstrap reliability.")
    }
    V_ref <- V_ref %||% .extract_V_list(fit)
    # Flatten list-of-matrices bootstrap into summarize_boot-like inputs
    boot_mean <- Reduce("+", boot_loadings) / length(boot_loadings)
    arr <- abind::abind(boot_loadings, along = 3)
    sd_mat <- apply(arr, c(1, 2), stats::sd)
    lower <- apply(arr, c(1, 2), stats::quantile, probs = alpha / 2, type = 1)
    upper <- apply(arr, c(1, 2), stats::quantile, probs = 1 - alpha / 2, type = 1)
    sign_cons <- apply(sign(arr), c(1, 2), function(v) {
      max(mean(v > 0), mean(v < 0))
    })
    return(list(
      mean = unname(boot_mean),
      sd = unname(sd_mat),
      lower = unname(lower),
      upper = unname(upper),
      sign_consistency = unname(sign_cons)
    ))
  }
  stop("Unknown method for loading_reliability.")
}

# --- internal helpers ---------------------------------------------------------
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
