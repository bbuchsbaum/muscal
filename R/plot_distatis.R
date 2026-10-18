#' Plot DISTATIS Object and Assessor Maps
#'
#' @param x,object A fitted `distatis` model.
#' @param dims Two distinct retained axes in the chosen map. Each map has its own rank.
#' @param type `"objects"` (alias `"roi"`) for the compromise map, or `"subjects"`
#'   for the assessor principal coordinates in RV space.
#' @param labels Logical; show object or assessor labels?
#' @param col,pch,cex,main Base graphics point color, symbol, size, and title.
#' @param color,alpha,size ggplot point colors, opacity, and size.
#' @param ... Further arguments passed to base graphics, or unused by `autoplot`.
#' @return `plot()` invisibly returns the displayed coordinates; `autoplot()`
#'   returns a ggplot. Axis percentages use the complete spectrum of the chosen
#'   map, including axes not retained in the scores.
#' @seealso [distatis()], [plot_partial_scores()], [plot_block_weights()],
#'   [plot_block_similarity()], [plot_variance()]
#' @md
#' @export
plot.distatis <- function(x, dims = c(1, 2), type = c("objects", "subjects", "roi"),
                          labels = FALSE, col = NULL, pch = 19, cex = 1,
                          main = NULL, ...) {
  map <- .distatis_plot_map(x, match.arg(type), dims)
  args <- list(...)
  args$asp <- args$asp %||% 1
  do.call(.plot_scores_base, c(list(S = map$scores, dims = dims, labels = labels,
                                   col = col, pch = pch, cex = cex,
                                   main = main %||% map$title,
                                   sdev = sqrt(pmax(map$roots, 0))), args))
}

#' @rdname plot.distatis
#' @md
#' @export
autoplot.distatis <- function(object, dims = c(1, 2), type = c("objects", "subjects", "roi"),
                              labels = FALSE, color = NULL, alpha = 0.8, size = 2, ...) {
  .distatis_unused(...)
  map <- .distatis_plot_map(object, match.arg(type), dims)
  .autoplot_scores(map$scores, dims, labels, color, alpha, size,
                   sqrt(pmax(map$roots, 0)), map$title) +
    ggplot2::coord_equal() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.12))
}

.distatis_plot_map <- function(x, type, dims) {
  subjects <- identical(type, "subjects")
  S <- if (subjects) x$subject_scores else x$object_scores
  .distatis_plot_dims(dims, ncol(S))
  list(scores = S,
       roots = if (subjects) x$interstructure_eigenvalues else x$eigenvalues,
       title = if (subjects) "DISTATIS Assessor Map" else "DISTATIS Object Map")
}

.distatis_plot_dims <- function(dims, n) {
  if (!is.numeric(dims) || length(dims) != 2L || any(!is.finite(dims)) ||
      any(dims != floor(dims) | dims < 1 | dims > n) || anyDuplicated(dims)) {
    stop(sprintf("`dims` must select two distinct retained axes; this map has %d.", n), call. = FALSE)
  }
}

#' @rdname plot.distatis
#' @param blocks Optional integer indices of assessors included in the partial map.
#' @param connect,show_consensus Connect partial points to the compromise and show
#'   compromise points?
#' @md
#' @export
plot_partial_scores.distatis <- function(x, dims = c(1, 2), blocks = NULL,
                                         connect = TRUE, show_consensus = TRUE,
                                         alpha = 0.6, ...) {
  .distatis_plot_dims(dims, ncol(x$object_scores))
  if (!is.null(blocks) && (!is.numeric(blocks) || !length(blocks) ||
      any(!is.finite(blocks)) || any(blocks != floor(blocks) | blocks < 1 |
      blocks > length(x$partial_scores)) || anyDuplicated(blocks))) {
    stop("`blocks` must select distinct assessor indices.", call. = FALSE)
  }
  p <- plot_partial_scores.covstatis(x, dims, blocks, connect, show_consensus, alpha, ...)
  p + ggplot2::coord_equal() + ggplot2::labs(title = "DISTATIS Partial Object Maps",
                    subtitle = "Open circles = compromise objects", color = "Assessor")
}

#' @rdname plot.distatis
#' @md
#' @export
plot_block_weights.distatis <- function(x, ...) {
  .plot_block_weights_impl(setNames(x$alpha, x$block_labels),
                           title = "DISTATIS Assessor Weights", ...)
}

#' @rdname plot.distatis
#' @param method Similarity measure. Only `"RV"` is supported for DISTATIS.
#' @md
#' @export
plot_block_similarity.distatis <- function(x, method = c("RV", "RV2"), ...) {
  if (match.arg(method) != "RV") stop("DISTATIS block similarity uses RV, not RV2.", call. = FALSE)
  .plot_heatmap(x$C, title = "DISTATIS Assessor RV Similarities", ...)
}

#' DISTATIS Inertia Plot
#'
#' @param x A fitted `distatis` model.
#' @param type `"bar"` or `"line"` (with cumulative inertia).
#' @param axes Plot `"objects"` or `"subjects"` inertia.
#' @param ... Additional arguments passed to [plot_variance()].
#' @return A ggplot of the complete positive spectrum, including axes not retained
#'   in the score maps. Call through `plot_variance(fit, axes = "subjects")`.
#' @md
#' @export
plot_variance.distatis <- function(x, type = c("bar", "line"),
                                   axes = c("objects", "subjects"), ...) {
  axes <- match.arg(axes)
  roots <- if (axes == "subjects") x$interstructure_eigenvalues else x$eigenvalues
  roots <- roots[roots > max(roots) * x$tol]
  plot_variance(list(sdev = sqrt(roots)), type = match.arg(type), ...) +
    ggplot2::labs(title = paste("DISTATIS", if (axes == "subjects") "Assessor" else "Object", "Inertia"),
                  y = "% Inertia")
}

#' @md
#' @export
print.distatis <- function(x, ...) {
  cat(sprintf("DISTATIS: %d objects, %d assessors\n", nrow(x$v), length(x$alpha)))
  cat(sprintf("Retained axes: %d object, %d assessor\n", ncol(x$v), ncol(x$subject_scores)))
  cat(sprintf("Normalization: %s; compromise quality: %.1f%%\n",
              x$norm_method, 100 * x$compromise_quality))
  invisible(x)
}
