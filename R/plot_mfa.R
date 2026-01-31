#' @title Plotting Functions for MFA and Linked MFA
#' @name mfa-plots
#' @description
#' Visualization methods for Multiple Factor Analysis (MFA) and Linked MFA results.
#' Includes score plots, partial factor scores, block contributions, convergence
#' diagnostics, and more.
#'
#' @importFrom graphics plot points lines legend abline barplot image axis box par text segments symbols
#' @importFrom stats cor sd aggregate
#' @importFrom grDevices hcl.colors adjustcolor
#' @importFrom rlang .data
NULL

# =============================================================================
# Core score plots (plot and autoplot)
# =============================================================================

#' Plot MFA Scores
#'
#' Displays observations in the MFA compromise space.
#'
#' @param x An object of class `mfa`.
#' @param dims Integer vector of length 2 specifying which components to plot.
#' @param labels Logical; if `TRUE`, label points with row names or indices.
#' @param col Point colors. If `NULL`, uses a default.
#' @param pch Point character.
#' @param cex Point size.
#' @param main Plot title.
#' @param ... Additional arguments passed to [plot()].
#'
#' @return Invisibly returns the score matrix (for the plotted dimensions).
#' @export
#' @examples
#' \donttest{
#' X <- replicate(3, matrix(rnorm(50 * 10), 50, 10), simplify = FALSE)
#' fit <- mfa(X, ncomp = 3)
#' plot(fit)
#' }
plot.mfa <- function(x, dims = c(1, 2), labels = FALSE, col = NULL,
                     pch = 19, cex = 1, main = "MFA Score Plot", ...) {

  S <- multivarious::scores(x)
  .plot_scores_base(S, dims = dims, labels = labels, col = col,
                    pch = pch, cex = cex, main = main,
                    sdev = x$sdev, ...)
}

#' Plot Linked MFA Scores
#'
#' Displays observations in the Linked MFA shared score space, optionally
#' indicating block coverage.
#'
#' @param x An object of class `linked_mfa`.
#' @param dims Integer vector of length 2 specifying which components to plot.
#' @param show_coverage Logical; if `TRUE`, color points by number of blocks
#'   contributing to each observation.
#' @param labels Logical; if `TRUE`, label points with row names or indices.
#' @param col Point colors. Overrides `show_coverage` coloring if provided.
#' @param pch Point character.
#' @param cex Point size.
#' @param main Plot title.
#' @param ... Additional arguments passed to [plot()].
#'
#' @return Invisibly returns the score matrix (for the plotted dimensions).
#' @export
plot.linked_mfa <- function(x, dims = c(1, 2), show_coverage = TRUE,
                            labels = FALSE, col = NULL, pch = 19, cex = 1,
                            main = "Linked MFA Score Plot", ...) {
  S <- multivarious::scores(x)
  N <- nrow(S)

  coverage <- NULL
  cov_levels <- NULL
  pal <- NULL
  show_cov_legend <- FALSE

  if (show_coverage && is.null(col)) {
    # Count how many X blocks cover each observation (Y counts as 1)
    coverage <- rep.int(1L, N)
    for (k in seq_along(x$row_index)) {
      idx <- unique(x$row_index[[k]])
      coverage[idx] <- coverage[idx] + 1L
    }
    cov_levels <- sort(unique(coverage))
    pal <- hcl.colors(length(cov_levels), palette = "viridis")
    col <- pal[match(coverage, cov_levels)]
    show_cov_legend <- TRUE
  }

  .plot_scores_base(S, dims = dims, labels = labels, col = col,
                    pch = pch, cex = cex, main = main,
                    sdev = x$sdev, ...)

  if (show_cov_legend) {
    legend(
      "topright",
      legend = paste0(cov_levels, " block(s)"),
      fill = pal,
      title = "Coverage",
      bty = "n",
      cex = 0.8
    )
  }
}

#' @keywords internal
.plot_scores_base <- function(S, dims, labels, col, pch, cex, main, sdev, ...) {
  if (is.null(col)) col <- "steelblue"
  d1 <- dims[1]
  d2 <- dims[2]

  var_exp <- if (!is.null(sdev)) {
    round(100 * sdev^2 / sum(sdev^2), 1)
  } else {
    rep(NA, ncol(S))
  }

  xlab <- if (!is.na(var_exp[d1])) {
    sprintf("Component %d (%.1f%%)", d1, var_exp[d1])
  } else {
    sprintf("Component %d", d1)
  }
  ylab <- if (!is.na(var_exp[d2])) {
    sprintf("Component %d (%.1f%%)", d2, var_exp[d2])
  } else {
    sprintf("Component %d", d2)
  }

  plot(S[, d1], S[, d2], col = col, pch = pch, cex = cex,
       xlab = xlab, ylab = ylab, main = main, ...)
  abline(h = 0, v = 0, col = "gray70", lty = 2)


  if (labels) {
    labs <- rownames(S)
    if (is.null(labs)) labs <- seq_len(nrow(S))
    text(S[, d1], S[, d2], labels = labs, pos = 3, cex = cex * 0.7)
  }

  invisible(S[, dims, drop = FALSE])
}

# =============================================================================
# autoplot methods (ggplot2-based)
# =============================================================================

#' Autoplot for MFA
#'
#' Creates a ggplot2 score plot for MFA results.
#'
#' @param object An object of class `mfa`.
#' @param dims Integer vector of length 2 specifying components to plot.
#' @param labels Logical; if `TRUE`, add text labels.
#' @param color Optional variable for point coloring (vector of length n).
#' @param alpha Point transparency.
#' @param size Point size.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 autoplot
autoplot.mfa <- function(object, dims = c(1, 2), labels = FALSE,
                         color = NULL, alpha = 0.8, size = 2, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot.", call. = FALSE)
  }

  S <- multivarious::scores(object)
  .autoplot_scores(S, dims = dims, labels = labels, color = color,
                   alpha = alpha, size = size, sdev = object$sdev,
                   title = "MFA Score Plot")
}

#' Autoplot for Linked MFA
#'
#' Creates a ggplot2 score plot for Linked MFA results.
#'
#' @param object An object of class `linked_mfa`.
#' @param dims Integer vector of length 2 specifying components to plot.
#' @param labels Logical; if `TRUE`, add text labels.
#' @param color Optional variable for point coloring. If `"coverage"`, colors
#'   by number of blocks contributing to each observation.
#' @param alpha Point transparency.
#' @param size Point size.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object.
#' @export
autoplot.linked_mfa <- function(object, dims = c(1, 2), labels = FALSE,
                                color = "coverage", alpha = 0.8, size = 2, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot.", call. = FALSE)
  }

  S <- multivarious::scores(object)
  N <- nrow(S)

  if (identical(color, "coverage")) {
    coverage <- rep.int(1L, N)
    for (k in seq_along(object$row_index)) {
      idx <- unique(object$row_index[[k]])
      coverage[idx] <- coverage[idx] + 1L
    }
    color <- factor(coverage, levels = sort(unique(coverage)))
  }

  .autoplot_scores(S, dims = dims, labels = labels, color = color,
                   alpha = alpha, size = size, sdev = object$sdev,
                   title = "Linked MFA Score Plot")
}

#' @keywords internal
.autoplot_scores <- function(S, dims, labels, color, alpha, size, sdev, title) {
  d1 <- dims[1]
  d2 <- dims[2]

  var_exp <- if (!is.null(sdev)) {
    round(100 * sdev^2 / sum(sdev^2), 1)
  } else {
    rep(NA, ncol(S))
  }

  xlab <- if (!is.na(var_exp[d1])) {
    sprintf("Component %d (%.1f%%)", d1, var_exp[d1])
  } else {
    sprintf("Component %d", d1)
  }
  ylab <- if (!is.na(var_exp[d2])) {
    sprintf("Component %d (%.1f%%)", d2, var_exp[d2])
  } else {
    sprintf("Component %d", d2)
  }

  df <- data.frame(
    x = S[, d1],
    y = S[, d2],
    label = if (!is.null(rownames(S))) rownames(S) else seq_len(nrow(S))
  )
  if (!is.null(color)) df$color <- color

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y))

  if (!is.null(color)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$color),
                                  alpha = alpha, size = size)
  } else {
    p <- p + ggplot2::geom_point(alpha = alpha, size = size, color = "steelblue")
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(x = xlab, y = ylab, title = title, color = NULL) +
    ggplot2::theme_minimal()

  if (labels) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$label),
                                 vjust = -0.5, size = 3)
  }

  p
}

# =============================================================================
# Partial factor scores
# =============================================================================

#' Plot Partial Factor Scores
#'
#' Displays how each block "sees" the observations. For MFA, this shows

#' block-specific projections. For Linked MFA, connects auxiliary block
#' observations to their reference positions.
#'
#' @param x A fitted `mfa` or `linked_mfa` object.
#' @param dims Integer vector of length 2 for component axes.
#' @param blocks Which blocks to include. `NULL` means all blocks.
#' @param connect Logical; if `TRUE`, draw lines connecting the same
#'   observation across blocks.
#' @param show_consensus Logical; if `TRUE`, show the global scores as larger points.
#' @param alpha Transparency for partial score points.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @details
#' For MFA: Each block's partial scores are computed by projecting that block
#' alone onto the compromise space. Observations where partial scores diverge
#' indicate block disagreement.
#'
#' For Linked MFA
#' the auxiliary block observations are shown at their
#' projected positions, with lines connecting to the reference score.
#'
#' @return A ggplot object (if ggplot2 available) or base R plot.
#' @export
plot_partial_scores <- function(x, dims = c(1, 2), blocks = NULL,
                                connect = TRUE, show_consensus = TRUE,
                                alpha = 0.6, ...) {
  UseMethod("plot_partial_scores")
}

#' @export
plot_partial_scores.mfa <- function(x, dims = c(1, 2), blocks = NULL,
                                    connect = TRUE, show_consensus = TRUE,
                                    alpha = 0.6, ...) {
  S <- multivarious::scores(x)
  n <- nrow(S)
  d1 <- dims[1]
  d2 <- dims[2]

  if (is.null(blocks)) blocks <- seq_along(x$block_indices)
  block_names <- x$names
  if (is.null(block_names)) block_names <- names(x$block_indices)
  if (is.null(block_names)) block_names <- paste0("Block", seq_along(x$block_indices))
  block_labels <- block_names[blocks]

  if (is.null(x$partial_scores)) {
    stop(
      "Partial scores are not available in this MFA object. ",
      "Refit with a newer muscal version to enable plot_partial_scores().",
      call. = FALSE
    )
  }
  partial_list <- x$partial_scores[blocks]

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R fallback
    return(.plot_partial_base(S, partial_list, dims, blocks, block_labels,
                              connect, show_consensus, alpha, ...))
  }

  # Build data frame
  df_list <- lapply(seq_along(blocks), function(i) {
    ps <- partial_list[[i]]
    data.frame(
      x = ps[, d1],
      y = ps[, d2],
      obs = seq_len(n),
      block = block_labels[i]
    )
  })
  df_partial <- do.call(rbind, df_list)

  df_global <- data.frame(
    x = S[, d1],
    y = S[, d2],
    obs = seq_len(n),
    block = "Consensus"
  )

  p <- ggplot2::ggplot()

  if (connect) {
    # Lines connecting partial to consensus
    for (i in seq_along(blocks)) {
      df_lines <- data.frame(
        x = partial_list[[i]][, d1],
        y = partial_list[[i]][, d2],
        xend = S[, d1],
        yend = S[, d2],
        block = block_labels[i]
      )
      p <- p + ggplot2::geom_segment(
        data = df_lines,
        ggplot2::aes(x = .data$x, y = .data$y,
                     xend = .data$xend, yend = .data$yend,
                     color = .data$block),
        alpha = alpha * 0.5, linewidth = 0.3
      )
    }
  }

  p <- p + ggplot2::geom_point(
    data = df_partial,
    ggplot2::aes(x = .data$x, y = .data$y, color = .data$block),
    alpha = alpha, size = 2
  )

  if (show_consensus) {
    p <- p + ggplot2::geom_point(
      data = df_global,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "black", size = 3, shape = 1
    )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      x = sprintf("Component %d", d1),
      y = sprintf("Component %d", d2),
      title = "Partial Factor Scores",
      color = "Block"
    ) +
    ggplot2::theme_minimal()

  p
}

#' @export
plot_partial_scores.linked_mfa <- function(x, dims = c(1, 2), blocks = NULL,
                                           connect = TRUE, show_consensus = TRUE,
                                           alpha = 0.6, ...) {
  S <- multivarious::scores(x)
  N <- nrow(S)
  d1 <- dims[1]
  d2 <- dims[2]

  if (is.null(blocks)) blocks <- seq_along(x$V_list)
  block_names <- names(x$V_list)
  if (is.null(block_names)) block_names <- paste0("X", blocks)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for plot_partial_scores.linked_mfa", call. = FALSE
)
  }

  if (is.null(x$partial_scores)) {
    stop(
      "Partial scores are not available in this linked_mfa object. ",
      "Refit with a newer muscal version to enable plot_partial_scores().",
      call. = FALSE
    )
  }

  # For each X block, plot its score estimate from X_k and connect it to
  # the reference score S[idx_k, ].
  df_list <- list()
  df_lines_list <- list()

  for (i in blocks) {
    nm <- block_names[i]
    idx_k <- if (!is.null(names(x$row_index)) && nm %in% names(x$row_index)) {
      x$row_index[[nm]]
    } else {
      x$row_index[[i]]
    }

    ps_hat <- x$partial_scores[[nm]]
    if (is.null(ps_hat)) {
      stop("No partial scores found for block '", nm, "'.", call. = FALSE)
    }

    df_list[[length(df_list) + 1L]] <- data.frame(
      x = ps_hat[, d1],
      y = ps_hat[, d2],
      ref_row = idx_k,
      block = nm
    )

    if (connect) {
      df_lines_list[[length(df_lines_list) + 1L]] <- data.frame(
        x = ps_hat[, d1],
        y = ps_hat[, d2],
        xend = S[idx_k, d1],
        yend = S[idx_k, d2],
        block = nm
      )
    }
  }

  df_partial <- do.call(rbind, df_list)
  df_lines <- if (connect) do.call(rbind, df_lines_list) else NULL

  df_global <- data.frame(
    x = S[, d1],
    y = S[, d2],
    obs = seq_len(N)
  )

  p <- ggplot2::ggplot()

  if (connect) {
    p <- p + ggplot2::geom_segment(
      data = df_lines,
      ggplot2::aes(
        x = .data$x, y = .data$y,
        xend = .data$xend, yend = .data$yend,
        color = .data$block
      ),
      alpha = alpha * 0.5, linewidth = 0.3
    )
  }

  p <- p + ggplot2::geom_point(
    data = df_partial,
    ggplot2::aes(x = .data$x, y = .data$y, color = .data$block),
    alpha = alpha, size = 2
  )

  if (show_consensus) {
    p <- p + ggplot2::geom_point(
      data = df_global,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "gray30", size = 2.5, shape = 1, alpha = 0.5
    )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      x = sprintf("Component %d", d1),
      y = sprintf("Component %d", d2),
      title = "Partial Factor Scores (Linked MFA)",
      subtitle = "Open circles = reference positions",
      color = "Block"
    ) +
    ggplot2::theme_minimal()

  p
}

#' @keywords internal
.plot_partial_base <- function(S, partial_list, dims, blocks, block_names,
                               connect, show_consensus, alpha, ...) {
  d1 <- dims[1]
  d2 <- dims[2]
  n <- nrow(S)

  all_x <- c(S[, d1], unlist(lapply(partial_list, function(p) p[, d1])))
  all_y <- c(S[, d2], unlist(lapply(partial_list, function(p) p[, d2])))

  n_blocks <- length(partial_list)
  cols <- hcl.colors(n_blocks, palette = "Set2")

  plot(NULL, xlim = range(all_x), ylim = range(all_y),
       xlab = sprintf("Component %d", d1),
       ylab = sprintf("Component %d", d2),
       main = "Partial Factor Scores")
  abline(h = 0, v = 0, col = "gray70", lty = 2)

  for (i in seq_len(n_blocks)) {
    ps <- partial_list[[i]]
    if (connect) {
      segments(ps[, d1], ps[, d2], S[, d1], S[, d2],
               col = adjustcolor(cols[i], alpha.f = alpha * 0.5), lwd = 0.5)
    }
    points(ps[, d1], ps[, d2], col = adjustcolor(cols[i], alpha.f = alpha), pch = 19)
  }

  if (show_consensus) {
    points(S[, d1], S[, d2], pch = 1, cex = 1.2, lwd = 1.5)
  }

  legend("topright", legend = block_names, fill = cols, bty = "n", cex = 0.8)
  invisible(NULL)
}

# =============================================================================
# Block weights
# =============================================================================

#' Plot Block Weights
#'
#' Bar chart showing the normalization weights (alpha) for each block.
#'
#' @param x A fitted `mfa` or `linked_mfa` object.
#' @param ... Additional arguments passed to [barplot()] or ggplot.
#'
#' @return A ggplot object or base R plot.
#' @export
plot_block_weights <- function(x, ...) {
  UseMethod("plot_block_weights")
}

#' @export
plot_block_weights.mfa <- function(x, ...) {
  alpha <- x$alpha
  names(alpha) <- x$names
  if (is.null(names(alpha))) names(alpha) <- paste0("Block", seq_along(alpha))
  .plot_block_weights_impl(alpha, title = "MFA Block Weights", ...)
}

#' @export
plot_block_weights.linked_mfa <- function(x, ...) {
  alpha <- x$alpha_blocks
  .plot_block_weights_impl(alpha, title = "Linked MFA Block Weights", ...)
}

#' @keywords internal
.plot_block_weights_impl <- function(alpha, title, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    barplot(alpha, main = title, ylab = "Weight", col = "steelblue",
            las = 2, ...)
    return(invisible(alpha))
  }

  df <- data.frame(
    block = factor(names(alpha), levels = names(alpha)),
    weight = as.numeric(alpha)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$block, y = .data$weight)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.8) +
    ggplot2::labs(x = "Block", y = "Weight", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p
}

# =============================================================================
# Block similarity (RV matrix)
# =============================================================================

#' Plot Block Similarity Matrix
#'
#' Heatmap of pairwise RV coefficients between blocks, showing which data
#' sources tell similar stories.
#'
#' @param x A fitted `mfa` object.
#' @param method Similarity method: `"RV"` or `"RV2"`.
#' @param ... Additional arguments.
#'
#' @return A ggplot object or base R image plot.
#' @export
plot_block_similarity <- function(x, method = c("RV", "RV2"), ...) {

  UseMethod("plot_block_similarity")
}
#' @export
plot_block_similarity.mfa <- function(x, method = c("RV", "RV2"), ...) {
  method <- match.arg(method)

  sim_mat <- if (method == "RV") x$rv else x$rv2
  if (is.null(sim_mat)) {
    stop(
      "RV matrix not available in this MFA object. ",
      "Refit with a newer muscal version to enable plot_block_similarity().",
      call. = FALSE
    )
  }

  .plot_heatmap(sim_mat, title = sprintf("MFA %s Matrix (Block Similarity)", method), ...)
}

#' @keywords internal
.plot_heatmap <- function(mat, title, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    image(t(mat[nrow(mat):1, ]), main = title,
          col = hcl.colors(50, "YlOrRd", rev = TRUE),
          axes = FALSE, ...)
    axis(1, at = seq(0, 1, length.out = ncol(mat)),
         labels = colnames(mat), las = 2, cex.axis = 0.8)
    axis(2, at = seq(0, 1, length.out = nrow(mat)),
         labels = rev(rownames(mat)), las = 1, cex.axis = 0.8)
    box()
    return(invisible(mat))
  }

  df <- expand.grid(
    row = factor(rownames(mat), levels = rownames(mat)),
    col = factor(colnames(mat), levels = colnames(mat))
  )
  df$value <- as.vector(mat)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$value)),
                        color = "white", size = 3) +
    ggplot2::scale_fill_viridis_c(option = "magma", direction = -1) +
    ggplot2::labs(x = NULL, y = NULL, title = title, fill = "Similarity") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed()

  p
}

# =============================================================================
# Variance / Scree plot
# =============================================================================

#' Plot Variance Explained (Scree Plot)
#'
#' Shows the proportion of variance captured by each component.
#'
#' @param x A fitted `mfa` or `linked_mfa` object.
#' @param type `"bar"` for bar chart, `"line"` for scree plot with cumulative.
#' @param ... Additional arguments.
#'
#' @return A ggplot object or base R plot.
#' @export
plot_variance <- function(x, type = c("bar", "line"), ...) {
  type <- match.arg(type)
  sdev <- x$sdev
  if (is.null(sdev) || length(sdev) == 0) {
    stop("No sdev found in the fitted object.", call. = FALSE)
  }

  var_exp <- sdev^2
  pct <- 100 * var_exp / sum(var_exp)
  cum_pct <- cumsum(pct)
  ncomp <- length(sdev)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    if (type == "bar") {
      barplot(pct, names.arg = seq_len(ncomp),
              main = "Variance Explained", xlab = "Component",
              ylab = "% Variance", col = "steelblue")
    } else {
      plot(seq_len(ncomp), pct, type = "b", pch = 19,
           main = "Scree Plot", xlab = "Component", ylab = "% Variance")
      lines(seq_len(ncomp), cum_pct, type = "b", pch = 17, col = "red")
      legend("right", legend = c("Individual", "Cumulative"),
             pch = c(19, 17), col = c("black", "red"), bty = "n")
    }
    return(invisible(data.frame(component = seq_len(ncomp), pct = pct, cumulative = cum_pct)))
  }

  df <- data.frame(
    component = seq_len(ncomp),
    pct = pct,
    cumulative = cum_pct
  )

  if (type == "bar") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(.data$component), y = .data$pct)) +
      ggplot2::geom_col(fill = "steelblue", alpha = 0.8) +
      ggplot2::labs(x = "Component", y = "% Variance", title = "Variance Explained") +
      ggplot2::theme_minimal()
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$component)) +
      ggplot2::geom_line(ggplot2::aes(y = .data$pct), color = "steelblue", linewidth = 1) +
      ggplot2::geom_point(ggplot2::aes(y = .data$pct), color = "steelblue", size = 3) +
      ggplot2::geom_line(ggplot2::aes(y = .data$cumulative), color = "firebrick",
                          linewidth = 1, linetype = "dashed") +
      ggplot2::geom_point(ggplot2::aes(y = .data$cumulative), color = "firebrick",
                           size = 2, shape = 17) +
      ggplot2::scale_x_continuous(breaks = seq_len(ncomp)) +
      ggplot2::labs(x = "Component", y = "% Variance", title = "Scree Plot") +
      ggplot2::annotate("text", x = ncomp, y = max(cum_pct) - 5,
                         label = "Cumulative", color = "firebrick", hjust = 1) +
      ggplot2::theme_minimal()
  }

  p
}

# =============================================================================
# Loadings plot
# =============================================================================

#' Plot Variable Loadings
#'
#' Displays variable contributions to components, optionally as a correlation
#' circle or bar chart.
#'
#' @param x A fitted `mfa` or `linked_mfa` object.
#' @param dims Integer vector of length 2 for component axes.
#' @param type `"circle"` for correlation circle, `"bar"` for bar chart of
#'   loadings on a single component.
#' @param component For `type = "bar"`, which component to show.
#' @param color_by For MFA, color variables by `"block"` or `NULL`.
#' @param top_n For `type = "bar"`, show only top N variables by absolute loading.
#' @param ... Additional arguments.
#'
#' @return A ggplot object or base R plot.
#' @export
plot_loadings <- function(x, dims = c(1, 2), type = c("circle", "bar"),
                          component = 1, color_by = "block", top_n = 20, ...) {
type <- match.arg(type)
  UseMethod("plot_loadings")
}

#' @export
plot_loadings.mfa <- function(x, dims = c(1, 2), type = c("circle", "bar"),
                              component = 1, color_by = "block", top_n = 20, ...) {
  type <- match.arg(type)
  V <- x$v
  if (type == "circle") {
    V <- .circle_loadings_matrix(x, fallback = V)
  }
  d1 <- dims[1]
  d2 <- dims[2]

  # Assign block labels to each variable
  block_labels <- character(nrow(V))
  block_names <- x$names
  if (is.null(block_names)) block_names <- paste0("Block", seq_along(x$block_indices))
  for (i in seq_along(x$block_indices)) {
    block_labels[x$block_indices[[i]]] <- block_names[i]
  }

  var_names <- rownames(V)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(nrow(V)))

  .plot_loadings_impl(V, dims, type, component, color_by, top_n,
                      block_labels, var_names, title_prefix = "MFA", ...)
}

#' @export
plot_loadings.linked_mfa <- function(x, dims = c(1, 2), type = c("circle", "bar"),
                                     component = 1, color_by = "block", top_n = 20, ...) {
  type <- match.arg(type)
  V <- x$v
  if (type == "circle") {
    V <- .circle_loadings_matrix(x, fallback = V)
  }
  d1 <- dims[1]
  d2 <- dims[2]

  # Assign block labels
  block_labels <- character(nrow(V))
  block_names <- names(x$block_indices)
  if (is.null(block_names)) block_names <- paste0("Block", seq_along(x$block_indices))
  for (i in seq_along(x$block_indices)) {
    block_labels[x$block_indices[[i]]] <- block_names[i]
  }

  var_names <- rownames(V)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(nrow(V)))

  .plot_loadings_impl(V, dims, type, component, color_by, top_n,
                      block_labels, var_names, title_prefix = "Linked MFA", ...)
}

#' @keywords internal
.circle_loadings_matrix <- function(x, fallback) {
  if (!is.null(x$cor_loadings)) {
    return(as.matrix(x$cor_loadings))
  }

  # For genpca-based objects (e.g., MFA), fall back to correlation-style
  # coordinates derived from the unscaled right vectors when available.
  if (!is.null(x$ov) && !is.null(x$sdev)) {
    n <- nrow(multivarious::scores(x))
    denom <- sqrt(max(1, n - 1))
    return(sweep(as.matrix(x$ov), 2, x$sdev / denom, `*`))
  }

  as.matrix(fallback)
}

#' @keywords internal
.plot_loadings_impl <- function(V, dims, type, component, color_by, top_n,
                                block_labels, var_names, title_prefix, ...) {
  d1 <- dims[1]
  d2 <- dims[2]

  V <- as.matrix(V)
  V[!is.finite(V)] <- 0

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    if (type == "circle") {
      coords <- V[, dims, drop = FALSE]
      norms <- sqrt(rowSums(coords^2))
      max_norm <- max(norms, na.rm = TRUE)
      if (is.finite(max_norm) && max_norm > 1) {
        coords <- coords / max_norm
      }

      plot(coords[, 1], coords[, 2], type = "n",
           xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), asp = 1,
           xlab = sprintf("Component %d", d1),
           ylab = sprintf("Component %d", d2),
           main = paste(title_prefix, "Loadings"))
      segments(0, 0, coords[, 1], coords[, 2],
               col = adjustcolor("steelblue", alpha.f = 0.7))
      symbols(0, 0, circles = 1, add = TRUE, inches = FALSE, fg = "gray")
      abline(h = 0, v = 0, col = "gray70", lty = 2)
    } else {
      loadings_comp <- V[, component]
      ord <- order(abs(loadings_comp), decreasing = TRUE)
      top_idx <- ord[seq_len(min(top_n, length(ord)))]
      barplot(loadings_comp[top_idx], names.arg = var_names[top_idx],
              main = sprintf("%s Loadings (Component %d)", title_prefix, component),
              las = 2, col = "steelblue")
    }
    return(invisible(NULL))
  }

  if (type == "circle") {
    # Correlation circle
    coords <- V[, dims, drop = FALSE]
    norms <- sqrt(rowSums(coords^2))
    max_norm <- max(norms, na.rm = TRUE)
    if (is.finite(max_norm) && max_norm > 1) {
      coords <- coords / max_norm
    }

    df <- data.frame(
      x = coords[, 1],
      y = coords[, 2],
      variable = var_names,
      block = block_labels
    )

    # Circle data
    theta <- seq(0, 2 * pi, length.out = 100)
    circle_df <- data.frame(x = cos(theta), y = sin(theta))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y))

    if (!is.null(color_by) && color_by == "block") {
      p <- p + ggplot2::geom_segment(
        ggplot2::aes(x = 0, y = 0, xend = .data$x, yend = .data$y, color = .data$block),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")),
        alpha = 0.7
      )
    } else {
      p <- p + ggplot2::geom_segment(
        ggplot2::aes(x = 0, y = 0, xend = .data$x, yend = .data$y),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")),
        color = "steelblue", alpha = 0.7
      )
    }

    p <- p +
      ggplot2::geom_path(data = circle_df, ggplot2::aes(x = .data$x, y = .data$y),
                          color = "gray50", linetype = "dashed") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
      ggplot2::coord_fixed() +
      ggplot2::labs(
        x = sprintf("Component %d", d1),
        y = sprintf("Component %d", d2),
        title = paste(title_prefix, "Loadings (Correlation Circle)"),
        color = if (!is.null(color_by) && color_by == "block") "Block" else NULL
      ) +
      ggplot2::theme_minimal()

  } else {
    # Bar chart
    loadings_comp <- V[, component]
    ord <- order(abs(loadings_comp), decreasing = TRUE)
    top_idx <- ord[seq_len(min(top_n, length(ord)))]

    df <- data.frame(
      variable = factor(var_names[top_idx], levels = var_names[top_idx]),
      loading = loadings_comp[top_idx],
      block = block_labels[top_idx]
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$variable, y = .data$loading))

    if (!is.null(color_by) && color_by == "block") {
      p <- p + ggplot2::geom_col(ggplot2::aes(fill = .data$block), alpha = 0.8)
    } else {
      p <- p + ggplot2::geom_col(fill = "steelblue", alpha = 0.8)
    }

    p <- p +
      ggplot2::geom_hline(yintercept = 0, color = "gray40") +
      ggplot2::labs(
        x = "Variable",
        y = "Loading",
        title = sprintf("%s Loadings (Component %d)", title_prefix, component),
        fill = if (!is.null(color_by) && color_by == "block") "Block" else NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  p
}

# =============================================================================
# Convergence trace (Linked MFA specific)
# =============================================================================

#' Plot Convergence Trace
#'
#' Shows the objective function value over iterations for iterative methods.
#'
#' @param x A fitted object with an `objective_trace` (e.g., `linked_mfa`,
#'   `bamfa`, `penalized_mfa`).
#' @param log_scale Logical; if `TRUE`, use log scale for y-axis.
#' @param ... Additional arguments.
#'
#' @return A ggplot object or base R plot.
#' @export
plot_convergence <- function(x, log_scale = FALSE, ...) {
  UseMethod("plot_convergence")
}

#' @rdname plot_convergence
#' @export
plot_convergence.linked_mfa <- function(x, log_scale = FALSE, ...) {
  obj <- x$objective_trace
  if (is.null(obj) || length(obj) == 0) {
    stop("No objective_trace found in the fitted object.", call. = FALSE)
  }

  n_iter <- length(obj)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    if (log_scale) {
      plot(seq_len(n_iter), log(obj), type = "b", pch = 19,
           xlab = "Iteration", ylab = "log(Objective)",
           main = "Linked MFA Convergence")
    } else {
      plot(seq_len(n_iter), obj, type = "b", pch = 19,
           xlab = "Iteration", ylab = "Objective",
           main = "Linked MFA Convergence")
    }
    return(invisible(obj))
  }

  df <- data.frame(iteration = seq_len(n_iter), objective = obj)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$objective)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(color = "steelblue", size = 2) +
    ggplot2::labs(
      x = "Iteration",
      y = if (log_scale) "log(Objective)" else "Objective",
      title = "Linked MFA Convergence"
    ) +
    ggplot2::theme_minimal()

  if (log_scale) {
    p <- p + ggplot2::scale_y_log10()
  }

  p
}

# =============================================================================
# Coverage map (Linked MFA specific)
# =============================================================================

#' Plot Block Coverage Map
#'
#' Heatmap showing which observations are present in which blocks.
#'
#' @param x A fitted `linked_mfa` object.
#' @param ... Additional arguments.
#'
#' @return A ggplot object or base R plot.
#' @export
plot_coverage <- function(x, ...) {
  if (!inherits(x, "linked_mfa")) {
    stop("plot_coverage is only available for linked_mfa objects.", call. = FALSE)
  }

  N <- nrow(multivarious::scores(x))
  n_blocks <- 1L + length(x$row_index)
  block_names <- names(x$block_indices)
  if (is.null(block_names)) block_names <- c("Y", paste0("X", seq_along(x$row_index)))

  # Build coverage matrix
  cov_mat <- matrix(0L, nrow = N, ncol = n_blocks)
  colnames(cov_mat) <- block_names
  rownames(cov_mat) <- seq_len(N)

  cov_mat[, 1] <- 1L
  for (k in seq_along(x$row_index)) {
    cov_mat[x$row_index[[k]], k + 1L] <- 1L
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    image(t(cov_mat), col = c("gray90", "steelblue"),
          main = "Block Coverage", axes = FALSE)
    axis(1, at = seq(0, 1, length.out = n_blocks), labels = block_names, las = 2)
    axis(2, at = seq(0, 1, length.out = N), labels = seq_len(N), las = 1, cex.axis = 0.5)
    box()
    return(invisible(cov_mat))
  }

  df <- expand.grid(
    observation = seq_len(N),
    block = factor(block_names, levels = block_names)
  )
  df$present <- as.vector(cov_mat)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$block, y = .data$observation,
                                         fill = factor(.data$present))) +
    ggplot2::geom_tile(color = "white", linewidth = 0.1) +
    ggplot2::scale_fill_manual(values = c("0" = "gray90", "1" = "steelblue"),
                                labels = c("Missing", "Present"),
                                name = NULL) +
    ggplot2::labs(x = "Block", y = "Observation", title = "Block Coverage Map") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = if (N > 50) ggplot2::element_blank() else ggplot2::element_text(size = 6),
      panel.grid = ggplot2::element_blank()
    )

  p
}

# =============================================================================
# Feature groups (Linked MFA specific)
# =============================================================================

#' Plot Feature Group Loadings
#'
#' For Linked MFA models with feature groups, shows how grouped features
#' align in the loading space.
#'
#' @param x A fitted `linked_mfa` object with `feature_groups`.
#' @param dims Integer vector of length 2 for component axes.
#' @param groups Which groups to show. `NULL` means all groups.
#' @param show_centers Logical; if `TRUE`, show group centers.
#' @param ... Additional arguments.
#'
#' @return A ggplot object or base R plot.
#' @export
plot_feature_groups <- function(x, dims = c(1, 2), groups = NULL,
                                show_centers = TRUE, ...) {
  if (!inherits(x, "linked_mfa")) {
    stop("plot_feature_groups is only available for linked_mfa objects.", call. = FALSE)
  }

  fg <- x$feature_groups
  if (is.null(fg) || identical(fg, FALSE) || x$feature_lambda == 0) {
    stop("This linked_mfa model was not fit with feature groups.", call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for plot_feature_groups.", call. = FALSE)
  }

  d1 <- dims[1]
  d2 <- dims[2]
  V_list <- x$V_list
  block_names <- names(V_list)
  if (is.null(block_names)) block_names <- paste0("X", seq_along(V_list))

  # Parse feature groups to find members
  # This depends on how feature_groups was stored
  if (is.data.frame(fg)) {
    group_col <- fg$group
    all_groups <- unique(group_col)
    if (!is.null(groups)) all_groups <- intersect(all_groups, groups)

    df_list <- list()
    for (g in all_groups) {
      rows_g <- fg[fg$group == g, , drop = FALSE]
      for (i in seq_len(nrow(rows_g))) {
        block_idx <- rows_g$block[i]
        if (is.character(block_idx)) block_idx <- match(block_idx, block_names)
        feat_idx <- rows_g$feature[i]
        if (is.character(feat_idx)) {
          feat_idx <- match(feat_idx, colnames(V_list[[block_idx]]))
        }
        v <- V_list[[block_idx]][feat_idx, dims]
        df_list[[length(df_list) + 1]] <- data.frame(
          x = v[1], y = v[2],
          group = g,
          block = block_names[block_idx],
          feature = feat_idx
        )
      }
    }
    df <- do.call(rbind, df_list)

  } else if (identical(fg, "colnames")) {
    # Reconstruct from colnames
    cn_all <- lapply(V_list, function(V) rownames(V))
    all_names <- unique(unlist(cn_all))

    df_list <- list()
    for (nm in all_names) {
      in_blocks <- which(sapply(cn_all, function(cn) nm %in% cn))
      if (length(in_blocks) < 2) next
      for (b in in_blocks) {
        idx <- match(nm, cn_all[[b]])
        v <- V_list[[b]][idx, dims]
        df_list[[length(df_list) + 1]] <- data.frame(
          x = v[1], y = v[2],
          group = nm,
          block = block_names[b],
          feature = idx
        )
      }
    }
    df <- do.call(rbind, df_list)
    if (!is.null(groups)) df <- df[df$group %in% groups, ]
  } else {
    stop("Unrecognized feature_groups format.", call. = FALSE)
  }

  if (nrow(df) == 0) {
    stop("No feature groups found to plot.", call. = FALSE)
  }

  # Compute centers
  centers <- aggregate(cbind(x, y) ~ group, data = df, FUN = mean)
  group_levels <- unique(df$group)
  df$group <- factor(df$group, levels = group_levels)
  centers$group <- factor(centers$group, levels = group_levels)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$group)) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$block), size = 2.5, alpha = 0.75)

  # Connect members of same group
  for (g in unique(df$group)) {
    df_g <- df[df$group == g, ]
    if (nrow(df_g) > 1) {
      cx <- centers$x[centers$group == g]
      cy <- centers$y[centers$group == g]
      p <- p + ggplot2::geom_segment(
        data = df_g,
        ggplot2::aes(xend = cx, yend = cy),
        alpha = 0.25, linetype = "dashed",
        show.legend = FALSE
      )
    }
  }

  if (show_centers) {
    p <- p + ggplot2::geom_point(
      data = centers,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$group),
      size = 4, shape = 4, stroke = 1.5,
      show.legend = FALSE
    )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      x = sprintf("Component %d", d1),
      y = sprintf("Component %d", d2),
      title = "Feature Group Loadings",
      color = "Group",
      shape = "Block"
    ) +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::theme_minimal()

  p
}

# =============================================================================
# Block fit / reconstruction error
# =============================================================================

#' Plot Block Reconstruction Error
#'
#' Bar chart showing per-block reconstruction error, indicating which blocks
#' are well-captured by the shared scores.
#'
#' @param x A fitted `linked_mfa` object.
#' @param normalize Logical; if `TRUE`, show relative (percentage) error.
#' @param ... Additional arguments.
#'
#' @return A ggplot object or base R plot.
#' @export
plot_block_fit <- function(x, normalize = TRUE, ...) {
  if (!inherits(x, "linked_mfa")) {
    stop("plot_block_fit is only available for linked_mfa objects.", call. = FALSE)
  }

  bf <- x$block_fit
  if (is.null(bf) || !is.data.frame(bf) || nrow(bf) == 0) {
    stop(
      "Block fit metrics are not available in this linked_mfa object. ",
      "Refit with a newer muscal version to enable plot_block_fit().",
      call. = FALSE
    )
  }

  block_levels <- bf$block
  metric <- if (isTRUE(normalize)) "r2" else "sse"
  y <- if (metric == "r2") 100 * bf$r2 else bf$sse

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    barplot(
      y,
      names.arg = block_levels,
      main = "Block Fit",
      ylab = if (metric == "r2") "R2 (%)" else "SSE",
      col = "steelblue",
      las = 2,
      ...
    )
    return(invisible(bf))
  }

  df <- data.frame(
    block = factor(block_levels, levels = block_levels),
    value = y
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$block, y = .data$value)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.8) +
    ggplot2::labs(
      x = "Block",
      y = if (metric == "r2") "R2 (%)" else "SSE",
      title = "Linked MFA Block Fit",
      subtitle = if (metric == "r2") "Explained variance by block (higher is better)" else "Reconstruction error by block"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p
}

# =============================================================================
# COVSTATIS plotting methods
# =============================================================================

#' Plot COVSTATIS Results
#'
#' Displays ROI scores or subject G-scores from a COVSTATIS analysis.
#'
#' @param x A fitted `covstatis` object.
#' @param dims Integer vector of length 2 specifying which components to plot.
#' @param type `"roi"` for ROI-level scores, `"subjects"` for subject G-scores.
#' @param labels Logical; if `TRUE`, label points.
#' @param col Point colors.
#' @param pch Point character.
#' @param cex Point size.
#' @param main Plot title.
#' @param ... Additional arguments passed to [plot()].
#'
#' @return Invisibly returns the plotted scores.
#' @export
plot.covstatis <- function(x, dims = c(1, 2), type = c("roi", "subjects"),
                           labels = FALSE, col = NULL, pch = 19, cex = 1,
                           main = NULL, ...) {
  type <- match.arg(type)

  if (type == "roi") {
    S <- x$roi_scores
    if (is.null(main)) main <- "COVSTATIS ROI Scores"
    point_labels <- x$labels
  } else {
    # G-scores: barycentric mean of partial scores per subject
    S <- do.call(rbind, lapply(x$partial_scores, colMeans))
    if (is.null(main)) main <- "COVSTATIS Subject G-Scores"
    point_labels <- x$block_labels
  }

  if (is.null(col)) col <- "steelblue"

  .plot_scores_base(S, dims = dims, labels = labels, col = col,
                    pch = pch, cex = cex, main = main,
                    sdev = x$sdev, ...)

  if (labels && !is.null(point_labels)) {
    text(S[, dims[1]], S[, dims[2]], labels = point_labels, pos = 3, cex = cex * 0.7)
  }

  invisible(S[, dims, drop = FALSE])
}

#' Autoplot for COVSTATIS
#'
#' Creates a ggplot2 visualization of COVSTATIS results.
#'
#' @param object A fitted `covstatis` object.
#' @param dims Integer vector of length 2 specifying components to plot.
#' @param type `"roi"` for ROI-level scores, `"subjects"` for subject G-scores.
#' @param labels Logical; if `TRUE`, add text labels.
#' @param color Optional variable for point coloring.
#' @param alpha Point transparency.
#' @param size Point size.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object.
#' @export
autoplot.covstatis <- function(object, dims = c(1, 2), type = c("roi", "subjects"),
                               labels = FALSE, color = NULL, alpha = 0.8, size = 2, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot.", call. = FALSE)
  }
  type <- match.arg(type)

  if (type == "roi") {
    S <- object$roi_scores
    title <- "COVSTATIS ROI Scores"
    point_labels <- object$labels
  } else {
    S <- do.call(rbind, lapply(object$partial_scores, colMeans))
    title <- "COVSTATIS Subject G-Scores"
    point_labels <- object$block_labels
  }

  .autoplot_scores(S, dims = dims, labels = labels, color = color,
                   alpha = alpha, size = size, sdev = object$sdev,
                   title = title)
}

#' @export
plot_block_weights.covstatis <- function(x, ...) {

  alpha <- x$alpha
  names(alpha) <- x$block_labels
  .plot_block_weights_impl(alpha, title = "COVSTATIS Block Weights (Alpha)", ...)
}

#' @export
plot_partial_scores.covstatis <- function(x, dims = c(1, 2), blocks = NULL,
                                          connect = TRUE, show_consensus = TRUE,
                                          alpha = 0.6, ...) {
  d1 <- dims[1]
  d2 <- dims[2]

  # Consensus: ROI scores
  S_consensus <- x$roi_scores

  if (is.null(blocks)) blocks <- seq_along(x$partial_scores)
  block_names <- x$block_labels[blocks]

  partial_list <- x$partial_scores[blocks]

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(.plot_partial_base(S_consensus, partial_list, dims, blocks, block_names,
                              connect, show_consensus, alpha, ...))
  }

  n <- nrow(S_consensus)
  df_list <- lapply(seq_along(blocks), function(i) {
    ps <- partial_list[[i]]
    data.frame(
      x = ps[, d1],
      y = ps[, d2],
      obs = seq_len(n),
      block = block_names[i]
    )
  })
  df_partial <- do.call(rbind, df_list)

  df_global <- data.frame(
    x = S_consensus[, d1],
    y = S_consensus[, d2],
    obs = seq_len(n)
  )

  p <- ggplot2::ggplot()

  if (connect) {
    for (i in seq_along(blocks)) {
      df_lines <- data.frame(
        x = partial_list[[i]][, d1],
        y = partial_list[[i]][, d2],
        xend = S_consensus[, d1],
        yend = S_consensus[, d2],
        block = block_names[i]
      )
      p <- p + ggplot2::geom_segment(
        data = df_lines,
        ggplot2::aes(x = .data$x, y = .data$y,
                     xend = .data$xend, yend = .data$yend,
                     color = .data$block),
        alpha = alpha * 0.5, linewidth = 0.3
      )
    }
  }

  p <- p + ggplot2::geom_point(
    data = df_partial,
    ggplot2::aes(x = .data$x, y = .data$y, color = .data$block),
    alpha = alpha, size = 2
  )

  if (show_consensus) {
    p <- p + ggplot2::geom_point(
      data = df_global,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "black", size = 3, shape = 1
    )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      x = sprintf("Component %d", d1),
      y = sprintf("Component %d", d2),
      title = "COVSTATIS Partial Factor Scores",
      subtitle = "Open circles = compromise (ROI scores)",
      color = "Subject"
    ) +
    ggplot2::theme_minimal()

  p
}

#' @export
plot_block_similarity.covstatis <- function(x, method = c("RV", "RV2"), ...) {
  # COVSTATIS already has the RV matrix computed as C
  sim_mat <- x$C
  block_names <- x$block_labels
  if (!is.null(block_names)) {
    rownames(sim_mat) <- colnames(sim_mat) <- block_names
  }

  .plot_heatmap(sim_mat, title = "COVSTATIS RV Matrix (Block Similarity)", ...)
}

#' @rdname plot_variance
#' @export
plot_variance.covstatis <- function(x, type = c("bar", "line"), ...) {
  type <- match.arg(type)
  sdev <- x$sdev
  if (is.null(sdev) || length(sdev) == 0) {
    stop("No sdev found in the fitted object.", call. = FALSE)
  }

  var_exp <- sdev^2
  pct <- 100 * var_exp / sum(var_exp)
  cum_pct <- cumsum(pct)
  ncomp <- length(sdev)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    if (type == "bar") {
      barplot(pct, names.arg = seq_len(ncomp),
              main = "COVSTATIS Variance Explained", xlab = "Component",
              ylab = "% Variance", col = "steelblue")
    } else {
      plot(seq_len(ncomp), pct, type = "b", pch = 19,
           main = "COVSTATIS Scree Plot", xlab = "Component", ylab = "% Variance")
      lines(seq_len(ncomp), cum_pct, type = "b", pch = 17, col = "red")
      legend("right", legend = c("Individual", "Cumulative"),
             pch = c(19, 17), col = c("black", "red"), bty = "n")
    }
    return(invisible(data.frame(component = seq_len(ncomp), pct = pct, cumulative = cum_pct)))
  }

  df <- data.frame(
    component = seq_len(ncomp),
    pct = pct,
    cumulative = cum_pct
  )

  if (type == "bar") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(.data$component), y = .data$pct)) +
      ggplot2::geom_col(fill = "steelblue", alpha = 0.8) +
      ggplot2::labs(x = "Component", y = "% Variance", title = "COVSTATIS Variance Explained") +
      ggplot2::theme_minimal()
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$component)) +
      ggplot2::geom_line(ggplot2::aes(y = .data$pct), color = "steelblue", linewidth = 1) +
      ggplot2::geom_point(ggplot2::aes(y = .data$pct), color = "steelblue", size = 3) +
      ggplot2::geom_line(ggplot2::aes(y = .data$cumulative), color = "firebrick",
                          linewidth = 1, linetype = "dashed") +
      ggplot2::geom_point(ggplot2::aes(y = .data$cumulative), color = "firebrick",
                           size = 2, shape = 17) +
      ggplot2::scale_x_continuous(breaks = seq_len(ncomp)) +
      ggplot2::labs(x = "Component", y = "% Variance", title = "COVSTATIS Scree Plot") +
      ggplot2::theme_minimal()
  }

  p
}

# =============================================================================
# BAMFA plotting methods
# =============================================================================

#' Plot BaMFA Results
#'
#' Displays global scores from a BaMFA analysis.
#'
#' @param x A fitted `bamfa` object.
#' @param dims Integer vector of length 2 specifying which components to plot.
#' @param block Which block's scores to show, or `"mean"` for averaged scores.
#' @param labels Logical; if `TRUE`, label points.
#' @param col Point colors.
#' @param pch Point character.
#' @param cex Point size.
#' @param main Plot title.
#' @param ... Additional arguments passed to [plot()].
#'
#' @return Invisibly returns the plotted scores.
#' @export
plot.bamfa <- function(x, dims = c(1, 2), block = "mean",
                       labels = FALSE, col = NULL, pch = 19, cex = 1,
                       main = "BaMFA Global Scores", ...) {

  if (identical(block, "mean")) {
    # Average scores across all blocks
    S <- Reduce("+", x$S_list) / length(x$S_list)
  } else {
    S <- x$S_list[[block]]
  }

  if (is.null(col)) col <- "steelblue"

  .plot_scores_base(S, dims = dims, labels = labels, col = col,
                    pch = pch, cex = cex, main = main,
                    sdev = NULL, ...)
}

#' Autoplot for BaMFA
#'
#' Creates a ggplot2 visualization of BaMFA global scores.
#'
#' @param object A fitted `bamfa` object.
#' @param dims Integer vector of length 2 specifying components to plot.
#' @param block Which block's scores to show, or `"mean"` for averaged scores.
#' @param labels Logical; if `TRUE`, add text labels.
#' @param color Optional variable for point coloring.
#' @param alpha Point transparency.
#' @param size Point size.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object.
#' @export
autoplot.bamfa <- function(object, dims = c(1, 2), block = "mean",
                           labels = FALSE, color = NULL, alpha = 0.8, size = 2, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot.", call. = FALSE)
  }

  if (identical(block, "mean")) {
    S <- Reduce("+", object$S_list) / length(object$S_list)
    title <- "BaMFA Global Scores (Mean)"
  } else {
    S <- object$S_list[[block]]
    title <- sprintf("BaMFA Global Scores (Block %s)", block)
  }

  .autoplot_scores(S, dims = dims, labels = labels, color = color,
                   alpha = alpha, size = size, sdev = NULL,
                   title = title)
}

#' @rdname plot_convergence
#' @export
plot_convergence.bamfa <- function(x, log_scale = FALSE, ...) {
  obj <- x$objective_trace
  if (is.null(obj) || length(obj) == 0) {
    stop("No objective_trace found in the fitted object.", call. = FALSE)
  }

  n_iter <- length(obj)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    if (log_scale) {
      plot(seq_len(n_iter), log(obj), type = "b", pch = 19,
           xlab = "Iteration", ylab = "log(MSE)",
           main = "BaMFA Convergence")
    } else {
      plot(seq_len(n_iter), obj, type = "b", pch = 19,
           xlab = "Iteration", ylab = "MSE",
           main = "BaMFA Convergence")
    }
    return(invisible(obj))
  }

  df <- data.frame(iteration = seq_len(n_iter), objective = obj)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$objective)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(color = "steelblue", size = 2) +
    ggplot2::labs(
      x = "Iteration",
      y = if (log_scale) "log(MSE)" else "MSE",
      title = "BaMFA Convergence"
    ) +
    ggplot2::theme_minimal()

  if (log_scale) {
    p <- p + ggplot2::scale_y_log10()
  }

  p
}

#' @export
plot_partial_scores.bamfa <- function(x, dims = c(1, 2), blocks = NULL,
                                      connect = TRUE, show_consensus = TRUE,
                                      alpha = 0.6, ...) {
  d1 <- dims[1]
  d2 <- dims[2]

  # Consensus: mean of global scores
  S_consensus <- Reduce("+", x$S_list) / length(x$S_list)
  n <- nrow(S_consensus)

  if (is.null(blocks)) blocks <- seq_along(x$S_list)
  block_names <- x$data_names[blocks]
  if (is.null(block_names)) block_names <- paste0("Block", blocks)

  partial_list <- x$S_list[blocks]

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(.plot_partial_base(S_consensus, partial_list, dims, blocks, block_names,
                              connect, show_consensus, alpha, ...))
  }

  df_list <- lapply(seq_along(blocks), function(i) {
    ps <- partial_list[[i]]
    data.frame(
      x = ps[, d1],
      y = ps[, d2],
      obs = seq_len(n),
      block = block_names[i]
    )
  })
  df_partial <- do.call(rbind, df_list)

  df_global <- data.frame(
    x = S_consensus[, d1],
    y = S_consensus[, d2],
    obs = seq_len(n)
  )

  p <- ggplot2::ggplot()

  if (connect) {
    for (i in seq_along(blocks)) {
      df_lines <- data.frame(
        x = partial_list[[i]][, d1],
        y = partial_list[[i]][, d2],
        xend = S_consensus[, d1],
        yend = S_consensus[, d2],
        block = block_names[i]
      )
      p <- p + ggplot2::geom_segment(
        data = df_lines,
        ggplot2::aes(x = .data$x, y = .data$y,
                     xend = .data$xend, yend = .data$yend,
                     color = .data$block),
        alpha = alpha * 0.5, linewidth = 0.3
      )
    }
  }

  p <- p + ggplot2::geom_point(
    data = df_partial,
    ggplot2::aes(x = .data$x, y = .data$y, color = .data$block),
    alpha = alpha, size = 2
  )

  if (show_consensus) {
    p <- p + ggplot2::geom_point(
      data = df_global,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "black", size = 3, shape = 1
    )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      x = sprintf("Component %d", d1),
      y = sprintf("Component %d", d2),
      title = "BaMFA Block Scores",
      subtitle = "Open circles = consensus (mean scores)",
      color = "Block"
    ) +
    ggplot2::theme_minimal()

  p
}

#' Plot Global vs Local Decomposition for BaMFA
#'
#' Visualizes how variance is partitioned between global and local components.
#'
#' @param x A fitted `bamfa` object.
#' @param ... Additional arguments.
#'
#' @return A ggplot object or base R plot.
#' @export
plot_global_local <- function(x, ...) {
  if (!inherits(x, "bamfa")) {
    stop("plot_global_local is only available for bamfa objects.", call. = FALSE)
  }

  # Compute variance explained by global vs local for each block
  n_blocks <- length(x$S_list)
  block_names <- x$data_names
  if (is.null(block_names)) block_names <- paste0("Block", seq_len(n_blocks))

  var_global <- numeric(n_blocks)
  var_local <- numeric(n_blocks)

  for (i in seq_len(n_blocks)) {
    G_i <- x$v[x$block_indices[[i]], , drop = FALSE]
    S_i <- x$S_list[[i]]
    U_i <- x$U_list[[i]]
    B_i <- x$B_list[[i]]

    # Variance from global reconstruction
    recon_global <- S_i %*% t(G_i)
    var_global[i] <- sum(recon_global^2)

    # Variance from local reconstruction
    if (!is.null(B_i) && ncol(B_i) > 0) {
      recon_local <- U_i %*% t(B_i)
      var_local[i] <- sum(recon_local^2)
    }
  }

  total <- var_global + var_local
  pct_global <- 100 * var_global / total
  pct_local <- 100 * var_local / total

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    mat <- rbind(Global = pct_global, Local = pct_local)
    colnames(mat) <- block_names
    barplot(mat, beside = FALSE, col = c("steelblue", "coral"),
            main = "BaMFA: Global vs Local Variance",
            ylab = "% Variance", las = 2)
    legend("topright", legend = c("Global", "Local"),
           fill = c("steelblue", "coral"), bty = "n")
    return(invisible(data.frame(block = block_names, global = pct_global, local = pct_local)))
  }

  df <- data.frame(
    block = rep(block_names, 2),
    type = rep(c("Global", "Local"), each = n_blocks),
    pct = c(pct_global, pct_local)
  )
  df$block <- factor(df$block, levels = block_names)
  df$type <- factor(df$type, levels = c("Local", "Global"))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$block, y = .data$pct, fill = .data$type)) +
    ggplot2::geom_col(position = "stack", alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c("Global" = "steelblue", "Local" = "coral")) +
    ggplot2::labs(
      x = "Block",
      y = "% Variance",
      title = "BaMFA: Global vs Local Variance Decomposition",
      fill = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p
}

# =============================================================================
# BADA plotting methods
# =============================================================================

#' Plot BaDA Results
#'
#' Displays scores from a Barycentric Discriminant Analysis, colored by class.
#'
#' @param x A fitted `bada` object.
#' @param dims Integer vector of length 2 specifying which components to plot.
#' @param labels Logical; if `TRUE`, label points.
#' @param show_barycenters Logical; if `TRUE`, show class centroids.
#' @param col Point colors. If `NULL`, colors by class label.
#' @param pch Point character.
#' @param cex Point size.
#' @param main Plot title.
#' @param ... Additional arguments passed to [plot()].
#'
#' @return Invisibly returns the plotted scores.
#' @export
plot.bada <- function(x, dims = c(1, 2), labels = FALSE, show_barycenters = TRUE,
                      col = NULL, pch = 19, cex = 1, main = "BaDA Score Plot", ...) {
  S <- x$s
  d1 <- dims[1]
  d2 <- dims[2]

  class_labels <- x$labels
  n_classes <- length(x$label_set)

  if (is.null(col)) {
    pal <- hcl.colors(n_classes, palette = "Set2")
    col <- pal[as.integer(class_labels)]
  }

  var_exp <- if (!is.null(x$sdev)) {
    round(100 * x$sdev^2 / sum(x$sdev^2), 1)
  } else {
    rep(NA, ncol(S))
  }

  xlab <- if (!is.na(var_exp[d1])) {
    sprintf("Component %d (%.1f%%)", d1, var_exp[d1])
  } else {
    sprintf("Component %d", d1)
  }
  ylab <- if (!is.na(var_exp[d2])) {
    sprintf("Component %d (%.1f%%)", d2, var_exp[d2])
  } else {
    sprintf("Component %d", d2)
  }

  plot(S[, d1], S[, d2], col = col, pch = pch, cex = cex,
       xlab = xlab, ylab = ylab, main = main, ...)
  abline(h = 0, v = 0, col = "gray70", lty = 2)

  if (show_barycenters) {
    fscores <- x$fscores
    pal <- hcl.colors(n_classes, palette = "Set2")
    points(fscores[, d1], fscores[, d2], pch = 17, cex = cex * 1.5, col = pal)
    text(fscores[, d1], fscores[, d2], labels = x$label_set, pos = 3, cex = cex * 0.8)
  }

  if (labels) {
    text(S[, d1], S[, d2], labels = seq_len(nrow(S)), pos = 3, cex = cex * 0.7)
  }

  legend("topright", legend = x$label_set, fill = hcl.colors(n_classes, palette = "Set2"),
         title = "Class", bty = "n", cex = 0.8)

  invisible(S[, dims, drop = FALSE])
}

#' Autoplot for BaDA
#'
#' Creates a ggplot2 visualization of BaDA scores, colored by class.
#'
#' @param object A fitted `bada` object.
#' @param dims Integer vector of length 2 specifying components to plot.
#' @param labels Logical; if `TRUE`, add text labels.
#' @param show_barycenters Logical; if `TRUE`, show class centroids.
#' @param alpha Point transparency.
#' @param size Point size.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object.
#' @export
autoplot.bada <- function(object, dims = c(1, 2), labels = FALSE,
                          show_barycenters = TRUE, alpha = 0.7, size = 2, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot.", call. = FALSE)
  }

  S <- object$s
  d1 <- dims[1]
  d2 <- dims[2]

  var_exp <- if (!is.null(object$sdev)) {
    round(100 * object$sdev^2 / sum(object$sdev^2), 1)
  } else {
    rep(NA, ncol(S))
  }

  xlab <- if (!is.na(var_exp[d1])) {
    sprintf("Component %d (%.1f%%)", d1, var_exp[d1])
  } else {
    sprintf("Component %d", d1)
  }
  ylab <- if (!is.na(var_exp[d2])) {
    sprintf("Component %d (%.1f%%)", d2, var_exp[d2])
  } else {
    sprintf("Component %d", d2)
  }

  df <- data.frame(
    x = S[, d1],
    y = S[, d2],
    class = object$labels
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$class)) +
    ggplot2::geom_point(alpha = alpha, size = size) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60")

  if (show_barycenters) {
    fscores <- object$fscores
    df_centers <- data.frame(
      x = fscores[, d1],
      y = fscores[, d2],
      class = object$label_set
    )
    p <- p +
      ggplot2::geom_point(data = df_centers, size = size * 2, shape = 17) +
      ggplot2::geom_text(data = df_centers, ggplot2::aes(label = .data$class),
                          vjust = -1, size = 3)
  }

  if (labels) {
    df$label <- seq_len(nrow(df))
    p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$label), vjust = -0.5, size = 2.5)
  }

  p <- p +
    ggplot2::labs(x = xlab, y = ylab, title = "BaDA Score Plot", color = "Class") +
    ggplot2::theme_minimal()

  p
}

#' @rdname plot_variance
#' @export
plot_variance.bada <- function(x, type = c("bar", "line"), ...) {
  type <- match.arg(type)
  sdev <- x$sdev
  if (is.null(sdev) || length(sdev) == 0) {
    stop("No sdev found in the fitted object.", call. = FALSE)
  }

  var_exp <- sdev^2
  pct <- 100 * var_exp / sum(var_exp)
  cum_pct <- cumsum(pct)
  ncomp <- length(sdev)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    if (type == "bar") {
      barplot(pct, names.arg = seq_len(ncomp),
              main = "BaDA Variance Explained", xlab = "Component",
              ylab = "% Variance", col = "steelblue")
    } else {
      plot(seq_len(ncomp), pct, type = "b", pch = 19,
           main = "BaDA Scree Plot", xlab = "Component", ylab = "% Variance")
      lines(seq_len(ncomp), cum_pct, type = "b", pch = 17, col = "red")
    }
    return(invisible(data.frame(component = seq_len(ncomp), pct = pct, cumulative = cum_pct)))
  }

  df <- data.frame(component = seq_len(ncomp), pct = pct, cumulative = cum_pct)

  if (type == "bar") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(.data$component), y = .data$pct)) +
      ggplot2::geom_col(fill = "steelblue", alpha = 0.8) +
      ggplot2::labs(x = "Component", y = "% Variance", title = "BaDA Variance Explained") +
      ggplot2::theme_minimal()
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$component)) +
      ggplot2::geom_line(ggplot2::aes(y = .data$pct), color = "steelblue", linewidth = 1) +
      ggplot2::geom_point(ggplot2::aes(y = .data$pct), color = "steelblue", size = 3) +
      ggplot2::geom_line(ggplot2::aes(y = .data$cumulative), color = "firebrick",
                          linewidth = 1, linetype = "dashed") +
      ggplot2::geom_point(ggplot2::aes(y = .data$cumulative), color = "firebrick", size = 2, shape = 17) +
      ggplot2::scale_x_continuous(breaks = seq_len(ncomp)) +
      ggplot2::labs(x = "Component", y = "% Variance", title = "BaDA Scree Plot") +
      ggplot2::theme_minimal()
  }

  p
}

#' Plot Class Barycenters for BaDA
#'
#' Visualizes class centroids with optional confidence ellipses.
#'
#' @param x A fitted `bada` object.
#' @param dims Integer vector of length 2 for component axes.
#' @param ellipses Logical; if `TRUE`, draw confidence ellipses around classes.
#' @param level Confidence level for ellipses (default 0.95).
#' @param ... Additional arguments.
#'
#' @return A ggplot object or base R plot.
#' @export
plot_barycenters <- function(x, dims = c(1, 2), ellipses = TRUE, level = 0.95, ...) {
  if (!inherits(x, "bada")) {
    stop("plot_barycenters is only available for bada objects.", call. = FALSE)
  }

  S <- x$s
  d1 <- dims[1]
  d2 <- dims[2]
  fscores <- x$fscores
  class_labels <- x$labels
  label_set <- x$label_set

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    n_classes <- length(label_set)
    pal <- hcl.colors(n_classes, palette = "Set2")

    plot(S[, d1], S[, d2], col = pal[as.integer(class_labels)],
         pch = 1, cex = 0.5, xlab = sprintf("Component %d", d1),
         ylab = sprintf("Component %d", d2), main = "BaDA Class Barycenters")
    points(fscores[, d1], fscores[, d2], pch = 17, cex = 2, col = pal)
    text(fscores[, d1], fscores[, d2], labels = label_set, pos = 3)
    abline(h = 0, v = 0, col = "gray70", lty = 2)
    legend("topright", legend = label_set, fill = pal, bty = "n")
    return(invisible(fscores[, dims]))
  }

  df <- data.frame(
    x = S[, d1],
    y = S[, d2],
    class = class_labels
  )

  df_centers <- data.frame(
    x = fscores[, d1],
    y = fscores[, d2],
    class = label_set
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$class))

  if (ellipses) {
    p <- p + ggplot2::stat_ellipse(level = level, linewidth = 0.8)
  }

  p <- p +
    ggplot2::geom_point(alpha = 0.3, size = 1.5) +
    ggplot2::geom_point(data = df_centers, size = 4, shape = 17) +
    ggplot2::geom_text(data = df_centers, ggplot2::aes(label = .data$class),
                        vjust = -1.2, fontface = "bold", show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      x = sprintf("Component %d", d1),
      y = sprintf("Component %d", d2),
      title = "BaDA Class Barycenters",
      color = "Class"
    ) +
    ggplot2::theme_minimal()

  p
}

#' @export
plot_partial_scores.bada <- function(x, dims = c(1, 2), blocks = NULL,
                                     connect = TRUE, show_consensus = TRUE,
                                     alpha = 0.6, ...) {
  d1 <- dims[1]
  d2 <- dims[2]

  # For bada, we need to extract block-specific projections
  # The consensus is x$s (all observations)
  S_consensus <- x$fscores  # Class barycenters as consensus
  n <- nrow(S_consensus)

  # Subjects/blocks
  subject_set <- names(x$block_indices)
  if (is.null(blocks)) blocks <- seq_along(x$block_indices)
  if (is.null(subject_set)) subject_set <- paste0("Subject", seq_along(x$block_indices))
  block_names <- subject_set[blocks]

  # Compute per-subject barycenters
  partial_list <- lapply(blocks, function(b) {
    # Get observations for this subject
    subject_mask <- x$subjects == subject_set[b]
    S_subj <- x$s[subject_mask, , drop = FALSE]
    labels_subj <- x$labels[subject_mask]
    # Compute barycenters per class for this subject
    bary <- matrix(NA, nrow = n, ncol = ncol(S_subj))
    for (i in seq_along(x$label_set)) {
      mask <- labels_subj == x$label_set[i]
      if (sum(mask) > 0) {
        bary[i, ] <- colMeans(S_subj[mask, , drop = FALSE])
      } else {
        bary[i, ] <- S_consensus[i, ]  # fallback
      }
    }
    bary
  })

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(.plot_partial_base(S_consensus, partial_list, dims, blocks, block_names,
                              connect, show_consensus, alpha, ...))
  }

  df_list <- lapply(seq_along(blocks), function(i) {
    ps <- partial_list[[i]]
    data.frame(
      x = ps[, d1],
      y = ps[, d2],
      obs = x$label_set,
      block = block_names[i]
    )
  })
  df_partial <- do.call(rbind, df_list)

  df_global <- data.frame(
    x = S_consensus[, d1],
    y = S_consensus[, d2],
    obs = x$label_set
  )

  p <- ggplot2::ggplot()

  if (connect) {
    for (i in seq_along(blocks)) {
      df_lines <- data.frame(
        x = partial_list[[i]][, d1],
        y = partial_list[[i]][, d2],
        xend = S_consensus[, d1],
        yend = S_consensus[, d2],
        block = block_names[i]
      )
      p <- p + ggplot2::geom_segment(
        data = df_lines,
        ggplot2::aes(x = .data$x, y = .data$y,
                     xend = .data$xend, yend = .data$yend,
                     color = .data$block),
        alpha = alpha * 0.5, linewidth = 0.3
      )
    }
  }

  p <- p + ggplot2::geom_point(
    data = df_partial,
    ggplot2::aes(x = .data$x, y = .data$y, color = .data$block),
    alpha = alpha, size = 2
  )

  if (show_consensus) {
    p <- p + ggplot2::geom_point(
      data = df_global,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "black", size = 3, shape = 17
    ) +
      ggplot2::geom_text(
        data = df_global,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$obs),
        vjust = -1, size = 3
      )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      x = sprintf("Component %d", d1),
      y = sprintf("Component %d", d2),
      title = "BaDA Subject Partial Scores",
      subtitle = "Triangles = group barycenters",
      color = "Subject"
    ) +
    ggplot2::theme_minimal()

  p
}

# =============================================================================
# Penalized MFA plotting methods
# =============================================================================

#' Plot Penalized MFA Results
#'
#' Displays consensus scores from a Penalized MFA analysis.
#'
#' @param x A fitted `penalized_mfa` object.
#' @param dims Integer vector of length 2 specifying which components to plot.
#' @param block Which block's scores to show, or `"consensus"` for averaged.
#' @param labels Logical; if `TRUE`, label points.
#' @param col Point colors.
#' @param pch Point character.
#' @param cex Point size.
#' @param main Plot title.
#' @param ... Additional arguments passed to [plot()].
#'
#' @return Invisibly returns the plotted scores.
#' @export
plot.penalized_mfa <- function(x, dims = c(1, 2), block = "consensus",
                               labels = FALSE, col = NULL, pch = 19, cex = 1,
                               main = "Penalized MFA Scores", ...) {

  S_list <- x$scores_list
  if (is.null(S_list)) S_list <- attr(x, "scores_list")
  if (is.null(S_list)) {
    stop(
      "Scores are not available for this penalized_mfa object. ",
      "Provide `data=` to compute them (see ?autoplot.penalized_mfa) or refit with a newer muscal version.",
      call. = FALSE
    )
  }

  S <- if (identical(block, "consensus")) {
    Reduce(`+`, S_list) / length(S_list)
  } else {
    S_list[[block]]
  }

  if (is.null(col)) col <- "steelblue"

  .plot_scores_base(S, dims = dims, labels = labels, col = col,
                    pch = pch, cex = cex, main = main,
                    sdev = x$sdev, ...)
}

#' Autoplot for Penalized MFA
#'
#' Creates a ggplot2 visualization of Penalized MFA scores.
#'
#' @param object A fitted `penalized_mfa` object.
#' @param dims Integer vector of length 2 specifying components to plot.
#' @param block Which block's scores to show, or `"consensus"` for averaged.
#' @param labels Logical; if `TRUE`, add text labels.
#' @param color Optional variable for point coloring.
#' @param alpha Point transparency.
#' @param size Point size.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object.
#' @export
autoplot.penalized_mfa <- function(object, dims = c(1, 2), block = "consensus",
                                   labels = FALSE, color = NULL, alpha = 0.8, size = 2, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot.", call. = FALSE)
  }

  data <- list(...)$data

  S_list <- object$scores_list
  if (is.null(S_list)) S_list <- attr(object, "scores_list")

  if (is.null(S_list)) {
    if (is.null(data)) {
      stop(
        "Scores are not stored for penalized_mfa. ",
        "Pass `data=` (the original blocks) to autoplot() to compute them.",
        call. = FALSE
      )
    }

    if (inherits(data, "multiblock")) data <- as.list(data)

    if (!is.list(data) || length(data) == 0) {
      stop("`data` must be the original list/multiblock of data matrices.", call. = FALSE)
    }

    # Recompute per-block scores as Xp_k %*% V_k (projected after preprocessing)
    V_list <- attr(object, "V_list")
    if (is.null(V_list)) stop("Missing V_list in penalized_mfa object.", call. = FALSE)

    # Use stored block-wise preprocessors when available
    preproc_list <- object$preproc$proclist %||% NULL
    if (is.null(preproc_list)) {
      # Fallback: apply the concat preprocessor then slice by block indices.
      Xp_concat <- multivarious::transform(object$preproc, do.call(cbind, data))
      Xp_split <- lapply(seq_along(object$block_indices), function(i) {
        Xp_concat[, object$block_indices[[i]], drop = FALSE]
      })
      names(Xp_split) <- names(object$block_indices)
    } else {
      Xp_split <- lapply(seq_along(data), function(i) {
        multivarious::transform(preproc_list[[i]], data[[i]])
      })
      names(Xp_split) <- names(data)
    }

    S_list <- lapply(seq_along(V_list), function(i) {
      Xp_split[[i]] %*% V_list[[i]]
    })
    names(S_list) <- names(V_list)
    attr(object, "scores_list") <- S_list
  }

  if (identical(block, "consensus")) {
    S <- Reduce(`+`, S_list) / length(S_list)
    title <- "Penalized MFA Consensus Scores"
  } else {
    S <- S_list[[block]]
    title <- sprintf("Penalized MFA Scores (Block %s)", block)
  }

  .autoplot_scores(S, dims = dims, labels = labels, color = color,
                   alpha = alpha, size = size, sdev = object$sdev,
                   title = title)
}

#' @rdname plot_convergence
#' @export
plot_convergence.penalized_mfa <- function(x, log_scale = FALSE, ...) {
  obj <- x$objective_trace
  if (is.null(obj)) obj <- attr(x, "obj_values")
  if (is.null(obj) || length(obj) == 0) {
    stop("No objective_trace found in the fitted object.", call. = FALSE)
  }

  n_iter <- length(obj)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    if (log_scale) {
      plot(seq_len(n_iter), log(obj), type = "b", pch = 19,
           xlab = "Iteration", ylab = "log(Objective)",
           main = "Penalized MFA Convergence")
    } else {
      plot(seq_len(n_iter), obj, type = "b", pch = 19,
           xlab = "Iteration", ylab = "Objective",
           main = "Penalized MFA Convergence")
    }
    return(invisible(obj))
  }

  df <- data.frame(iteration = seq_len(n_iter), objective = obj)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$objective)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(color = "steelblue", size = 2) +
    ggplot2::labs(
      x = "Iteration",
      y = if (log_scale) "log(Objective)" else "Objective",
      title = "Penalized MFA Convergence"
    ) +
    ggplot2::theme_minimal()

  if (log_scale) {
    p <- p + ggplot2::scale_y_log10()
  }

  p
}

#' @rdname plot_variance
#' @export
plot_variance.penalized_mfa <- function(x, type = c("bar", "line"), ...) {
  type <- match.arg(type)
  sdev <- x$sdev
  if (is.null(sdev) || length(sdev) == 0) {
    stop("No sdev found in the fitted object.", call. = FALSE)
  }

  var_exp <- sdev^2
  pct <- 100 * var_exp / sum(var_exp)
  cum_pct <- cumsum(pct)
  ncomp <- length(sdev)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    if (type == "bar") {
      barplot(pct, names.arg = seq_len(ncomp),
              main = "Penalized MFA Variance Explained", xlab = "Component",
              ylab = "% Variance", col = "steelblue")
    } else {
      plot(seq_len(ncomp), pct, type = "b", pch = 19,
           main = "Penalized MFA Scree Plot", xlab = "Component", ylab = "% Variance")
      lines(seq_len(ncomp), cum_pct, type = "b", pch = 17, col = "red")
    }
    return(invisible(data.frame(component = seq_len(ncomp), pct = pct, cumulative = cum_pct)))
  }

  df <- data.frame(component = seq_len(ncomp), pct = pct, cumulative = cum_pct)

  if (type == "bar") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(.data$component), y = .data$pct)) +
      ggplot2::geom_col(fill = "steelblue", alpha = 0.8) +
      ggplot2::labs(x = "Component", y = "% Variance", title = "Penalized MFA Variance Explained") +
      ggplot2::theme_minimal()
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$component)) +
      ggplot2::geom_line(ggplot2::aes(y = .data$pct), color = "steelblue", linewidth = 1) +
      ggplot2::geom_point(ggplot2::aes(y = .data$pct), color = "steelblue", size = 3) +
      ggplot2::geom_line(ggplot2::aes(y = .data$cumulative), color = "firebrick",
                          linewidth = 1, linetype = "dashed") +
      ggplot2::geom_point(ggplot2::aes(y = .data$cumulative), color = "firebrick", size = 2, shape = 17) +
      ggplot2::scale_x_continuous(breaks = seq_len(ncomp)) +
      ggplot2::labs(x = "Component", y = "% Variance", title = "Penalized MFA Scree Plot") +
      ggplot2::theme_minimal()
  }

  p
}

#' @export
plot_partial_scores.penalized_mfa <- function(x, dims = c(1, 2), blocks = NULL,
                                              connect = TRUE, show_consensus = TRUE,
                                              alpha = 0.6, ...) {
  d1 <- dims[1]
  d2 <- dims[2]

  # Consensus: mean of block scores
  S_list <- x$scores_list
  if (is.null(S_list)) S_list <- attr(x, "scores_list")
  if (is.null(S_list)) {
    stop(
      "Scores are not available for this penalized_mfa object. ",
      "Pass `data=` to autoplot.penalized_mfa first (or refit with a newer muscal version).",
      call. = FALSE
    )
  }

  S_consensus <- Reduce(`+`, S_list) / length(S_list)
  n <- nrow(S_consensus)

  if (is.null(blocks)) blocks <- seq_along(S_list)
  block_names <- names(S_list)[blocks]
  if (is.null(block_names)) block_names <- paste0("Block", blocks)

  partial_list <- S_list[blocks]

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(.plot_partial_base(S_consensus, partial_list, dims, blocks, block_names,
                              connect, show_consensus, alpha, ...))
  }

  df_list <- lapply(seq_along(blocks), function(i) {
    ps <- partial_list[[i]]
    data.frame(
      x = ps[, d1],
      y = ps[, d2],
      obs = seq_len(n),
      block = block_names[i]
    )
  })
  df_partial <- do.call(rbind, df_list)

  df_global <- data.frame(
    x = S_consensus[, d1],
    y = S_consensus[, d2],
    obs = seq_len(n)
  )

  p <- ggplot2::ggplot()

  if (connect) {
    for (i in seq_along(blocks)) {
      df_lines <- data.frame(
        x = partial_list[[i]][, d1],
        y = partial_list[[i]][, d2],
        xend = S_consensus[, d1],
        yend = S_consensus[, d2],
        block = block_names[i]
      )
      p <- p + ggplot2::geom_segment(
        data = df_lines,
        ggplot2::aes(x = .data$x, y = .data$y,
                     xend = .data$xend, yend = .data$yend,
                     color = .data$block),
        alpha = alpha * 0.5, linewidth = 0.3
      )
    }
  }

  p <- p + ggplot2::geom_point(
    data = df_partial,
    ggplot2::aes(x = .data$x, y = .data$y, color = .data$block),
    alpha = alpha, size = 2
  )

  if (show_consensus) {
    p <- p + ggplot2::geom_point(
      data = df_global,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "black", size = 3, shape = 1
    )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      x = sprintf("Component %d", d1),
      y = sprintf("Component %d", d2),
      title = "Penalized MFA Block Scores",
      subtitle = "Open circles = consensus (mean)",
      color = "Block"
    ) +
    ggplot2::theme_minimal()

  p
}

#' Plot Loading Alignment for Penalized MFA
#'
#' Visualizes how well block-specific loadings align with each other,
#' showing the effect of the penalty parameter.
#'
#' @param x A fitted `penalized_mfa` object.
#' @param dims Integer vector of length 2 for component axes.
#' @param top_n Number of top variables to show per block.
#' @param ... Additional arguments.
#'
#' @return A ggplot object or base R plot.
#' @export
plot_loading_alignment <- function(x, dims = c(1, 2), top_n = 10, ...) {
  if (!inherits(x, "penalized_mfa")) {
    stop("plot_loading_alignment is only available for penalized_mfa objects.", call. = FALSE)
  }

  V_list <- x$V_list
  if (is.null(V_list)) V_list <- attr(x, "V_list")
  if (is.null(V_list)) stop("Missing V_list in penalized_mfa object.", call. = FALSE)
  n_blocks <- length(V_list)
  block_names <- names(V_list)
  if (is.null(block_names)) block_names <- paste0("Block", seq_len(n_blocks))

  d1 <- dims[1]
  d2 <- dims[2]

  p_dims <- vapply(V_list, nrow, integer(1))
  if (!all(p_dims == p_dims[1])) {
    stop(
      "plot_loading_alignment() requires blocks with the same feature dimension (same number of rows in each V_k).",
      call. = FALSE
    )
  }

  # Compute mean loadings
  V_mean <- Reduce("+", V_list) / n_blocks

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R: just plot the mean loadings
    norms <- sqrt(rowSums(V_mean[, dims, drop = FALSE]^2))
    top_idx <- order(norms, decreasing = TRUE)[seq_len(min(top_n, nrow(V_mean)))]

    plot(V_mean[top_idx, d1], V_mean[top_idx, d2], pch = 19,
         xlab = sprintf("Component %d", d1),
         ylab = sprintf("Component %d", d2),
         main = "Penalized MFA Mean Loadings")
    abline(h = 0, v = 0, col = "gray70", lty = 2)
    return(invisible(V_mean))
  }

  # Build data frame with all block loadings
  df_list <- list()
  for (i in seq_len(n_blocks)) {
    V_i <- V_list[[i]]
    norms <- sqrt(rowSums(V_i[, dims, drop = FALSE]^2))
    top_idx <- order(norms, decreasing = TRUE)[seq_len(min(top_n, nrow(V_i)))]
    df_list[[i]] <- data.frame(
      x = V_i[top_idx, d1],
      y = V_i[top_idx, d2],
      var_idx = top_idx,
      block = block_names[i]
    )
  }
  df <- do.call(rbind, df_list)

  # Mean loadings
  norms_mean <- sqrt(rowSums(V_mean[, dims, drop = FALSE]^2))
  top_idx_mean <- order(norms_mean, decreasing = TRUE)[seq_len(min(top_n, nrow(V_mean)))]
  df_mean <- data.frame(
    x = V_mean[top_idx_mean, d1],
    y = V_mean[top_idx_mean, d2],
    var_idx = top_idx_mean
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = df,
                         ggplot2::aes(x = .data$x, y = .data$y, color = .data$block),
                         alpha = 0.5, size = 2) +
    ggplot2::geom_point(data = df_mean,
                         ggplot2::aes(x = .data$x, y = .data$y),
                         color = "black", size = 3, shape = 4, stroke = 1.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      x = sprintf("Component %d", d1),
      y = sprintf("Component %d", d2),
      title = "Penalized MFA Loading Alignment",
      subtitle = sprintf("Top %d variables per block; crosses = mean", top_n),
      color = "Block"
    ) +
    ggplot2::theme_minimal()

  p
}

#' @export
plot_loadings.penalized_mfa <- function(x, dims = c(1, 2), type = c("circle", "bar"),
                                        component = 1, color_by = "block", top_n = 20, ...) {
  type <- match.arg(type)

  V <- if (type == "circle") {
    .circle_loadings_matrix(x, fallback = x$v)
  } else {
    as.matrix(x$v)
  }

  # Block labels from concatenated indices
  block_labels <- character(nrow(V))
  block_names <- names(x$block_indices)
  if (is.null(block_names)) block_names <- paste0("Block", seq_along(x$block_indices))
  for (i in seq_along(x$block_indices)) {
    block_labels[x$block_indices[[i]]] <- block_names[i]
  }

  var_names <- rownames(V)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(nrow(V)))

  .plot_loadings_impl(V, dims, type, component, color_by = color_by, top_n,
                      block_labels, var_names, title_prefix = "Penalized MFA", ...)
}
