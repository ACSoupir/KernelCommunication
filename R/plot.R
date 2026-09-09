#' Plot Multi-Scale Moran's I Results
#'
#' Creates ggplot2 visualization of Moran's I across bandwidths with
#' permutation-based confidence bands.
#'
#' The shaded band is the \emph{null envelope}: quantiles of the permutation
#' distribution, showing where the statistic would fall under complete spatial
#' randomness. It is not a confidence interval for the observed curve. Its width is
#' set by \code{envelope_level} in \code{\link{compute_significance}}, not here.
#'
#' @param result Object of class \code{kernelcomm_result} (must have significance computed)
#' @param type Character. One of:
#'   \itemize{
#'     \item "observed": observed Moran's I with the null envelope (default)
#'     \item "difference": observed minus the permutation mean
#'     \item "zscore": deviation standardised by the permutation standard deviation
#'   }
#' @param show_permutations Logical. Draw individual permutation curves behind the
#'   observed one (default FALSE). This is the most direct way to judge whether the
#'   observed curve is unusual.
#' @param show_permutations_n Integer, how many permutation curves to draw when
#'   \code{show_permutations = TRUE} (default 100). Drawing all of them is slow and
#'   illegible for large \code{perms}.
#' @param ribbon_alpha Numeric transparency of the envelope (default 0.1). Named to
#'   avoid confusion with the significance threshold, which this is not.
#' @param point_size Size of points (default: 2)
#' @param line_width Width of lines (default: 1)
#' @param colors Named list with colors: list(line = "...", ribbon = "...", null = "...")
#'
#' @return ggplot2::ggplot object
#' @export
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point labs theme_minimal theme element_blank
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' coords <- matrix(runif(2 * n), ncol = 2)
#' g1 <- rnorm(n)
#' sm <- exp(-as.matrix(dist(coords))^2 / (2 * 0.1^2))
#' g2 <- as.vector(sm %*% g1) / rowSums(sm) + rnorm(n, sd = 0.5)
#'
#' res <- multiScaleAC(coords, g1, g2, h = c(0.05, 0.1, 0.2), perms = 199)
#' res <- compute_significance(res)
#'
#' plot_scale_curve(res)
#' plot_scale_curve(res, type = "difference")
#' plot_scale_curve(res, show_permutations = TRUE, show_permutations_n = 50)
plot_scale_curve <- function(result,
                             type = "observed",
                             show_permutations = FALSE,
                             show_permutations_n = 100,
                             ribbon_alpha = 0.1,
                             point_size = 2,
                             line_width = 1,
                             colors = NULL) {
  
  # Validate input
  if (!inherits(result, "kernelcomm_result")) {
    stop("result must be a kernelcomm_result object")
  }
  
  if (length(result@p_values) == 0) {
    stop("Must call compute_significance() before plotting")
  }
  
  type <- match.arg(type, c("observed", "difference", "zscore"))
  
  # Default colors
  default_colors <- list(
    line = "#2E86AB",
    ribbon = "#A8D5E2",
    null = "#E63946",
    point = "#2E86AB"
  )
  
  if (is.null(colors)) {
    colors <- default_colors
  } else {
    colors <- utils::modifyList(default_colors, colors)
  }
  
  # Prepare data
  h <- result@h
  n_perms <- ncol(result@permutations)
  
  if (type == "observed") {
    y_values <- result@observed
    y_label <- "Moran's I"
    ribbon_min <- result@null_lower
    ribbon_max <- result@null_upper
    null_line <- result@perm_mean
    
  } else if (type == "difference") {
    y_values <- result@observed - result@perm_mean
    y_label <- "Observed - Mean(Permuted)"
    ribbon_min <- result@null_lower - result@perm_mean
    ribbon_max <- result@null_upper - result@perm_mean
    null_line <- rep(0, length(h))
    
  } else if (type == "zscore") {
    y_values <- result@z_scores
    y_label <- "Standardised deviation from the permutation mean"
    # The envelope comes from the permutation distribution itself, standardised the same way as the
    # observed value -- not from +-1.96, which is what this drew before.
    #
    # A normal band is wrong here. The permutation null of this statistic is skewed whenever the
    # expression data is skewed, which for zero-inflated single-cell data is the usual case, so a
    # symmetric +-1.96 band misstates both tails and would disagree with the p-values shown
    # elsewhere in the same object. See compute_significance() for the measured size of that error.
    sd_safe <- ifelse(is.finite(result@perm_sd) & result@perm_sd > 0, result@perm_sd, NA_real_)
    zperm <- sweep(sweep(result@permutations, 1, result@perm_mean, "-"), 1, sd_safe, "/")
    ribbon_min <- apply(zperm, 1, stats::quantile, probs = 0.025, names = FALSE, na.rm = TRUE)
    ribbon_max <- apply(zperm, 1, stats::quantile, probs = 0.975, names = FALSE, na.rm = TRUE)
    null_line <- rep(0, length(h))
  }
  
  # Create data frame
  df <- data.frame(
    h = h,
    y = y_values,
    null_lower = ribbon_min,
    null_upper = ribbon_max,
    null = null_line
  )
  
  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Please install it.")
  }
  
  # Build plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = h)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = null_lower, ymax = null_upper),
                fill = colors$ribbon, alpha = ribbon_alpha) +
    ggplot2::geom_line(ggplot2::aes(y = y), color = colors$line, linewidth = line_width) +
    ggplot2::geom_point(ggplot2::aes(y = y), color = colors$point, size = point_size) +
    NULL

  # Draw individual permutation curves behind the observed one.
  #
  # This argument was documented and had no implementation at all. Showing the null curves is the most
  # direct way to see whether the observed curve is unusual, so it is worth having rather than
  # removing. Layers are inserted before the existing ones so they sit underneath.
  if (isTRUE(show_permutations) && n_perms > 0) {
    n_show <- min(as.integer(show_permutations_n), n_perms)
    idx <- if (n_perms <= n_show) seq_len(n_perms) else round(seq(1, n_perms, length.out = n_show))
    pdat <- result@permutations[, idx, drop = FALSE]
    if (type == "difference") pdat <- pdat - result@perm_mean
    if (type == "zscore") {
      sd_safe <- ifelse(is.finite(result@perm_sd) & result@perm_sd > 0, result@perm_sd, NA_real_)
      pdat <- sweep(sweep(pdat, 1, result@perm_mean, "-"), 1, sd_safe, "/")
    }
    perm_df <- data.frame(
      h = rep(h, times = ncol(pdat)),
      y = as.vector(pdat),
      perm = rep(seq_len(ncol(pdat)), each = length(h))
    )
    p <- p + ggplot2::geom_line(
      data = perm_df,
      mapping = ggplot2::aes(x = h, y = y, group = perm),
      inherit.aes = FALSE, colour = "grey60", linewidth = 0.25, alpha = 0.5
    )
  }

  p <- p +
    ggplot2::geom_line(ggplot2::aes(y = null), linetype = "dashed", color = colors$null, linewidth = line_width * 0.8) +
    ggplot2::labs(x = "Bandwidth (h)", y = y_label,
         title = paste0("Multi-Scale Moran's I (", result@gene1_name, " x ", result@gene2_name, ")"),
         subtitle = paste0(n_perms, " permutations, ", result@n_cells, " cells")) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none"
    )
  
  # Add significance indicators
  if (type != "zscore") {
    sig_points <- df[result@significant, ]
    if (nrow(sig_points) > 0) {
      p <- p + ggplot2::geom_point(data = sig_points, ggplot2::aes(x = h, y = y),
                         color = colors$null, size = point_size * 1.5, shape = 1)
    }
  }
  
  # No +-1.96 reference lines on the z-score plot.
  #
  # They were drawn here, and they invite exactly the mistake this package warns against: reading a
  # normal-theory threshold off a permutation null that is skewed whenever the expression data is. On
  # real data at small bandwidths the true 2.5% and 97.5% points of this null are nowhere near
  # +-1.96, so the lines would have implied significance boundaries that the accompanying p-values
  # contradict. The shaded envelope, taken from the permutations themselves, carries the same
  # information honestly.

  return(p)
}


#' Plot Permutation Distribution
#'
#' Shows histogram of permutation values with observed statistic marked.
#'
#' @param result kernelcomm_result object
#' @param h_index Which bandwidth to show (default: 1)
#' @param bins Number of histogram bins (default: 30)
#'
#' @return ggplot2::ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' result <- multiScaleAC(coords, gene1, gene2, h = c(10, 20, 40))
#' result <- compute_significance(result)
#' plot_permutation_distribution(result, h_index = 1)
#' }
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline labs theme_minimal
plot_permutation_distribution <- function(result, h_index = 1, bins = 30) {
  
  if (!inherits(result, "kernelcomm_result")) {
    stop("result must be a kernelcomm_result object")
  }
  
  if (length(result@p_values) == 0) {
    stop("Must call compute_significance() first")
  }
  
  if (h_index < 1 || h_index > length(result@h)) {
    stop("h_index out of range")
  }
  
  # Check ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting")
  }
  
  # Extract permutation values for this bandwidth
  perm_vals <- result@permutations[h_index, ]
  observed <- result@observed[h_index]
  h_value <- result@h[h_index]
  p_val <- result@p_values[h_index]
  
  # Create data frame
  df <- data.frame(perm = perm_vals)
  
  # Build plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = perm)) +
    ggplot2::geom_histogram(bins = bins, fill = "#A8D5E2", color = "#2E86AB", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = observed, color = "#E63946", linewidth = 1.2, linetype = "solid") +
    ggplot2::geom_vline(xintercept = result@perm_mean[h_index], color = "#2E86AB", linewidth = 1, linetype = "dashed") +
    ggplot2::labs(x = "Permuted Moran's I", y = "Frequency",
         title = paste0("Permutation Distribution (h = ", h_value, ")"),
         subtitle = paste0("Observed = ", format(observed, digits = 4),
                          ", p = ", format(p_val, digits = 3))) +
    ggplot2::theme_minimal(base_size = 12)
  
  return(p)
}


#' Plot Batch Comparison Heatmap
#'
#' Creates a heatmap comparing significance across multiple samples and bandwidths.
#' Useful for visualizing batch processing results.
#'
#' @param results List of \code{kernelcomm_result} objects (from multiScaleAC_batch)
#' @param value_var Character. Which value to plot: "z_score", "p_value", "observed" (default: "z_score")
#' @param cluster_rows Logical. Cluster samples by similarity? (default: TRUE)
#' @param color_low Character. Color for low values (default: blue)
#' @param color_high Character. Color for high values (default: red)
#' @param na_value Character. Color for NA values (default: gray90)
#' @param show_values Logical. Show numeric values in cells? (default: FALSE)
#' @param value_size Numeric. Text size for cell values (default: 3)
#'
#' @return ggplot2::ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Process multiple samples
#' data_list <- list(sample1 = ccRCCExample, sample2 = ccRCCExample)
#' results <- multiScaleAC_batch(data_list, h_vec = c(10, 20, 40), perms = 99)
#' results <- lapply(results, compute_significance)
#'
#' # Plot heatmap of z-scores
#' plot_batch_comparison(results, value_var = "z_score")
#'
#' # Plot heatmap of p-values
#' plot_batch_comparison(results, value_var = "p_value")
#' }
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 labs theme_minimal theme element_text element_blank
plot_batch_comparison <- function(results,
                                   value_var = "z_score",
                                   cluster_rows = TRUE,
                                   color_low = "#2166AC",
                                   color_high = "#B2182B",
                                   na_value = "gray90",
                                   show_values = FALSE,
                                   value_size = 3) {
  
  # Validate input
  if (!is.list(results) || length(results) == 0) {
    stop("results must be a non-empty list")
  }
  
  if (!all(sapply(results, inherits, "kernelcomm_result"))) {
    stop("All elements in results must be kernelcomm_result objects")
  }
  
  # Check significance computed for all
  if (!all(sapply(results, function(r) length(r@p_values) > 0))) {
    stop("All results must have significance computed. Call compute_significance() on each.")
  }
  
  # Validate value_var
  value_var <- match.arg(value_var, c("z_score", "p_value", "observed"))
  
  # Check ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting")
  }
  
  # Labels come from the list names when present, falling back to the sample_id slot.
  #
  # Results built by multiScaleAC() directly all carry the default sample_id "sample1", so taking the
  # slot alone gave duplicated factor levels and the plot failed with "factor level [2] is
  # duplicated". Only multiScaleAC_batch() sets distinct ids, and a user assembling a list by hand
  # will naturally name it instead.
  sample_ids <- .result_labels(results)
  h_values <- results[[1]]@h

  # Check all results have same h values
  if (!all(sapply(results, function(r) all(r@h == h_values)))) {
    stop("All results must have the same bandwidth values")
  }

  # Build data frame
  df_list <- lapply(seq_along(results), function(i) {
    r <- results[[i]]
    data.frame(
      sample_id = sample_ids[i],
      h = r@h,
      value = switch(value_var,
                     z_score = r@z_scores,
                     p_value = r@p_values,
                     observed = r@observed),
      stringsAsFactors = FALSE
    )
  })
  
  df <- do.call(rbind, df_list)
  
  # Order samples by similarity of their curves.
  #
  # This argument is documented as "cluster samples by similarity" and previously did no clustering:
  # it sorted by max absolute z-score, and only when value_var happened to be "z_score". It now
  # performs actual hierarchical clustering on the per-sample curves, for any value_var.
  if (cluster_rows && length(sample_ids) > 2) {
    mat <- matrix(NA_real_, nrow = length(sample_ids), ncol = length(h_values),
                  dimnames = list(sample_ids, as.character(h_values)))
    for (sid in sample_ids) {
      sl <- df[df$sample_id == sid, ]
      mat[sid, as.character(sl$h)] <- sl$value
    }
    # Distances need complete rows; fall back to the original order if the data is too sparse or
    # degenerate to cluster, rather than failing inside hclust.
    ord <- tryCatch({
      keep <- colSums(is.na(mat)) == 0
      if (sum(keep) < 2) stop("insufficient complete bandwidths")
      stats::hclust(stats::dist(mat[, keep, drop = FALSE]))$order
    }, error = function(e) seq_along(sample_ids))
    df$sample_id <- factor(df$sample_id, levels = sample_ids[ord])
  } else {
    df$sample_id <- factor(df$sample_id, levels = sample_ids)
  }
  
  # Format h as factor for discrete axis
  df$h <- factor(df$h, levels = h_values)
  
  # Build plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = h, y = sample_id, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(low = color_low, high = color_high, mid = "white",
                         midpoint = 0, na.value = na_value,
                         name = value_var) +
    ggplot2::labs(x = "Bandwidth (h)", y = "Sample",
         title = paste0("Batch Comparison (", length(results), " samples)"),
         subtitle = paste0("Value: ", value_var, " across ", length(h_values), " bandwidths")) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 9),
      legend.position = "right"
    )
  
  # Mark significant cells with an asterisk.
  #
  # inherit.aes = FALSE is essential and was missing. The main plot maps fill = value, so with
  # inheritance on, ggplot evaluated `value` against sig_df, which has no such column, and
  # plot_batch_comparison() errored on its own default value_var as soon as any cell was significant.
  sig_list <- lapply(seq_along(results), function(i) {
    r <- results[[i]]
    if (length(r@significant) == 0) return(NULL)
    data.frame(sample_id = sample_ids[i], h = as.character(r@h),
               significant = as.logical(r@significant), stringsAsFactors = FALSE)
  })
  sig_df <- do.call(rbind, sig_list[!vapply(sig_list, is.null, logical(1))])
  if (!is.null(sig_df) && any(sig_df$significant)) {
    sig_df <- sig_df[sig_df$significant, , drop = FALSE]
    sig_df$sample_id <- factor(sig_df$sample_id, levels = levels(df$sample_id))
    sig_df$h <- factor(sig_df$h, levels = as.character(h_values))
    p <- p + ggplot2::geom_text(
      data = sig_df,
      mapping = ggplot2::aes(x = h, y = sample_id),
      label = "*", inherit.aes = FALSE,
      color = "white", size = value_size * 1.2, vjust = 0.75
    )
  }
  
  # Add numeric values if requested
  if (show_values) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = format(value, digits = 2)),
                       size = value_size, color = "black")
  }
  
  return(p)
}


#' Plot Significance Summary Across Samples
#'
#' Creates a bar plot showing number of significant bandwidths per sample.
#'
#' @param results List of \code{kernelcomm_result} objects
#' @param alpha Significance threshold (default: 0.05)
#' @param sort_by Character. Sort bars by: "count", "sample_id" (default: "count")
#'
#' @return ggplot2::ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' results <- multiScaleAC_batch(data_list, h_vec = seq(10, 80, 10))
#' results <- lapply(results, compute_significance)
#' plot_significance_summary(results)
#' }
#' @importFrom ggplot2 ggplot aes geom_bar geom_hline labs theme_minimal theme element_text element_blank
plot_significance_summary <- function(results,
                                       alpha = 0.05,
                                       sort_by = "count") {
  
  # Validate input
  if (!is.list(results) || length(results) == 0) {
    stop("results must be a non-empty list")
  }
  
  if (!all(sapply(results, inherits, "kernelcomm_result"))) {
    stop("All elements must be kernelcomm_result objects")
  }
  
  # Check significance computed
  if (!all(sapply(results, function(r) length(r@p_values) > 0))) {
    stop("All results must have significance computed")
  }
  
  # Check ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting")
  }
  
  sort_by <- match.arg(sort_by, c("count", "sample_id"))

  # `alpha` is applied here rather than only appearing in the plot title. Previously the bars came
  # from each result's `significant` slot, computed with whatever alpha was passed to
  # compute_significance(), while the title announced this function's alpha -- so the label could
  # contradict the bars.
  summary_data <- data.frame(
    sample_id = .result_labels(results),
    n_significant = vapply(results, function(r) sum(r@p_values < alpha, na.rm = TRUE), numeric(1)),
    n_total = vapply(results, function(r) length(r@h), numeric(1)),
    stringsAsFactors = FALSE
  )

  # Both orderings are implemented. "sample_id" was documented but had no branch, so it silently
  # fell through to whatever order the list happened to be in.
  summary_data <- switch(
    sort_by,
    count     = summary_data[order(summary_data$n_significant, decreasing = TRUE), ],
    sample_id = summary_data[order(summary_data$sample_id), ]
  )
  summary_data$sample_id <- factor(summary_data$sample_id, levels = summary_data$sample_id)
  
  # Build plot
  p <- ggplot2::ggplot(summary_data, ggplot2::aes(x = sample_id, y = n_significant)) +
    ggplot2::geom_bar(stat = "identity", fill = "#2E86AB") +
    ggplot2::geom_hline(yintercept = summary_data$n_total[1], linetype = "dashed", color = "gray50") +
    ggplot2::labs(x = "Sample", y = "Number of Significant Bandwidths",
         title = paste0("Significance Summary (alpha = ", alpha, ")"),
         subtitle = paste0("Total bandwidths tested: ", summary_data$n_total[1])) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(p)
}
