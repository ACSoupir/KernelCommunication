#' Print Method for kernelcomm_result
#'
#' @param x kernelcomm_result object
#' @param ... Additional arguments (unused)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- multiScaleAC(coords, gene1, gene2, h = c(10, 20, 40))
#' print(result)
#' }
print.kernelcomm_result <- function(x, ...) {
  cat("KernelCommunication Result\n")
  cat("==========================\n")
  cat("Sample:", x@sample_id, "\n")
  cat("Gene pair:", x@gene1_name, "x", x@gene2_name, "\n")
  cat("Bandwidths:", length(x@h), "values from", 
      format(min(x@h), digits = 3), "to", 
      format(max(x@h), digits = 3), "\n")
  cat("Permutations:", ncol(x@permutations), "\n")
  cat("Cells:", x@n_cells, "\n")
  
  if (length(x@p_values) > 0) {
    cat("Significant bandwidths (alpha=0.05):", sum(x@p_values < 0.05), "/", length(x@h), "\n")
    
    # Show best result
    best_idx <- which.min(x@p_values)
    cat("Most significant: h =", x@h[best_idx], 
        ", p =", format(x@p_values[best_idx], digits = 3),
        ", z =", format(x@z_scores[best_idx], digits = 2), "\n")
  } else {
    cat("Significance: Not computed (use compute_significance())\n")
  }
}


#' Summary Method for kernelcomm_result
#'
#' @param object kernelcomm_result object
#' @param alpha Significance level (default: 0.05)
#' @param ... Additional arguments (unused)
#'
#' @return data.frame with detailed statistics
#' @export
#'
#' @examples
#' \dontrun{
#' result <- multiScaleAC(coords, gene1, gene2, h = c(10, 20, 40))
#' result <- compute_significance(result)
#' summary(result)
#' }
summary.kernelcomm_result <- function(object, alpha = 0.05, ...) {
  if (length(object@p_values) == 0) {
    stop("Must call compute_significance() first")
  }
  
  # Handle p_adjusted - use NA if not computed (length 0)
  p_adj <- if (length(object@p_adjusted) > 0) object@p_adjusted else rep(NA_real_, length(object@h))
  
  # Create summary data.frame
  df <- data.frame(
    h = object@h,
    observed = object@observed,
    perm_mean = object@perm_mean,
    perm_sd = object@perm_sd,
    z_score = object@z_scores,
    p_value = object@p_values,
    p_adjusted = p_adj,
    significant = object@p_values < alpha,
    ci_lower = object@ci_lower,
    ci_upper = object@ci_upper
  )
  
  return(df)
}


#' Plot Method for kernelcomm_result
#'
#' Dispatches to appropriate plot function based on type.
#'
#' @param x kernelcomm_result object
#' @param y Unused
#' @param type Character. One of "observed", "difference", "zscore"
#' @param ... Additional arguments passed to plot_scale_curve
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' result <- multiScaleAC(coords, gene1, gene2, h = c(10, 20, 40))
#' result <- compute_significance(result)
#' plot(result)  # Default: observed with CI
#' plot(result, type = "difference")  # Deviation from null
#' }
plot.kernelcomm_result <- function(x, y, type = "observed", ...) {
  # Check if significance is computed
  if (length(x@p_values) == 0) {
    stop("Must call compute_significance() before plotting")
  }
  
  # Dispatch to plot function
  plot_scale_curve(x, type = type, ...)
}


#' Helper: NULL coalescing operator
#'
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
