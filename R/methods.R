#' Print Method for kernelcomm_result
#'
#' @param x kernelcomm_result object
#' @param ... Additional arguments (unused)
#'
#' @return \code{x}, invisibly.
#' @method print kernelcomm_result
#' @export
#'
#' @examples
#' \dontrun{
#' result <- multiScaleAC(coords, gene1, gene2, h = c(0.02, 0.05, 0.1))
#' print(result)
#' }
print.kernelcomm_result <- function(x, ...) {
  cat("multiScaleAC result\n")
  cat("===================\n")
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
  invisible(x)
}


#' Summary Method for kernelcomm_result
#'
#' @param object kernelcomm_result object
#' @param alpha Significance level (default: 0.05)
#' @param ... Additional arguments (unused)
#'
#' @return data.frame with detailed statistics
#' @method summary kernelcomm_result
#' @export
#'
#' @examples
#' \dontrun{
#' result <- multiScaleAC(coords, gene1, gene2, h = c(0.02, 0.05, 0.1))
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
#' @method plot kernelcomm_result
#' @export
#'
#' @examples
#' \dontrun{
#' result <- multiScaleAC(coords, gene1, gene2, h = c(0.02, 0.05, 0.1))
#' result <- compute_significance(result)
#' plot(result)                        # observed curve with the null envelope
#' plot(result, type = "difference")   # deviation from the permutation mean
#' }
plot.kernelcomm_result <- function(x, y, type = "observed", ...) {
  # Check if significance is computed
  if (length(x@p_values) == 0) {
    stop("Must call compute_significance() before plotting")
  }
  
  # Dispatch to plot function
  plot_scale_curve(x, type = type, ...)
}


#' Show Method for kernelcomm_result
#'
#' S4 dispatch for automatic printing. Without this, evaluating a result at the
#' prompt falls through to the default S4 show method, which dumps every slot --
#' including the full n_h x perms permutation matrix and the N x 2 coordinate
#' matrix -- across the console. The S3 \code{print} method only fires on an
#' explicit \code{print(x)} call, so it did not cover the common case.
#'
#' @param object kernelcomm_result object
#' @return \code{object}, invisibly.
#' @keywords internal
setMethod("show", "kernelcomm_result", function(object) {
  print.kernelcomm_result(object)
  invisible(object)
})


#' Helper: NULL coalescing operator
#'
#' Defined for R < 4.4, where base R does not provide it.
#'
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
