#' Compute Significance Statistics from Permutation Results
#'
#' Adds empirical p-values, z-scores, and confidence intervals to
#' a kernelcomm_result object.
#'
#' @param result Object of class \code{kernelcomm_result}
#' @param alternative Character. One of "two.sided", "greater", "less"
#' @param alpha Significance level for CI (default: 0.05)
#'
#' @return Same object with added slots:
#'   \item{p_values}{Empirical p-values}
#'   \item{z_scores}{Standardized effect sizes}
#'   \item{perm_mean}{Mean of permutations}
#'   \item{perm_sd}{SD of permutations}
#'   \item{ci_lower}{Lower CI bound}
#'   \item{ci_upper}{Upper CI bound}
#'   \item{significant}{Logical. p < alpha}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- multiScaleAC(coords, gene1, gene2, h = c(10, 20, 40))
#' result <- compute_significance(result)
#' summary(result)  # Now shows p-values
#' }
compute_significance <- function(result, 
                                  alternative = "two.sided",
                                  alpha = 0.05) {
  
  # Validate input
  if (!inherits(result, "kernelcomm_result")) {
    stop("result must be a kernelcomm_result object")
  }
  
  # Check if significance already computed
  if (!length(result@p_values) == 0) {
    warning("Significance already computed. Recomputing with new parameters.")
  }
  
  # Extract permutation matrix and observed values
  perms <- result@permutations
  observed <- result@observed
  n_perms <- ncol(perms)
  
  # Compute permutation statistics
  perm_mean <- rowMeans(perms)
  perm_sd <- apply(perms, 1, sd)
  
  # Empirical p-values (with correction for finite permutations)
  # p = (r + 1) / (n + 1) where r = rank of observed among permutations
  p_values <- sapply(1:length(observed), function(i) {
    # Count how many permutation values are >= observed
    r <- sum(perms[i, ] >= observed[i])
    (r + 1) / (n_perms + 1)
  })
  
  # Adjust for one-sided tests
  if (alternative == "less") {
    p_values <- 1 - p_values
  } else if (alternative == "greater") {
    # Already correct for "greater" (testing if observed > permuted)
  } else if (alternative != "two.sided") {
    stop("alternative must be 'two.sided', 'greater', or 'less'")
  }
  
  # Z-scores (standardized effect sizes)
  # Handle case where perm_sd is 0 or NA
  z_scores <- rep(NA_real_, length(observed))
  valid_sd <- !is.na(perm_sd) & perm_sd > 0
  z_scores[valid_sd] <- (observed[valid_sd] - perm_mean[valid_sd]) / perm_sd[valid_sd]
  
  # Confidence intervals from permutation distribution
  ci_lower <- apply(perms, 1, quantile, probs = alpha / 2, names = FALSE)
  ci_upper <- apply(perms, 1, quantile, probs = 1 - alpha / 2, names = FALSE)
  
  # Significance indicator
  significant <- p_values < alpha
  
  # Add slots to object
  result@p_values <- p_values
  result@z_scores <- z_scores
  result@perm_mean <- perm_mean
  result@perm_sd <- perm_sd
  result@ci_lower <- ci_lower
  result@ci_upper <- ci_upper
  result@significant <- significant
  
  return(result)
}


#' Adjust P-values for Multiple Testing
#'
#' Applies multiple testing correction across bandwidths.
#'
#' @param result Object of class \code{kernelcomm_result} (must have significance computed)
#' @param method Character. One of "BH", "BY", "bonferroni", "holm", "hochberg", "hommel", "fdr", "none"
#'
#' @return Object with added slot: p_adjusted
#' @export
#'
#' @examples
#' \dontrun{
#' result <- multiScaleAC(coords, gene1, gene2, h = c(10, 20, 40))
#' result <- compute_significance(result)
#' result <- adjust_pvalues(result, method = "BH")
#' summary(result)  # Shows adjusted p-values
#' }
adjust_pvalues <- function(result, method = "BH") {
  
  # Validate input
  if (!inherits(result, "kernelcomm_result")) {
    stop("result must be a kernelcomm_result object")
  }
  
  if (length(result@p_values) == 0) {
    stop("Must call compute_significance() before adjust_pvalues()")
  }
  
  # Validate method
  valid_methods <- c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel", "fdr", "none")
  method <- match.arg(method, valid_methods)
  
  # Adjust p-values
  result@p_adjusted <- stats::p.adjust(result@p_values, method = method)
  
  return(result)
}
