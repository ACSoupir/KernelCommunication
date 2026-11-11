#' Export to mxfda Format
#'
#' Converts a result to a data frame shaped for
#' \code{mxfda::add_summary_function()}, so the interaction curves can go into a
#' functional data analysis, as the companion manuscript does with functional
#' principal component analysis.
#'
#' @section Bandwidth in the radius column:
#' mxfda's spatial summary functions are indexed by a radius \code{r}. This export
#' places the kernel bandwidth \code{h} in that column. \strong{They are not the
#' same quantity.} A radius is a hard cutoff: pairs inside contribute, pairs
#' outside do not. A Gaussian bandwidth is a scale parameter with unbounded
#' support: every pair contributes, weighted by distance, and a bandwidth of
#' \code{h} draws most of its weight from within roughly \code{h} but not all of
#' it.
#'
#' The substitution is reasonable, because both index a family of nested spatial
#' scales and that is all the downstream functional analysis needs. But curves
#' exported this way are not interchangeable with radius-based summary functions
#' such as Ripley's K or a pair correlation function, and should not be pooled with
#' them. This is stated rather than assumed because the previous version of this
#' function made the substitution silently.
#'
#' @section Verification status:
#' The column names here follow mxfda's documented contract for a summary function
#' table. \code{mxfda} could not be installed in the environment where this was
#' developed, so the round trip through
#' \code{mxfda::add_summary_function()} has \strong{not} been executed here. The
#' test suite contains a check that runs it and is skipped when mxfda is absent;
#' run the tests with mxfda installed to confirm.
#'
#' @param result Object of class \code{kernelcomm_result} (must have significance computed)
#' @param sample_id Character. Sample identifier (overrides object slot)
#' @param interaction_name Character. Gene pair name (e.g., "CD274_PDCD1")
#'
#' @return A data frame with one row per bandwidth and columns \code{sample_id},
#'   \code{interaction_name}, \code{r} (the bandwidth), \code{sumfun} (the observed
#'   statistic), \code{csr} (its value under complete spatial randomness),
#'   \code{fundiff} (\code{sumfun - csr}, which is what the manuscript plots),
#'   \code{p_value}, \code{z_score}, \code{null_lower}, \code{null_upper} and
#'   \code{p_adjusted}. \code{p_adjusted} is always present, as \code{NA} when not
#'   computed, so exports from different results always share a column set.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 150
#' coords <- matrix(runif(2 * n), ncol = 2)
#' res <- multiScaleAC(coords, rnorm(n), rnorm(n), h = c(0.05, 0.1, 0.2), perms = 99)
#' res <- compute_significance(res)
#'
#' mx <- export_mxfda(res, sample_id = "RCC4_17")
#' head(mx)
#'
#' \dontrun{
#' # Downstream, in mxfda:
#' # mxfda_obj <- mxfda::add_summary_function(mxfda_obj, mx, metric = "bivariate")
#' }
export_mxfda <- function(result, 
                         sample_id = NULL,
                         interaction_name = NULL) {
  
  # Validate input
  if (!inherits(result, "kernelcomm_result")) {
    stop("result must be a kernelcomm_result object")
  }
  
  # Check that significance was computed
  if (length(result@p_values) == 0) {
    stop("Must call compute_significance() before export_mxfda()")
  }
  
  # Use provided sample_id or default from object
  sample_id <- sample_id %||% result@sample_id
  
  # Create interaction name if not provided
  if (is.null(interaction_name)) {
    interaction_name <- paste0(result@gene1_name, "_", result@gene2_name)
  }
  
  # Build data.frame in mxfda format
  df <- data.frame(
    sample_id = rep(sample_id, length(result@h)),
    interaction_name = rep(interaction_name, length(result@h)),
    r = result@h,
    sumfun = result@observed,
    csr = result@perm_mean,
    fundiff = result@observed - result@perm_mean,
    p_value = result@p_values,
    z_score = result@z_scores,
    null_lower = result@null_lower,
    null_upper = result@null_upper,
    # Always present, NA when not computed. Previously this column was added only when adjusted
    # p-values existed, so exports from different results had different column counts and
    # export_mxfda_batch()'s rbind() failed -- a very likely state after calling
    # compute_significance() on a list without also calling adjust_pvalues().
    p_adjusted = if (length(result@p_adjusted) > 0) {
      result@p_adjusted
    } else {
      rep(NA_real_, length(result@h))
    },
    stringsAsFactors = FALSE
  )

  df
}

# The canonical column set, so every export shares it regardless of what was computed.
.mxfda_columns <- c("sample_id", "interaction_name", "r", "sumfun", "csr", "fundiff",
                    "p_value", "z_score", "null_lower", "null_upper", "p_adjusted")


#' Export Batch Results to mxfda Format
#'
#' Converts a list of kernelcomm_result objects to a single data.frame
#' compatible with mxfda::add_summary_function().
#'
#' @param results List of kernelcomm_result objects (from multiScaleAC_batch)
#' @param sample_ids Character vector of sample identifiers (optional, uses object slots if NULL)
#' @param interaction_names Character vector of interaction names (optional, auto-generated if NULL)
#'
#' @return A single data frame stacking the per-sample exports, with the column
#'   set documented in \code{\link{export_mxfda}}. Binding is on a fixed column
#'   set, so results with and without adjusted p-values combine without error.
#' @export
#'
#' @seealso \code{\link{export_mxfda}}, in particular its note on the bandwidth
#'   occupying mxfda's radius column.
#'
#' @examples
#' set.seed(1)
#' n <- 120
#' mk <- function(seed) {
#'   set.seed(seed)
#'   co <- matrix(runif(2 * n), ncol = 2)
#'   compute_significance(
#'     multiScaleAC(co, rnorm(n), rnorm(n), h = c(0.05, 0.1), perms = 99, verbose = FALSE)
#'   )
#' }
#' results <- list(A = mk(1), B = mk(2))
#'
#' # One has adjusted p-values, the other does not; the bind still succeeds.
#' results$A <- adjust_pvalues(results$A)
#' mx <- export_mxfda_batch(results)
#' str(mx)
export_mxfda_batch <- function(results, 
                                sample_ids = NULL,
                                interaction_names = NULL) {
  
  # Validate input
  if (!is.list(results)) {
    stop("results must be a list")
  }
  
  if (length(results) == 0) {
    stop("results list is empty")
  }
  
  # Check all elements are kernelcomm_result
  if (!all(sapply(results, inherits, "kernelcomm_result"))) {
    stop("All elements in results must be kernelcomm_result objects")
  }
  
  # Check significance computed for all
  if (!all(sapply(results, function(r) length(r@p_values) > 0))) {
    stop("All results must have significance computed. Call compute_significance() on each.")
  }
  
  # Handle sample_ids. Prefers list names over the slot, since results built by multiScaleAC()
  # directly all share the default sample_id and would otherwise be indistinguishable in the output.
  if (is.null(sample_ids)) {
    sample_ids <- .result_labels(results)
  } else {
    if (length(sample_ids) != length(results)) {
      stop("sample_ids length must match results length")
    }
  }
  
  # Handle interaction_names
  if (is.null(interaction_names)) {
    interaction_names <- sapply(results, function(r) {
      paste0(r@gene1_name, "_", r@gene2_name)
    })
  } else {
    if (length(interaction_names) != length(results)) {
      stop("interaction_names length must match results length")
    }
  }
  
  # Export each result and combine
  df_list <- mapply(
    function(result, sid, iname) {
      export_mxfda(result, sample_id = sid, interaction_name = iname)
    },
    results,
    sample_ids,
    interaction_names,
    SIMPLIFY = FALSE
  )
  
  # Bind on a fixed column set. export_mxfda() now always emits the same columns, so this is belt
  # and braces, but it makes the guarantee explicit rather than relying on every caller path.
  df_list <- lapply(df_list, function(d) d[, .mxfda_columns, drop = FALSE])
  combined_df <- do.call(rbind, df_list)
  rownames(combined_df) <- NULL
  combined_df
}
