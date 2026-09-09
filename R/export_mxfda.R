#' Export to mxfda Format
#'
#' Converts kernelcomm_result to data.frame compatible with
#' mxfda::add_summary_function().
#'
#' @param result Object of class \code{kernelcomm_result} (must have significance computed)
#' @param sample_id Character. Sample identifier (overrides object slot)
#' @param interaction_name Character. Gene pair name (e.g., "CD274_PDCD1")
#'
#' @return data.frame with columns: sample_id, interaction_name, r, sumfun, csr, fundiff,
#'         p_value, z_score, null_lower, null_upper
#' @export
#'
#' @examples
#' \dontrun{
#' result <- multiScaleAC(coords, gene1, gene2, h = c(10, 20, 40))
#' result <- compute_significance(result)
#' mxfda_data <- export_mxfda(result, sample_id = "RCC4_17")
#' # mxfda_obj <- mxfda::add_summary_function(mxfda_obj, mxfda_data, metric = "bivariate")
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
    stringsAsFactors = FALSE
  )
  
  # Add adjusted p-values if available
  if (length(result@p_adjusted) > 0) {
    df$p_adjusted <- result@p_adjusted
  }
  
  return(df)
}


#' Export Batch Results to mxfda Format
#'
#' Converts a list of kernelcomm_result objects to a single data.frame
#' compatible with mxfda::add_summary_function().
#'
#' @param results List of kernelcomm_result objects (from multiScaleAC_batch)
#' @param sample_ids Character vector of sample identifiers (optional, uses object slots if NULL)
#' @param interaction_names Character vector of interaction names (optional, auto-generated if NULL)
#'
#' @return data.frame with columns: sample_id, interaction_name, r, sumfun, csr, fundiff,
#'         p_value, z_score, null_lower, null_upper
#' @export
#'
#' @examples
#' \dontrun{
#' # Process multiple samples
#' results <- multiScaleAC_batch(data_list, h_vec = c(10, 20, 40))
#' mxfda_data <- export_mxfda_batch(results)
#' # mxfda_obj <- mxfda::add_summary_function(mxfda_obj, mxfda_data, metric = "bivariate")
#' }
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
  
  # Handle sample_ids
  if (is.null(sample_ids)) {
    sample_ids <- sapply(results, function(r) r@sample_id)
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
  
  # Combine all data.frames
  combined_df <- do.call(rbind, df_list)
  rownames(combined_df) <- NULL
  
  return(combined_df)
}
