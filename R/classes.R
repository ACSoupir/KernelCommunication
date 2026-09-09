#' KernelCommunication Result Object
#'
#' S4 class for storing multi-scale Moran's I computation results.
#'
#' @slot h Numeric vector of bandwidths
#' @slot observed Numeric vector of observed Moran's I values
#' @slot permutations Numeric matrix [n_h × n_perms] of permutation results
#' @slot coords Matrix of original coordinates (for plotting)
#' @slot gene1_name Character. Name of gene 1
#' @slot gene2_name Character. Name of gene 2
#' @slot sample_id Character. Sample identifier
#' @slot seed Numeric. Random seed used
#' @slot n_cells Integer. Number of cells
#' @slot call Language. Function call that created this object
#' @slot p_values Numeric vector of empirical p-values (NULL until computed)
#' @slot z_scores Numeric vector of z-scores (NULL until computed)
#' @slot perm_mean Numeric vector, mean of permutations (NULL until computed)
#' @slot perm_sd Numeric vector, SD of permutations (NULL until computed)
#' @slot ci_lower Numeric vector, lower CI bound (NULL until computed)
#' @slot ci_upper Numeric vector, upper CI bound (NULL until computed)
#' @slot significant Logical vector, p < alpha (NULL until computed)
#' @slot p_adjusted Numeric vector, adjusted p-values (NULL until computed)
#'
#' @aliases kernelcomm_result-class
#' @exportClass kernelcomm_result
setClass("kernelcomm_result",
  slots = list(
    # Required slots (always present)
    h = "numeric",
    observed = "numeric",
    permutations = "matrix",
    coords = "matrix",
    gene1_name = "character",
    gene2_name = "character",
    sample_id = "character",
    seed = "numeric",
    n_cells = "integer",
    call = "language",
    
    # Optional slots (NULL until compute_significance() is called)
    # Using "vector" type to allow NULL or numeric/logic
    p_values = "vector",
    z_scores = "vector",
    perm_mean = "vector",
    perm_sd = "vector",
    ci_lower = "vector",
    ci_upper = "vector",
    significant = "vector",
    p_adjusted = "vector"
  )
)

#' Validity check for kernelcomm_result
#'
#' Ensures object structure is consistent.
#'
#' @name kernelcomm_result-validity
#' @param object A kernelcomm_result object
#' @return TRUE if valid, error message otherwise
#' @keywords internal
setValidity("kernelcomm_result", function(object) {
  # Required slot validation
  if (length(object@h) != length(object@observed)) {
    return("'h' and 'observed' must have same length")
  }
  
  if (!all(dim(object@permutations) == c(length(object@h), ncol(object@permutations)))) {
    return("'permutations' rows must match length of 'h'")
  }
  
  if (nrow(object@coords) != object@n_cells) {
    return("coords rows should match n_cells")
  }
  
  # Optional slot validation: check if significance has been computed
  # Slots initialized as empty vectors (length 0) are treated as NULL
  sig_slots <- c("p_values", "z_scores", "perm_mean", "perm_sd", 
                 "ci_lower", "ci_upper", "significant")
  
  # Check which slots have actual data (length > 0)
  sig_has_data <- sapply(sig_slots, function(s) length(slot(object, s)) > 0)
  
  # If some have data and some don't, that's an error
  if (any(sig_has_data) && !all(sig_has_data)) {
    return("All significance slots must be empty or all must have data")
  }
  
  # If significance is computed, check lengths match
  if (length(object@p_values) > 0 && length(object@p_values) != length(object@h)) {
    return("Significance slots must match length of h")
  }
  
  TRUE
})

#' Create a kernelcomm_result object
#'
#' Helper function to create a valid kernelcomm_result object.
#'
#' @param h Numeric vector of bandwidths
#' @param observed Numeric vector of observed Moran's I
#' @param permutations Numeric matrix of permutation results
#' @param coords Matrix of coordinates
#' @param gene1_name Character. Gene 1 name
#' @param gene2_name Character. Gene 2 name
#' @param sample_id Character. Sample identifier
#' @param seed Numeric. Random seed
#' @param n_cells Integer. Number of cells
#' @param call Language. Function call
#'
#' @return kernelcomm_result object
#' @export
create_kernelcomm_result <- function(h, observed, permutations, coords,
                                      gene1_name, gene2_name, sample_id,
                                      seed, n_cells, call) {
  new("kernelcomm_result",
    h = h,
    observed = observed,
    permutations = permutations,
    coords = coords,
    gene1_name = gene1_name,
    gene2_name = gene2_name,
    sample_id = sample_id,
    seed = seed,
    n_cells = as.integer(n_cells),
    call = call
  )
}
