#' multiScaleAC Result Object
#'
#' S4 class for storing multi-scale Moran's I results.
#'
#' @slot h Numeric vector of bandwidths
#' @slot observed Numeric vector of observed Moran's I values
#' @slot permutations Numeric matrix (n_h x n_perms) of permuted values. May have
#'   zero columns when the result was computed with \code{perms = 0}.
#' @slot coords Matrix of original coordinates (for plotting)
#' @slot gene1_name Character. Name of gene 1
#' @slot gene2_name Character. Name of gene 2
#' @slot sample_id Character. Sample identifier
#' @slot seed Numeric. Random seed used
#' @slot n_cells Integer. Number of cells
#' @slot call Language. Function call that created this object
#' @slot p_values Numeric vector of p-values
#' @slot z_scores Numeric vector of standardised effect sizes. These are effect
#'   sizes only; see \code{\link{compute_significance}} for why converting them
#'   to p-values with \code{pnorm()} is unsound for this statistic.
#' @slot perm_mean Numeric vector, mean of the permutation distribution
#' @slot perm_sd Numeric vector, standard deviation of the permutation
#'   distribution
#' @slot null_lower Numeric vector, lower bound of the null envelope. This is a
#'   quantile of the permutation distribution, describing where the statistic
#'   would fall under complete spatial randomness. It is not a confidence
#'   interval for the observed value, which is what the former name
#'   \code{ci_lower} wrongly implied.
#' @slot null_upper Numeric vector, upper bound of the null envelope
#' @slot significant Logical vector, p below the chosen threshold
#' @slot p_adjusted Numeric vector, p-values after multiple-testing correction
#' @slot sup_statistic Numeric scalar, whole-curve supremum statistic
#' @slot sup_p_value Numeric scalar, p-value for the whole curve
#'
#' @section Uncomputed slots:
#' Optional slots are zero-length vectors until computed, not \code{NA}. So
#' \code{length(x@p_values) == 0} means "not computed", which is distinguishable
#' from a computed-but-missing value.
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

    # Populated by compute_significance(). Typed "vector" so they can be empty.
    p_values = "vector",
    z_scores = "vector",
    perm_mean = "vector",
    perm_sd = "vector",
    null_lower = "vector",
    null_upper = "vector",
    significant = "vector",
    p_adjusted = "vector",

    # Whole-curve test: scalars, not per-bandwidth
    sup_statistic = "vector",
    sup_p_value = "vector"
  )
)

# Slots that must be populated together, checked as a group. Grouping them rather than requiring all
# optional slots to be populated at once matters because not every inference route fills every slot:
# a moment-based result has no permutation distribution to take an envelope from. The previous rule
# was all-or-nothing across seven slots, which made such a result impossible to construct.
.sig_slot_groups <- list(
  permutation_inference = c("p_values", "z_scores", "perm_mean", "perm_sd",
                            "null_lower", "null_upper", "significant"),
  supremum = c("sup_statistic", "sup_p_value")
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
  n_h <- length(object@h)

  if (n_h != length(object@observed)) {
    return("'h' and 'observed' must have the same length")
  }
  if (n_h < 1L) {
    return("'h' must contain at least one bandwidth")
  }
  # The previous form was all(dim(permutations) == c(length(h), ncol(permutations))), whose second
  # element compares ncol to itself and is therefore always TRUE. This is what it meant to say.
  if (nrow(object@permutations) != n_h) {
    return("'permutations' must have one row per bandwidth")
  }
  if (nrow(object@coords) != object@n_cells) {
    return("'coords' must have one row per cell, matching 'n_cells'")
  }
  if (length(object@sample_id) != 1L) {
    return("'sample_id' must be a single string")
  }

  for (grp in names(.sig_slot_groups)) {
    slots_in_grp <- .sig_slot_groups[[grp]]
    filled <- vapply(slots_in_grp, function(s) length(slot(object, s)) > 0L, logical(1))
    if (any(filled) && !all(filled)) {
      return(paste0("slots in the '", grp, "' group must be either all empty or all populated; ",
                    "missing: ", paste(slots_in_grp[!filled], collapse = ", ")))
    }
  }

  # Per-bandwidth slots must be per-bandwidth.
  for (s in c("p_values", "z_scores", "perm_mean", "perm_sd",
              "null_lower", "null_upper", "significant", "p_adjusted")) {
    v <- slot(object, s)
    if (length(v) > 0L && length(v) != n_h) {
      return(paste0("'", s, "' must have one value per bandwidth"))
    }
  }
  # Whole-curve slots are scalars.
  for (s in c("sup_statistic", "sup_p_value")) {
    v <- slot(object, s)
    if (length(v) > 0L && length(v) != 1L) {
      return(paste0("'", s, "' describes the whole curve and must be a single value"))
    }
  }

  # [0, 1] rather than (0, 1]: the manuscript's r/n estimator can legitimately return exactly 0 when
  # no permutation is as extreme as the observed value. That is a real hazard downstream -- it breaks
  # log transforms and FDR procedures -- but it is a property of an estimator this package offers, so
  # the class must be able to hold it. compute_significance() warns when it happens instead.
  if (length(object@p_values) > 0L) {
    p <- object@p_values[!is.na(object@p_values)]
    if (length(p) && (any(p < 0) || any(p > 1))) {
      return("'p_values' must lie in [0, 1]")
    }
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
#'
#' @examples
#' res <- create_kernelcomm_result(
#'   h = c(0.05, 0.1),
#'   observed = c(0.02, 0.01),
#'   permutations = matrix(rnorm(20), nrow = 2),
#'   coords = matrix(runif(20), ncol = 2),
#'   gene1_name = "CD274", gene2_name = "PDCD1",
#'   sample_id = "example", seed = 1, n_cells = 10,
#'   call = quote(create_kernelcomm_result())
#' )
#' res
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
