#' listw to Sparse Matrix
#'
#' @param listw list object of neighbor weights
#' @param n total number of cells in the sample
#'
#' @returns a sparse matrix representing \code{spdep} listw
#' @export
#'
listw_to_sparse <- function(listw, n) {
  i <- unlist(mapply(rep, 1:n, sapply(listw$neighbours, length)))
  j <- unlist(listw$neighbours)
  x <- unlist(listw$weights)

  #make sparse matrix for light weight
  W_sparse <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n))
  return(W_sparse)
}

#' Lag Gene 2 with Sparse Weights
#'
#' @param W_sparse weight matrix in sparse format
#' @param scaled_gene2 gene2 scaled to mean 0 and unit variance
#'
#' @returns a matrix of lagged gene 2 values based on the sparse weight matrix
#' @export
#'
lag_sparse = function(W_sparse, scaled_gene2){
  as.matrix(W_sparse %*% scaled_gene2)
}

# No `lag.sparse` compatibility shim is provided, deliberately.
#
# This function was originally named `lag.sparse`, and the pre-refactor NAMESPACE contained only
# `S3method(lag, sparse)` with no corresponding `export()`. It was therefore never callable as
# `lag.sparse(...)`, and the S3 route could not dispatch either: Matrix sparse matrices are S4
# objects of class "dgCMatrix", not "sparse". The function was unreachable, so renaming it to
# `lag_sparse` broke no caller and a deprecation shim would only reintroduce the misleading
# registration that the rename removed.

#' Unique labels for a list of results
#'
#' Internal. Prefers the list names, falls back to the \code{sample_id} slot, and de-duplicates
#' whatever it ends up with.
#'
#' Results produced by \code{multiScaleAC()} directly all carry the default \code{sample_id}
#' "sample1". Functions that used the slot alone therefore built duplicated factor levels from an
#' unnamed list and failed. Only \code{multiScaleAC_batch()} assigns distinct ids.
#'
#' @param results List of kernelcomm_result objects.
#' @return Character vector of unique labels, one per result.
#' @keywords internal
#' @noRd
.result_labels <- function(results) {
  lab <- names(results)
  if (is.null(lab) || any(!nzchar(lab))) {
    slot_ids <- vapply(results, function(r) r@sample_id, character(1))
    if (is.null(lab)) lab <- slot_ids else lab[!nzchar(lab)] <- slot_ids[!nzchar(lab)]
  }
  if (anyDuplicated(lab)) {
    # Suffix rather than error: a user comparing several results built ad hoc has done nothing wrong,
    # and silently collapsing them would be worse than labelling them 1..n.
    dup <- duplicated(lab) | duplicated(lab, fromLast = TRUE)
    lab[dup] <- paste0(lab[dup], "_", seq_len(length(lab))[dup])
  }
  lab
}

# Global variables for ggplot2 and R CMD check
utils::globalVariables(c(
  "null_lower", "null_upper", "h", "n_significant", "null", "perm", 
  "sample_id", "value", "y"
))

