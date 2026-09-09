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

# Global variables for ggplot2 and R CMD check
utils::globalVariables(c(
  "ci_lower", "ci_upper", "h", "n_significant", "null", "perm", 
  "sample_id", "value", "y"
))

