#' listw to Sparse Matrix
#'
#' @param listw list object of neighbor weights
#' @param n total number of cells in the sample
#'
#' @returns a sparse matrix representing {spdep} listw
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
lag.sparse = function(W_sparse, scaled_gene2){
  as.matrix(W_sparse %*% scaled_gene2)
}

