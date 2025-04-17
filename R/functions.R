message("Weight list to sparse matrix")
listw_to_sparse <- function(listw, n) {
  i <- unlist(mapply(rep, 1:n, sapply(listw$neighbours, length)))
  j <- unlist(listw$neighbours)
  x <- unlist(listw$weights)

  #make sparse matrix for light weight
  W_sparse <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n))
  return(W_sparse)
}

message("Spatial lag")
lag.sparse = function(W_sparse, scaled_gene2){
  as.matrix(W_sparse %*% scaled_gene2)
}

#' Kernel-Based Bivariate Moran's I
#'
#' Computes a bivariate kernel-based Moran's I statistic between two gene expression vectors using a Gaussian kernel applied to a squared distance matrix.
#' The function returns the observed statistic and a distribution of permuted statistics to assess significance.
#'
#' @param h Numeric bandwidth parameter for the Gaussian kernel (in same units as `dist_matrix`).
#' @param dist_matrix A square numeric matrix of pairwise distances between observations (e.g., cells).
#' @param gene1 A numeric vector of expression values for the first gene (length equal to number of rows in `dist_matrix`).
#' @param gene2 A numeric vector of expression values for the second gene.
#' @param perms Number of permutations to compute the null distribution.
#'
#' @return A list with two components:
#' \describe{
#'   \item{observed}{Observed Moran's I value.}
#'   \item{permutations}{Vector of permuted Moran's I values.}
#' }
#'
#' @examples
#' set.seed(42)
#' d <- as.matrix(dist(matrix(rnorm(20), ncol = 2)))
#' g1 <- rnorm(10)
#' g2 <- rnorm(10)
#' result <- soupir_moran_bv(h = 1.5, dist_matrix = d, gene1 = g1, gene2 = g2, perms = 100)
#' result$observed
#' hist(result$permutations)
#'
#' @export
soupir_moran_bv <- function(h, dist_matrix, gene1, gene2, perms) {
  dist_sq = dist_matrix^2
  g1_s = scale(gene1)
  g2_s = scale(gene2)
  denom = sum(g1_s^2)
  P = replicate(perms,
                 sample(g2_s, length(g2_s), replace = FALSE))
  K = function(d, h) {
    exp(- (as.matrix(d)) / (2 * h^2)) * (1 - diag(dim(d)[1]))
  }  #gaussian kernel

  listw_kernel = K(dist_sq, h)
  #diag(listw_kernel) = 0

  nz_idx = which(as.matrix(listw_kernel) != 0, arr.ind = TRUE)
  neighbors <- split(nz_idx[,2], nz_idx[,1])
  class(neighbors) = c("nb")
  weights = lapply(1:nrow(listw_kernel), function(r){
    v = listw_kernel[r,neighbors[[r]]]
    v/sum(v)
  })
  lw_kernel = list(style = "W",
                   neighbours = neighbors,
                   weights = weights)
  class(lw_kernel) = c("listw", "nb")
  Ws = listw_to_sparse(lw_kernel, nrow(dist_sq))

  u_obs = Ws %*% g2_s
  obs_val = sum(u_obs * g1_s) / denom

  U_perm = sparse_dense_mult(Ws, P)
  perm_vals = colSums(U_perm * matrix(g1_x, nrow = length(g1_s), ncol = ncol(U_perm))) / denom

  return(list(observed = obs_val,
       permutations = perm_vals))
}
