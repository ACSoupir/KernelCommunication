#' Kernel-Based Bivariate Moran's I
#'
#' Computes a bivariate kernel-based Moran's I statistic between two gene expression vectors using a Gaussian kernel applied to a squared distance matrix.
#' The function returns the observed statistic and a distribution of permuted statistics to assess significance.
#'
#' @param h Numeric bandwidth parameter for the Gaussian kernel (in same units as `dist_matrix`).
#' @param coords A matrix with 2 columns (x and y locations) and rows as cells.
#' @param gene1 A numeric vector of expression values for the first gene (length equal to number of rows in `dist_matrix`).
#' @param gene2 A numeric vector of expression values for the second gene.
#' @param perms integer of permutations to compute the null distribution.
#' @param eps numeric value for the lower end threshold of the gaussian kernel. Values below this will be ignored when calculating neighbor weights. Default is `0`, i.e. inclusion of all neighbors in weight matrix.
#'
#' @return A list with two components:
#' \describe{
#'   \item{observed}{Observed Moran's I value.}
#'   \item{permutations}{Vector of permuted Moran's I values.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' coor <- matrix(rnorm(20), ncol = 2)
#' g1 <- rnorm(10)
#' g2 <- rnorm(10)
#' result <- soupir_moran_bv(h = 1.5, coords = coor, gene1 = g1, gene2 = g2, perms = 100)
#' result$observed
#' hist(result$permutations)
#' }
#'
#' @export
soupir_moran_bv <- function(h, coords, gene1, gene2, perms, eps = 0) {
  if(length(h) != 1) stop("`h` must be single value")
  g1_s = scale(gene1)
  g2_s = scale(gene2)
  denom = sum(g1_s^2)
  P = replicate(perms,
                 sample(g2_s, length(g2_s), replace = FALSE))

  Ws = gaussian_kernel_weights(coords, h =h, eps = eps)

  u_obs = Ws %*% g2_s
  obs_val = sum(u_obs * g1_s) / denom

  U_perm = sparse_dense_mult(Ws, P)
  perm_vals = colSums(U_perm * matrix(g1_s, nrow = length(g1_s), ncol = ncol(U_perm))) / denom

  return(list(observed = obs_val,
       permutations = perm_vals))
}
