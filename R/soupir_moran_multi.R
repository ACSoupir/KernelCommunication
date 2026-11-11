#' Kernel‑based bivariate Moran's I curve (classical definition)
#'
#' @param gene1 numeric vector of expression values for the first gene
#' @param gene2 numeric vector of expression values for the second gene
#' @param coords numeric matrix of spatial coordinates (rows = cells, 2 columns)
#' @param h_vec numeric vector of Gaussian bandwidths
#' @param perms integer number of random permutations (same set used for all bandwidths)
#' @param seed integer random seed (single seed controls all permutations)
#'
#' @return A list with components
#'   `h` – the bandwidth vector,
#'   `observed` – Moran's I for each bandwidth,
#'   `permutations` – matrix `length(h) × perms` of permuted I values,
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' coords <- matrix(runif(2*n, -50, 50), ncol = 2)
#' g1 <- rnorm(n)
#' g2 <- rnorm(n)
#' hvec <- c(5, 15, 30, 60)
#' res <- soupir_moran_multi(coords, g1, g2, h_vec = hvec,
#'                           perms = 500, seed = 123)
#' plot(hvec, res$observed, type = "b")
#'
#' @export
soupir_moran_multi <- function(coords, gene1, gene2, h_vec,
                                   perms, seed=333) {
  N <- nrow(coords)
  g1_s <- scale(gene1)
  g2_s <- scale(gene2)
  denom <- sum(g1_s^2)

  dist_sq <- as.matrix(dist(coords))^2
  set.seed(seed)
  Y_perm <- replicate(perms, sample(g2_s, N, replace = FALSE))

  results <- list()

  for(h in h_vec) {
    #the kernel
    K <- exp(-dist_sq / (2 * h^2))
    #remove self-weight - doens't make biological sense for single cell
    diag(K) <- 0

    # C. Row sums for normalization
    r_sums <- rowSums(K)

    # D. Optimization: x^T W  ==  K %*% (x / r_sums)
    x_weighted <- g1_s / r_sums
    v_h <- as.vector(K %*% x_weighted)

    perm_vals <- (v_h %*% Y_perm) / denom
    obs_val   <- (v_h %*% g2_s) / denom

    results[[paste0("h_", h)]] <- list(observed = obs_val, permutations = as.vector(perm_vals))
  }

  out = list()
  out$h = h_vec
  out$observed = unname(sapply(results, function(x) as.numeric(x$observed)))
  out$observed[is.nan(out$observed)] = 0
  out$permutations = t(sapply(results, function(x) as.numeric(x$permutations))) #rows are h, columns are perms
  out$permutations[is.nan(out$permutations)] = 0
  return(out)
}
