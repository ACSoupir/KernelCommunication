#' Validate Input Data
#'
#' Checks for common issues before computation.
#'
#' @param coords Coordinates matrix
#' @param gene1 Gene 1 expression vector
#' @param gene2 Gene 2 expression vector
#'
#' @return List with components:
#'   \item{errors}{Character vector of critical errors}
#'   \item{warnings}{Character vector of warnings}
#'   \item{valid}{Logical. TRUE if no errors}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' validation <- validate_inputs(coords, gene1, gene2)
#' if (!validation$valid) {
#'   stop("Validation failed: ", paste(validation$errors, collapse = ", "))
#' }
#' if (length(validation$warnings) > 0) {
#'   warning(paste(validation$warnings, collapse = ", "))
#' }
#' }
validate_inputs <- function(coords, gene1, gene2) {
  issues <- list(errors = character(), warnings = character())
  
  # Length checks
  if (!is.matrix(coords) || ncol(coords) != 2) {
    issues$errors <- c(issues$errors, "coords must be a matrix with 2 columns")
  }
  
  if (length(gene1) != nrow(coords)) {
    issues$errors <- c(issues$errors, "gene1 length must match number of coordinates")
  }
  
  if (length(gene2) != nrow(coords)) {
    issues$errors <- c(issues$errors, "gene2 length must match number of coordinates")
  }
  
  # Return early if critical errors
  if (length(issues$errors) > 0) {
    issues$valid <- FALSE
    return(issues)
  }
  
  # Check for duplicate coordinates
  if (any(duplicated(coords))) {
    dup_count <- sum(duplicated(coords))
    issues$warnings <- c(issues$warnings,
                         paste0("Found ", dup_count, " duplicate coordinates. ",
                               "Consider adding small jitter."))
  }
  
  # Zero inflation check
  zero_frac1 <- mean(gene1 == 0)
  zero_frac2 <- mean(gene2 == 0)
  
  if (zero_frac1 > 0.9) {
    issues$warnings <- c(issues$warnings,
                         paste0("gene1 is ", round(zero_frac1 * 100),
                               "% zeros. Consider filtering low-expression genes."))
  }
  
  if (zero_frac2 > 0.9) {
    issues$warnings <- c(issues$warnings,
                         paste0("gene2 is ", round(zero_frac2 * 100),
                               "% zeros. Consider filtering low-expression genes."))
  }
  
  # Variance check
  if (sd(gene1) == 0 || all(gene1 == gene1[1])) {
    issues$errors <- c(issues$errors, "gene1 has zero variance")
  }
  
  if (sd(gene2) == 0 || all(gene2 == gene2[1])) {
    issues$errors <- c(issues$errors, "gene2 has zero variance")
  }
  
  # Missing values
  if (any(is.na(gene1)) || any(is.na(gene2))) {
    issues$errors <- c(issues$errors, "Gene expression contains NA values")
  }
  
  if (any(is.na(coords))) {
    issues$errors <- c(issues$errors, "Coordinates contain NA values")
  }
  
  issues$valid <- length(issues$errors) == 0
  return(issues)
}


#' Suggest Bandwidth Values
#'
#' Recommends bandwidth values based on nearest neighbor distances.
#'
#' @param coords Coordinate matrix (N x 2)
#' @param n_bandwidths Number of bandwidth values to suggest
#' @param range Character. One of "local", "medium", "global"
#'   \itemize{
#'     \item "local": 10th-50th percentile of NN distances (fine-scale)
#'     \item "medium": 25th-75th percentile (intermediate scale)
#'     \item "global": 50th-90th percentile (broad-scale)
#'   }
#'
#' @return Numeric vector of suggested h values
#' @export
#'
#' @examples
#' \dontrun{
#' data(ccRCCExample)
#' coords <- as.matrix(ccRCCExample[, c("x", "y")])
#' h_suggested <- suggest_bandwidths(coords, n_bandwidths = 5, range = "medium")
#' print(h_suggested)
#' }
suggest_bandwidths <- function(coords, n_bandwidths = 5, range = "medium") {
  
  # Validate inputs
  if (!is.matrix(coords) || ncol(coords) != 2) {
    stop("coords must be a matrix with 2 columns")
  }
  
  if (n_bandwidths < 1) {
    stop("n_bandwidths must be positive")
  }
  
  range <- match.arg(range, c("local", "medium", "global"))
  
  # Compute nearest neighbor distances
  dist_matrix <- as.matrix(dist(coords))
  diag(dist_matrix) <- Inf  # Exclude self-distances
  nn1 <- apply(dist_matrix, 1, min)  # 1st nearest neighbor
  
  # Choose quantiles based on range
  if (range == "local") {
    probs <- seq(0.10, 0.50, length.out = n_bandwidths)
  } else if (range == "medium") {
    probs <- seq(0.25, 0.75, length.out = n_bandwidths)
  } else {  # global
    probs <- seq(0.50, 0.90, length.out = n_bandwidths)
  }
  
  # Return quantiles
  h_values <- as.numeric(quantile(nn1, probs = probs, names = FALSE))
  
  return(h_values)
}


#' Compute Neighborhood Quality Metrics
#'
#' Reports statistics about neighborhood structure at a given bandwidth.
#'
#' @param coords Coordinates matrix
#' @param h Bandwidth value
#'
#' @return List with neighborhood statistics:
#'   \item{mean_neighbors}{Average number of neighbors per cell}
#'   \item{sd_neighbors}{Standard deviation of neighbors}
#'   \item{min_neighbors}{Minimum neighbors}
#'   \item{max_neighbors}{Maximum neighbors}
#'   \item{isolated_cells}{Number of cells with no neighbors}
#'   \item{percent_isolated}{Percentage of isolated cells}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(ccRCCExample)
#' coords <- as.matrix(ccRCCExample[, c("x", "y")])
#' stats <- qc_neighborhoods(coords, h = 20)
#' print(stats)
#' }
qc_neighborhoods <- function(coords, h) {
  
  # Compute distance matrix
  dist_matrix <- as.matrix(dist(coords))
  
  # Create adjacency matrix (within bandwidth h)
  W <- (dist_matrix <= h) * 1
  diag(W) <- 0  # Exclude self
  
  # Count neighbors
  n_neighbors <- rowSums(W)
  
  # Compute statistics
  list(
    mean_neighbors = mean(n_neighbors),
    sd_neighbors = sd(n_neighbors),
    min_neighbors = min(n_neighbors),
    max_neighbors = max(n_neighbors),
    isolated_cells = sum(n_neighbors == 0),
    percent_isolated = round(mean(n_neighbors == 0) * 100, 2)
  )
}
