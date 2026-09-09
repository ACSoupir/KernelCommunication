#' Multi-scale Spatial Autocorrelation Analysis
#'
#' Computes bivariate Moran's I across multiple kernel bandwidth scales using
#' RcppEigen for efficient sparse matrix operations. Uses permutation-based
#' inference with optimized handling for large datasets.
#'
#' @param coords numeric matrix of spatial coordinates (N x 2)
#' @param gene1 numeric vector of expression values for the first variable (length N)
#' @param gene2 numeric vector of expression values for the second variable (length N)
#' @param h numeric vector of kernel bandwidths (must be > 0)
#' @param perms integer number of permutations (default: 1000)
#' @param seed integer random seed for reproducibility (default: 333)
#' @param eps numeric epsilon threshold for weight pruning (default: 0 = no pruning)
#' @param batch_size integer permutations per batch for memory optimization (default: 100)
#' @param verbose logical show progress bar (default: NULL = auto for large datasets)
#' @param large_mode_n integer N threshold for automatic sparse mode (default: 10000)
#' @param kernel character or function. If character: built-in kernel ("gaussian" or "uniform").
#'   If function: custom kernel function with signature \code{kernel(dist_matrix, h, ...)} where
#'   \code{dist_matrix} is an N x N matrix of pairwise distances and \code{h} is the bandwidth.
#'   The function must return an N x N weight matrix (dense or sparse Matrix::dsCMatrix).
#' @param distance character or function. If character: built-in distance method ("euclidean" or "manhattan").
#'   If function: custom distance function with signature \code{distance(coords, ...)} where
#'   \code{coords} is an N x 2 matrix. The function must return an N x N distance matrix.
#' @param ... additional parameters passed to distance and kernel functions (e.g., \code{p} for Minkowski distance)
#'
#' @return An S4 object of class \linkS4class{kernelcomm_result} with slots:
#' \describe{
#'   \item{h}{numeric vector of bandwidth values used}
#'   \item{observed}{numeric vector of Moran's I for each bandwidth}
#'   \item{permutations}{numeric matrix (n_h x perms) of permuted I values}
#'   \item{coords}{matrix of original coordinates}
#'   \item{gene1_name}{character. Name of gene 1}
#'   \item{gene2_name}{character. Name of gene 2}
#'   \item{sample_id}{character. Sample identifier}
#'   \item{seed}{numeric. Random seed used}
#'   \item{n_cells}{integer. Number of cells}
#'   \item{call}{language. Function call that created this object}
#'   \item{p_values, z_scores, perm_mean, perm_sd, ci_lower, ci_upper, significant, p_adjusted}{NULL (computed by \code{\link{compute_significance}})}
#' }
#'
#' @details
#' \strong{Built-in Kernels:}
#' \itemize{
#'   \item \code{"gaussian"}: \eqn{w_{ij} = \exp(-d_{ij}^2 / (2h^2))} where \eqn{d_{ij}} is the distance
#'   \item \code{"uniform"}: \eqn{w_{ij} = 1} if \eqn{d_{ij} \le h}, else 0
#' }
#'
#' \strong{Built-in Distance Methods:}
#' \itemize{
#'   \item \code{"euclidean"}: \eqn{\sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}}
#'   \item \code{"manhattan"}: \eqn{|x_i - x_j| + |y_i - y_j|}
#' }
#'
#' \strong{Custom Kernel Functions:}
#' Custom kernel functions must have signature: \code{kernel(dist_matrix, h, ...)} where:
#' \itemize{
#'   \item \code{dist_matrix}: N x N matrix of pairwise distances
#'   \item \code{h}: scalar bandwidth for current computation
#'   \item \code{...}: additional parameters passed from multiScaleAC
#' }
#' The function must return an N x N weight matrix (can be dense matrix or sparse Matrix::dsCMatrix).
#'
#' \strong{Custom Distance Functions:}
#' Custom distance functions must have signature: \code{distance(coords, ...)} where:
#' \itemize{
#'   \item \code{coords}: N x 2 matrix of coordinates
#'   \item \code{...}: additional parameters (e.g., \code{p} for Minkowski distance)
#' }
#' The function must return an N x N distance matrix.
#'
#' \strong{Performance Notes:}
#' Custom kernel/distance functions are computed in R and disable C++ optimizations.
#' For large datasets (N > 10,000), using built-in kernels and distances is recommended.
#'
#' @examples
#' \dontrun{
#' # Basic usage with permutation test
#' set.seed(1)
#' n <- 500
#' coords <- matrix(runif(2 * n, -50, 50), ncol = 2)
#' g1 <- rnorm(n)
#' g2 <- rnorm(n)
#' h_vec <- c(5, 10, 20, 40)
#'
#' result <- multiScaleAC(
#'   coords = coords,
#'   gene1 = g1,
#'   gene2 = g2,
#'   h = h_vec,
#'   perms = 999,
#'   seed = 123
#' )
#'
#' # Custom kernel (exponential decay)
#' exp_kernel <- function(dist_matrix, h) {
#'   exp(-dist_matrix / h)
#' }
#' result_custom <- multiScaleAC(
#'   coords = coords,
#'   gene1 = g1,
#'   gene2 = g2,
#'   h = h_vec,
#'   kernel = exp_kernel
#' )
#'
#' # Custom distance (Minkowski with p=3) and custom kernel
#' minkowski_dist <- function(coords, p = 3) {
#'   as.matrix(dist(coords, method = "minkowski", p = p))
#' }
#' result_mink <- multiScaleAC(
#'   coords = coords,
#'   gene1 = g1,
#'   gene2 = g2,
#'   h = h_vec,
#'   distance = minkowski_dist,
#'   p = 3
#' )
#' }
#'
#' @seealso \code{\link{soupir_moran_bv}} for single-bandwidth analysis,
#' \code{\link{soupir_moran_multi}} for the deprecated multi-bandwidth function
#'
#' @export
multiScaleAC <- function(coords,
                         gene1,
                         gene2,
                         h,
                         perms = 1000,
                         seed = 333,
                         eps = 0,
                         batch_size = 100,
                         verbose = NULL,
                         large_mode_n = 10000,
                         kernel = c("gaussian", "uniform"),
                         distance = c("euclidean", "manhattan"),
                         ...) {

  # Capture additional parameters
  extra_params <- list(...)
  
  # Validate and normalize kernel parameter
  kernel_is_custom <- is.function(kernel)
  if (!kernel_is_custom) {
    kernel <- match.arg(kernel)
  }
  
  # Validate and normalize distance parameter
  distance_is_custom <- is.function(distance)
  if (!distance_is_custom) {
    distance <- match.arg(distance)
  }

  # Validate inputs
  .validate_inputs(coords, gene1, gene2, h, perms,
                   eps, batch_size, large_mode_n, seed,
                   kernel, distance, kernel_is_custom, distance_is_custom)

  # Get N
  N <- nrow(coords)

  # Sanitize h (filter negative/zero values with warnings)
  h_clean <- .sanitize_h(h)

  # Determine large mode (auto-disable for custom functions)
  large_mode <- N > large_mode_n && !kernel_is_custom && !distance_is_custom

  # Auto-enable verbose for large datasets
  if (is.null(verbose)) {
    verbose <- large_mode
  }

  # Warn about large mode disable for custom functions
  if (N > large_mode_n && kernel_is_custom) {
    warning("Custom kernel disables large dataset optimization; computation may be slow")
  }
  if (N > large_mode_n && distance_is_custom) {
    warning("Custom distance method disables large dataset optimization; computation may be slow")
  }

  # Center and scale genes
  g1_s <- as.vector(scale(gene1))
  g2_s <- as.vector(scale(gene2))

  # Generate permutation matrix ONCE (for all paths to ensure RNG compatibility)
  Y_perm <- NULL
  if (perms > 0) {
    set.seed(seed)
    Y_perm <- replicate(perms, sample(g2_s, N, replace = FALSE))
  }

  # Call C++ function for built-in kernels/distances
  if (!kernel_is_custom && !distance_is_custom) {
    result <- moran_multi_cpp(
      coords = coords,
      gene1 = g1_s,
      gene2 = g2_s,
      h_vec = h_clean,
      Y_perm = Y_perm,
      perms = perms,
      seed = seed,
      eps = eps,
      batch_size = batch_size,
      verbose = verbose,
      large_mode_n = large_mode_n,
      kernel = kernel,
      dist_method = distance
    )
  } else {
    # Custom kernel/distance path (R-based)
    result <- .compute_custom_kernel_moran(
      coords = coords,
      gene1 = g1_s,
      gene2 = g2_s,
      h_vec = h_clean,
      Y_perm = Y_perm,
      perms = perms,
      eps = eps,
      distance = distance,
      kernel = kernel,
      extra_params = extra_params,
      verbose = verbose
    )
  }

  # Create S4 object with raw results only
  .create_result_object(result, coords, gene1, gene2, seed, call = match.call())
}


# Internal helper functions -------------------------------------------------

#' Validate input parameters
#'
#' @noRd
.validate_inputs <- function(coords, gene1, gene2, h, perms,
                             eps, batch_size, large_mode_n, seed,
                             kernel, distance, kernel_is_custom, distance_is_custom) {

  # Validate coords
  if (!is.matrix(coords) || ncol(coords) != 2) {
    stop("coords must be a numeric matrix with 2 columns")
  }
  if (!is.numeric(coords)) {
    stop("coords must be numeric")
  }
  if (any(is.na(coords)) || any(is.infinite(coords))) {
    stop("Coordinates contain missing or infinite values")
  }

  # Validate gene1 and gene2
  if (!is.numeric(gene1) || !is.numeric(gene2)) {
    stop("gene1 and gene2 must be numeric vectors")
  }
  if (any(is.na(gene1)) || any(is.infinite(gene1)) ||
      any(is.na(gene2)) || any(is.infinite(gene2))) {
    stop("Gene expression values contain missing or infinite values")
  }

  # Check lengths match
  N <- nrow(coords)
  if (length(gene1) != N || length(gene2) != N) {
    stop("coords, gene1, and gene2 must have same number of observations")
  }

  # Check minimum N
  if (N < 3) {
    stop("At least 3 observations required")
  }

  # Validate h
  if (!is.numeric(h) || length(h) < 1) {
    stop("h must contain at least one bandwidth value")
  }

  # Validate perms
  if (!is.numeric(perms) || perms < 1 || perms != as.integer(perms)) {
    stop("perms must be a positive integer")
  }

  # Validate seed
  if (!is.numeric(seed)) {
    stop("seed must be an integer")
  }

  # Validate eps
  if (!is.numeric(eps) || eps < 0) {
    stop("eps must be a non-negative numeric value")
  }

  # Validate batch_size
  if (!is.numeric(batch_size) || batch_size < 1 || batch_size != as.integer(batch_size)) {
    stop("batch_size must be a positive integer")
  }

  # Validate large_mode_n
  if (!is.numeric(large_mode_n) || large_mode_n < 3 || large_mode_n != as.integer(large_mode_n)) {
    stop("large_mode_n must be an integer >= 3")
  }

  # Validate kernel
  if (!kernel_is_custom) {
    if (!kernel %in% c("gaussian", "uniform")) {
      stop("kernel must be 'gaussian', 'uniform', or a function")
    }
  } else {
    # kernel is already verified to be a function
  }

  # Validate distance
  if (!distance_is_custom) {
    if (!distance %in% c("euclidean", "manhattan")) {
      stop("distance must be 'euclidean', 'manhattan', or a function")
    }
  } else {
    # distance is already verified to be a function
  }
}


#' Sanitize bandwidth vector (filter negative/zero values)
#'
#' @noRd
.sanitize_h <- function(h) {
  h_clean <- h

  # Remove negative values
  neg_count <- sum(h < 0)
  if (neg_count > 0) {
    warning(paste0("Removed ", neg_count, " negative bandwidth values"))
    h_clean <- h_clean[h_clean >= 0]
  }

  # Remove zero values
  zero_count <- sum(h_clean == 0)
  if (zero_count > 0) {
    warning(paste0("Removed ", zero_count, " zero bandwidth values (undefined kernel)"))
    h_clean <- h_clean[h_clean > 0]
  }

  # Check at least one valid
  if (length(h_clean) < 1) {
    stop("No valid bandwidth values remaining")
  }

  return(h_clean)
}


#' Compute Moran's I with custom kernel/distance (R-based)
#'
#' @noRd
.compute_custom_kernel_moran <- function(coords, gene1, gene2, h_vec, Y_perm,
                                         perms, eps,
                                         distance, kernel, extra_params,
                                         verbose) {

  N <- nrow(coords)
  g1_s <- gene1
  g2_s <- gene2
  denom <- sum(g1_s^2)

  # Compute distance matrix
  if (is.function(distance)) {
    # Custom distance function: pass coords and extra_params
    # Filter extra_params to only those that distance function accepts
    dist_formals <- names(formals(distance))
    dist_params <- extra_params[names(extra_params) %in% dist_formals]
    dist_matrix <- do.call(distance, c(list(coords), dist_params))
  } else {
    # Built-in distance method
    if (distance == "euclidean") {
      dist_matrix <- as.matrix(dist(coords))
    } else if (distance == "manhattan") {
      dist_matrix <- as.matrix(dist(coords, method = "manhattan"))
    }
  }
  
  # Validate distance matrix
  if (!is.matrix(dist_matrix) || nrow(dist_matrix) != N || ncol(dist_matrix) != N) {
    stop("distance function must return a square matrix with dimensions N x N")
  }

  n_h <- length(h_vec)
  observed <- numeric(n_h)
  permutations <- matrix(0, nrow = n_h, ncol = perms)

  # Progress reporting
  if (verbose && n_h > 1) {
    progress_interval <- max(1, n_h %/% 20)
  }

  for (h_idx in seq_along(h_vec)) {
    h <- h_vec[h_idx]

    # Compute kernel weights
    if (is.function(kernel)) {
      # Custom kernel function: pass dist_matrix, h, and only params kernel accepts
      kernel_formals <- names(formals(kernel))
      kernel_params <- extra_params[names(extra_params) %in% kernel_formals]
      W <- do.call(kernel, c(list(dist_matrix, h), kernel_params))
    } else {
      # Built-in kernel (should not reach here for custom path, but for safety)
      stop("Built-in kernels should use C++ path")
    }

    # Convert to sparse if needed
    if (!inherits(W, "dsCMatrix")) {
      W <- Matrix::Matrix(W, sparse = TRUE)
    }

    # Row-normalize
    row_sums <- Matrix::rowSums(W)
    row_sums[row_sums == 0] <- 1  # Avoid division by zero
    W_norm <- W / row_sums

    # Compute observed Moran's I
    lag_g2 <- as.vector(W_norm %*% g2_s)
    observed[h_idx] <- sum(g1_s * lag_g2) / denom

    # Compute permutations
    if (perms > 0) {
      # Batch processing
      n_batches <- ceiling(perms / 100)
      batch_size <- 100

      for (b in seq_len(n_batches)) {
        batch_start <- (b - 1) * batch_size + 1
        batch_end <- min(b * batch_size, perms)
        current_batch_size <- batch_end - batch_start + 1

        Y_batch <- Y_perm[, batch_start:batch_end, drop = FALSE]
        U_batch <- W_norm %*% Y_batch
        perm_vals <- as.vector(crossprod(g1_s, U_batch)) / denom
        permutations[h_idx, batch_start:batch_end] <- perm_vals
      }
    }

    # Progress reporting
    if (verbose && n_h > 1) {
      if (h_idx %% progress_interval == 0 || h_idx == n_h) {
        pct <- h_idx * 100 / n_h
        cat(sprintf("\r  Computing bandwidths: %d%% (%d/%d)", pct, h_idx, n_h))
        if (h_idx == n_h) cat("\n")
      }
    }
  }

  list(
    h = h_vec,
    observed = observed,
    permutations = permutations
  )
}


# #' Compute analytical variance using numerically stable formula
# #'
# #' Formula: Var(I) = ||W'x||² / (N-1)²
# #' where W'x = t(W) %*% x, avoiding catastrophic cancellation
# #'
# #' @param W_norm row-normalized weight matrix (N x N, dense or sparse)
# #' @param N number of observations
# #' @param g1 standardized gene1 vector (length N)
# #' @return scalar variance value
# #'
# #' @noRd
# .compute_analytical_variance_direct <- function(W_norm, N, g1) {
#   # Direct computation: Var(I) = ||W'x||² / (N-1)²
#   
#   # For sparse matrices, use Matrix::t() for transpose
#   if (inherits(W_norm, "Matrix")) {
#     # Sparse or Matrix class
#     Wtx <- as.vector(Matrix::t(W_norm) %*% g1)
#   } else {
#     # Regular dense matrix
#     Wtx <- as.vector(t(W_norm) %*% g1)
#   }
#   
#   # Squared norm and variance
#   Wtx_sq_norm <- sum(Wtx^2)
#   variance <- Wtx_sq_norm / ((N - 1)^2)
#   
#   return(max(variance, 0))  # Ensure non-negative
# }


#' Create kernelcomm_result S4 object
#'
#' Helper function to create a valid kernelcomm_result object from computation results.
#'
#' @param result List with h, observed, permutations from computation
#' @param coords Matrix of coordinates
#' @param gene1 Gene 1 expression vector
#' @param gene2 Gene 2 expression vector
#' @param seed Random seed used
#' @param call Function call that created this object
#'
#' @return kernelcomm_result S4 object
#' @noRd
.create_result_object <- function(result, coords, gene1, gene2, seed, call) {
  # Get gene names from expressions (deparse substitute)
  gene1_name <- tryCatch(deparse(substitute(gene1)), error = function(e) "gene1")
  gene2_name <- tryCatch(deparse(substitute(gene2)), error = function(e) "gene2")
  
  # Create S4 object
  create_kernelcomm_result(
    h = result$h,
    observed = result$observed,
    permutations = result$permutations,
    coords = as.matrix(coords),
    gene1_name = gene1_name,
    gene2_name = gene2_name,
    sample_id = "sample1",  # Default, can be overridden later
    seed = seed,
    n_cells = nrow(coords),
    call = call
  )
}


#' Format output list (deprecated - kept for backward compatibility)
#'
#' @noRd
.format_output <- function(result) {
  # This function is deprecated. Use .create_result_object() instead.
  
  output <- list(
    h = result$h,
    observed = result$observed,
    permutations = result$permutations,
    method = "permutation"
  )

  return(output)
}


#' Batch Processing for Multiple Samples
#'
#' Processes multiple samples in parallel or sequentially, returning a list of
#' kernelcomm_result objects. Supports checkpointing for long-running batches.
#'
#' @param data_list List of data frames, each containing coordinates and gene expressions.
#'   Each data frame must have columns specified in coords_cols, gene1_col, gene2_col.
#' @param sample_ids Character vector of sample identifiers (optional, auto-generated if NULL)
#' @param h_vec Numeric vector of bandwidths
#' @param perms Integer number of permutations
#' @param seed Integer random seed
#' @param n_cores Integer number of CPU cores for parallel processing (default: 1 = sequential)
#' @param checkpoint_file Character path to checkpoint file (default: NULL = no checkpointing)
#' @param coords_cols Character vector of two column names for coordinates (default: c("x", "y"))
#' @param gene1_col Character column name for gene 1 (default: "CD274")
#' @param gene2_col Character column name for gene 2 (default: "PDCD1")
#' @param ... Additional parameters passed to multiScaleAC()
#'
#' @return List of \linkS4class{kernelcomm_result} objects, one per sample
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare data for multiple samples
#' data_list <- list(
#'   sample1 = ccRCCExample,
#'   sample2 = ccRCCExample  # In practice, different samples
#' )
#'
#' # Process all samples (sequential)
#' results <- multiScaleAC_batch(
#'   data_list = data_list,
#'   h_vec = c(10, 20, 40),
#'   perms = 99
#' )
#'
#' # Process with parallel computing (4 cores)
#' results_parallel <- multiScaleAC_batch(
#'   data_list = data_list,
#'   h_vec = c(10, 20, 40),
#'   perms = 99,
#'   n_cores = 4
#' )
#'
#' # With checkpointing for very large batches
#' results_checkpoint <- multiScaleAC_batch(
#'   data_list = data_list,
#'   h_vec = c(10, 20, 40),
#'   checkpoint_file = "batch_checkpoint.rds"
#' )
#' }
#'
#' @seealso \code{\link{export_mxfda_batch}} for exporting batch results
multiScaleAC_batch <- function(data_list,
                                sample_ids = NULL,
                                h_vec,
                                perms = 1000,
                                seed = 333,
                                n_cores = 1,
                                checkpoint_file = NULL,
                                coords_cols = c("x", "y"),
                                gene1_col = "CD274",
                                gene2_col = "PDCD1",
                                ...) {
  
  # Validate input list
  if (!is.list(data_list) || length(data_list) == 0) {
    stop("data_list must be a non-empty list")
  }
  
  # Auto-generate sample_ids if not provided
  if (is.null(sample_ids)) {
    sample_ids <- names(data_list)
    if (is.null(sample_ids) || all(sample_ids == "")) {
      sample_ids <- paste0("sample", seq_along(data_list))
    }
  } else {
    if (length(sample_ids) != length(data_list)) {
      stop("sample_ids length must match data_list length")
    }
  }
  
  # Validate checkpoint_file if provided
  if (!is.null(checkpoint_file)) {
    if (!is.character(checkpoint_file) || length(checkpoint_file) != 1) {
      stop("checkpoint_file must be a single character path")
    }
  }
  
  # Load checkpoint if exists
  checkpoint <- NULL
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    cat("Loading checkpoint from", checkpoint_file, "\n")
    checkpoint <- readRDS(checkpoint_file)
    if (!is.list(checkpoint) || length(checkpoint) == 0) {
      checkpoint <- NULL
      cat("Checkpoint file is empty, starting fresh\n")
    } else {
      cat("Resuming from sample", length(checkpoint) + 1, "of", length(data_list), "\n")
    }
  }
  
  # Define worker function
  process_sample <- function(i) {
    sample_id <- sample_ids[i]
    data <- data_list[[i]]
    
    # Check if already processed (checkpoint)
    if (!is.null(checkpoint) && i <= length(checkpoint)) {
      cat("Skipping sample", i, "(", sample_id, "): already done\n")
      return(checkpoint[[i]])
    }
    
    cat("Processing sample", i, "of", length(data_list), ":", sample_id, "\n")
    
    # Extract coordinates and genes
    coords <- as.matrix(data[, coords_cols, drop = FALSE])
    gene1 <- data[[gene1_col]]
    gene2 <- data[[gene2_col]]
    
    # Check required columns exist
    if (!all(coords_cols %in% colnames(data))) {
      stop("Sample ", sample_id, " missing coordinate columns: ", 
           paste(setdiff(coords_cols, colnames(data)), collapse = ", "))
    }
    if (!gene1_col %in% colnames(data)) {
      stop("Sample ", sample_id, " missing gene1 column: ", gene1_col)
    }
    if (!gene2_col %in% colnames(data)) {
      stop("Sample ", sample_id, " missing gene2 column: ", gene2_col)
    }
    
    # Run multiScaleAC
    result <- multiScaleAC(
      coords = coords,
      gene1 = gene1,
      gene2 = gene2,
      h = h_vec,
      perms = perms,
      seed = seed + i,  # Different seed per sample
      ...
    )
    
    # Update sample_id in result
    result@sample_id <- sample_id
    
    # Save checkpoint if enabled
    if (!is.null(checkpoint_file)) {
      # Append to checkpoint
      if (is.null(checkpoint)) {
        checkpoint <- list()
      }
      checkpoint[[i]] <<- result
      saveRDS(checkpoint, checkpoint_file)
    }
    
    # Clear memory
    gc()
    
    return(result)
  }
  
  # Process samples
  n_samples <- length(data_list)
  results <- vector("list", n_samples)
  
  if (n_cores > 1 && .Platform$OS.type == "unix") {
    # Parallel processing with mclapply (Unix/Mac only)
    cat("Processing", n_samples, "samples in parallel with", n_cores, "cores\n")
    
    # Determine which samples to process
    if (is.null(checkpoint)) {
      indices <- seq_len(n_samples)
    } else {
      indices <- seq(length(checkpoint) + 1, n_samples)
      # Copy checkpoint results
      for (i in seq_along(checkpoint)) {
        results[[i]] <- checkpoint[[i]]
      }
    }
    
    # Run parallel processing
    parallel_results <- parallel::mclapply(
      indices,
      process_sample,
      mc.cores = n_cores,
      mc.preschedule = TRUE
    )
    
    # Store results
    for (j in seq_along(indices)) {
      results[[indices[j]]] <- parallel_results[[j]]
    }
    
  } else {
    # Sequential processing
    if (n_cores > 1) {
      warning("Parallel processing only supported on Unix/Mac. Using sequential processing.")
    }
    
    cat("Processing", n_samples, "samples sequentially\n")
    
    # Determine starting index
    start_idx <- 1
    if (!is.null(checkpoint)) {
      start_idx <- length(checkpoint) + 1
      # Copy checkpoint results
      for (i in seq_along(checkpoint)) {
        results[[i]] <- checkpoint[[i]]
      }
    }
    
    # Run sequential processing
    for (i in start_idx:n_samples) {
      results[[i]] <- process_sample(i)
    }
  }
  
  # Remove checkpoint file on successful completion
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    file.remove(checkpoint_file)
    cat("Checkpoint file removed (batch complete)\n")
  }
  
  # Name the results list
  names(results) <- sample_ids
  
  cat("Batch processing complete:", length(results), "samples processed\n")
  
  return(results)
}
