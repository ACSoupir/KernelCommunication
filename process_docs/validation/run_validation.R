#' Validation Test Suite: multiScaleAC() vs soupir_moran_multi()
#'
#' This script validates that multiScaleAC() produces observed values
#' identical to the legacy soupir_moran_multi() function within
#' floating-point precision (1e-14).

# Load development version
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}
devtools::load_all()

# =============================================================================
# Configuration
# =============================================================================

VALIDATION_THRESHOLD <- 1e-14
WARNING_THRESHOLD <- 1e-10

# =============================================================================
# Helper Functions
# =============================================================================

#' Compare observed values between legacy and new implementations
#'
#' @param result_legacy Result from soupir_moran_multi()
#' @param result_new Result from multiScaleAC()
#' @param test_id Test identifier
#' @return List with comparison results
compare_observed <- function(result_legacy, result_new, test_id) {
  obs_legacy <- result_legacy$observed
  obs_new <- result_new$observed

  abs_diff <- abs(obs_legacy - obs_new)
  rel_diff <- abs((obs_legacy - obs_new) / obs_legacy)
  rel_diff[is.infinite(rel_diff)] <- NA

  list(
    test_id = test_id,
    n_h = length(obs_legacy),
    obs_legacy = obs_legacy,
    obs_new = obs_new,
    abs_diff = abs_diff,
    rel_diff = rel_diff,
    max_abs_diff = max(abs_diff),
    mean_abs_diff = mean(abs_diff),
    pass = max(abs_diff) < VALIDATION_THRESHOLD,
    warning = max(abs_diff) >= VALIDATION_THRESHOLD & max(abs_diff) < WARNING_THRESHOLD
  )
}

#' Compare permutation values
#'
#' @param result_legacy Result from soupir_moran_multi()
#' @param result_new Result from multiScaleAC()
#' @param test_id Test identifier
#' @return List with permutation comparison results
compare_permutations <- function(result_legacy, result_new, test_id) {
  perm_legacy <- result_legacy$permutations
  perm_new <- result_new$permutations

  if (!all(dim(perm_legacy) == dim(perm_new))) {
    return(list(
      test_id = test_id,
      error = "Dimension mismatch",
      dim_legacy = dim(perm_legacy),
      dim_new = dim(perm_new)
    ))
  }

  abs_diff <- abs(perm_legacy - perm_new)

  list(
    test_id = test_id,
    n_h = nrow(perm_legacy),
    n_perms = ncol(perm_legacy),
    max_abs_diff = max(abs_diff),
    mean_abs_diff = mean(abs_diff),
    pass = max(abs_diff) < VALIDATION_THRESHOLD
  )
}

#' Generate synthetic test data
#'
#' @param N Number of observations
#' @param seed Random seed
#' @return List with coords, gene1, gene2
generate_synthetic_data <- function(N, seed = 1) {
  set.seed(seed)
  coords <- matrix(runif(2 * N, -50, 50), ncol = 2)
  g1 <- rnorm(N)
  g2 <- rnorm(N)
  list(coords = coords, gene1 = g1, gene2 = g2)
}

#' Run single validation test
#'
#' @param N Dataset size
#' @param h_vec Bandwidth vector
#' @param perms Number of permutations
#' @param test_id Test identifier
#' @param seed Random seed for data generation
#' @param legacy_seed Seed for legacy function
#' @param new_seed Seed for new function
#' @return List with comparison results
run_validation_test <- function(N, h_vec, perms, test_id,
                                 seed = 1, legacy_seed = 333, new_seed = 333) {
  cat(sprintf("\n=== %s (N=%d, h=[%s], perms=%d) ===\n",
              test_id, N, paste(h_vec, collapse = ","), perms))

  # Generate data
  data <- generate_synthetic_data(N, seed)

  # Run legacy function
  cat("Running legacy soupir_moran_multi()... ")
  start <- Sys.time()
  result_legacy <- tryCatch({
    soupir_moran_multi(
      coords = data$coords,
      gene1 = data$gene1,
      gene2 = data$gene2,
      h_vec = h_vec,
      perms = perms,
      seed = legacy_seed
    )
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    return(NULL)
  })
  time_legacy <- difftime(Sys.time(), start, units = "secs")
  cat(sprintf("Done (%.3fs)\n", time_legacy))

  if (is.null(result_legacy)) {
    return(list(test_id = test_id, error = "Legacy function failed"))
  }

  # Run new function
  cat("Running new multiScaleAC()... ")
  start <- Sys.time()
  result_new <- tryCatch({
    multiScaleAC(
      coords = data$coords,
      gene1 = data$gene1,
      gene2 = data$gene2,
      h = h_vec,
      perms = perms,
      seed = new_seed
    )
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    return(NULL)
  })
  time_new <- difftime(Sys.time(), start, units = "secs")
  cat(sprintf("Done (%.3fs)\n", time_new))

  if (is.null(result_new)) {
    return(list(test_id = test_id, error = "New function failed"))
  }

  # Compare observed
  obs_comparison <- compare_observed(result_legacy, result_new, test_id)

  # Compare permutations
  perm_comparison <- compare_permutations(result_legacy, result_new, test_id)

  # Report
  cat(sprintf("Observed: max_abs_diff = %.2e (%s)\n",
              obs_comparison$max_abs_diff,
              ifelse(obs_comparison$pass, "PASS",
                     ifelse(obs_comparison$warning, "WARNING", "FAIL"))))

  if (is.null(perm_comparison$error)) {
    cat(sprintf("Permutations: max_abs_diff = %.2e (%s)\n",
                perm_comparison$max_abs_diff,
                ifelse(perm_comparison$pass, "PASS", "FAIL")))
  }

  list(
    test_id = test_id,
    N = N,
    h_vec = h_vec,
    perms = perms,
    time_legacy = time_legacy,
    time_new = time_new,
    observed = obs_comparison,
    permutations = perm_comparison
  )
}

# =============================================================================
# Test Suite
# =============================================================================

#' Run complete validation test suite
#'
#' @return List with all test results
run_validation_suite <- function() {
  results <- list()

  # Test 1: Small dataset, single bandwidth
  results[[1]] <- run_validation_test(
    N = 100, h_vec = c(10), perms = 10,
    test_id = "T1_small_single_h"
  )

  # Test 2: Small dataset, multiple bandwidths
  results[[2]] <- run_validation_test(
    N = 100, h_vec = c(5, 10, 20), perms = 10,
    test_id = "T2_small_multi_h"
  )

  # Test 3: Medium dataset
  results[[3]] <- run_validation_test(
    N = 500, h_vec = c(10, 20, 40), perms = 50,
    test_id = "T3_medium"
  )

  # Test 4: Large dataset
  results[[4]] <- run_validation_test(
    N = 1000, h_vec = c(20, 40, 80), perms = 100,
    test_id = "T4_large"
  )

  # Test 5: ccRCCExample dataset
  if (requireNamespace("KernelCommunication", quietly = TRUE)) {
    data(ccRCCExample, package = "KernelCommunication")
    cat("\n=== T5_ccRCCExample (N=7122) ===\n")

    # Use subset for speed
    set.seed(42)
    idx <- sample(1:nrow(ccRCCExample), 500)
    coords <- as.matrix(ccRCCExample[idx, c("x", "y")])
    g1 <- ccRCCExample$CD274[idx]
    g2 <- ccRCCExample$PDCD1[idx]

    cat("Running legacy... ")
    start <- Sys.time()
    result_legacy <- soupir_moran_multi(coords, g1, g2, c(0.05, 0.1), perms = 20, seed = 333)
    time_legacy <- difftime(Sys.time(), start, units = "secs")
    cat(sprintf("Done (%.3fs)\n", time_legacy))

    cat("Running new... ")
    start <- Sys.time()
    result_new <- multiScaleAC(coords, g1, g2, c(0.05, 0.1), perms = 20, seed = 333)
    time_new <- difftime(Sys.time(), start, units = "secs")
    cat(sprintf("Done (%.3fs)\n", time_new))

    obs_comp <- compare_observed(result_legacy, result_new, "T5_ccRCCExample")
    perm_comp <- compare_permutations(result_legacy, result_new, "T5_ccRCCExample")

    results[[5]] <- list(
      test_id = "T5_ccRCCExample",
      N = 500,
      observed = obs_comp,
      permutations = perm_comp,
      time_legacy = time_legacy,
      time_new = time_new
    )

    cat(sprintf("Observed: max_abs_diff = %.2e\n", obs_comp$max_abs_diff))
  }

  # Summary
  cat("\n\n========== VALIDATION SUMMARY ==========\n")
  n_pass <- sum(sapply(results, function(r) {
    if (!is.null(r$observed)) r$observed$pass else FALSE
  }))
  n_total <- length(results)
  cat(sprintf("Passed: %d/%d tests\n", n_pass, n_total))

  if (n_pass == n_total) {
    cat("STATUS: ALL TESTS PASSED\n")
  } else {
    cat("STATUS: SOME TESTS FAILED - See details above\n")
  }

  invisible(results)
}

# =============================================================================
# Issue Investigation Helpers
# =============================================================================

#' Extract intermediate values for debugging
#'
#' @param coords Coordinates matrix
#' @param gene1 Gene 1 expression
#' @param gene2 Gene 2 expression
#' @param h Bandwidth
#' @return List with intermediate values
debug_intermediates <- function(coords, gene1, gene2, h) {
  N <- nrow(coords)
  g1_s <- as.vector(scale(gene1))
  g2_s <- as.vector(scale(gene2))
  denom <- sum(g1_s^2)

  # Legacy approach
  dist_sq_legacy <- as.matrix(dist(coords))^2
  K_legacy <- exp(-dist_sq_legacy / (2 * h^2))
  diag(K_legacy) <- 0
  r_sums_legacy <- rowSums(K_legacy)
  x_weighted_legacy <- g1_s / r_sums_legacy
  v_h_legacy <- as.vector(K_legacy %*% x_weighted_legacy)
  obs_legacy <- sum(v_h_legacy * g2_s) / denom

  # New approach (Gaussian kernel, Euclidean distance)
  coords_centered <- coords - colMeans(coords)
  dist_sq_new <- matrix(0, N, N)
  for (i in 1:N) {
    for (j in i:N) {
      dx <- coords_centered[i, 1] - coords_centered[j, 1]
      dy <- coords_centered[i, 2] - coords_centered[j, 2]
      d2 <- dx * dx + dy * dy
      dist_sq_new[i, j] <- d2
      dist_sq_new[j, i] <- d2
    }
  }
  K_new <- exp(-dist_sq_new / (2 * h^2))
  diag(K_new) <- 0
  r_sums_new <- rowSums(K_new)
  x_weighted_new <- g1_s / r_sums_new
  v_h_new <- as.vector(K_new %*% x_weighted_new)
  obs_new <- sum(v_h_new * g2_s) / denom

  list(
    dist_sq_max_diff = max(abs(dist_sq_legacy - dist_sq_new)),
    K_max_diff = max(abs(K_legacy - K_new)),
    r_sums_max_diff = max(abs(r_sums_legacy - r_sums_new)),
    v_h_max_diff = max(abs(v_h_legacy - v_h_new)),
    obs_legacy = obs_legacy,
    obs_new = obs_new,
    obs_diff = abs(obs_legacy - obs_new)
  )
}

# =============================================================================
# Main Execution
# =============================================================================

if (sys.nframe() == 0) {
  cat("Starting validation test suite...\n")
  results <- run_validation_suite()

  # Save results
  save(results, file = "docs/validation/validation_results.RData")
  cat("\nResults saved to docs/validation/validation_results.RData\n")
}
