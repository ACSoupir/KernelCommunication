#!/usr/bin/env Rscript
# Benchmark and validate multiScaleAC against soupir_moran_multi (ground truth)
# Tests both gaussian and uniform kernels

library(KernelCommunication)
library(Matrix)

# Load example data
data(ccRCCExample)
coords <- as.matrix(ccRCCExample[, c("x", "y")])
gene1 <- ccRCCExample$CD274
gene2 <- ccRCCExample$PDCD1

cat("Dataset: ccRCCExample (", nrow(coords), " cells)\n", sep="")

# Test configuration
h_vec <- c(10, 20, 40, 80)
perms <- 999
seed <- 333

cat("Bandwidths:", paste(h_vec, collapse=", "), "\n")
cat("Permutations:", perms, "\n")
cat("Seed:", seed, "\n\n")

# =============================================================================
# Test 1: Gaussian Kernel (multiScaleAC vs soupir_moran_multi)
# =============================================================================
cat("=", rep("=", 70), "\n", sep="")
cat("TEST 1: Gaussian Kernel Validation\n")
cat("=", rep("=", 70), "\n\n", sep="")

# Ground truth: soupir_moran_multi (deprecated but still accurate)
cat("Computing ground truth (soupir_moran_multi, gaussian)...\n")
set.seed(seed)
start_time <- Sys.time()
result_gt_gauss <- soupir_moran_multi(
  coords = coords,
  gene1 = gene1,
  gene2 = gene2,
  h_vec = h_vec,
  perms = perms,
  seed = seed
)
time_gt_gauss <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat("  Time:", round(time_gt_gauss, 3), "seconds\n\n")

# New implementation: multiScaleAC with gaussian kernel
cat("Computing multiScaleAC (gaussian)...\n")
set.seed(seed)
start_time <- Sys.time()
result_new_gauss <- multiScaleAC(
  coords = coords,
  gene1 = gene1,
  gene2 = gene2,
  h = h_vec,
  perms = perms,
  seed = seed,
  kernel = "gaussian",
  distance = "euclidean"
)
time_new_gauss <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat("  Time:", round(time_new_gauss, 3), "seconds\n\n")

# Compare gaussian results
cat("Gaussian Kernel Comparison:\n")
obs_diff_gauss <- max(abs(result_gt_gauss$observed - result_new_gauss$observed))
perm_diff_gauss <- max(abs(result_gt_gauss$permutations - result_new_gauss$permutations))

# Compute empirical p-values
pval_gt_gauss <- (rowSums(result_gt_gauss$permutations >= result_gt_gauss$observed) + 1) / (perms + 1)
pval_new_gauss <- (rowSums(result_new_gauss$permutations >= result_new_gauss$observed) + 1) / (perms + 1)
pval_diff_gauss <- max(abs(pval_gt_gauss - pval_new_gauss))

cat("  Observed Moran's I match:     ", ifelse(obs_diff_gauss < 1e-10, "YES", "NO"), 
    "(max diff:", format(obs_diff_gauss, scientific=TRUE), ")\n")
cat("  Permutation values match:     ", ifelse(perm_diff_gauss < 1e-10, "YES", "NO"),
    "(max diff:", format(perm_diff_gauss, scientific=TRUE), ")\n")
cat("  Empirical p-values match:     ", ifelse(pval_diff_gauss < 1e-10, "YES", "NO"),
    "(max diff:", format(pval_diff_gauss, scientific=TRUE), ")\n")
cat("  Speedup:                    ", round(time_gt_gauss / time_new_gauss, 2), "x\n\n")

gauss_pass <- (obs_diff_gauss < 1e-10) && (perm_diff_gauss < 1e-10) && (pval_diff_gauss < 1e-10)
cat("  GAUSSIAN TEST:", ifelse(gauss_pass, "PASSED", "FAILED"), "\n\n")

# =============================================================================
# Test 2: Uniform Kernel (built-in vs custom cross-validation)
# =============================================================================
cat("=", rep("=", 70), "\n", sep="")
cat("TEST 2: Uniform Kernel Validation\n")
cat("=", rep("=", 70), "\n\n", sep="")

# New implementation: multiScaleAC with uniform kernel (built-in)
cat("Computing multiScaleAC (uniform, built-in)...\n")
set.seed(seed)
start_time <- Sys.time()
result_new_unif <- multiScaleAC(
  coords = coords,
  gene1 = gene1,
  gene2 = gene2,
  h = h_vec,
  perms = perms,
  seed = seed,
  kernel = "uniform",
  distance = "euclidean"
)
time_new_unif <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat("  Time:", round(time_new_unif, 3), "seconds\n\n")

# Custom uniform kernel for cross-validation
uniform_kernel <- function(dist_matrix, h) {
  W <- (dist_matrix <= h) * 1.0
  diag(W) <- 0
  W
}

cat("Computing multiScaleAC (uniform, custom kernel)...\n")
set.seed(seed)
start_time <- Sys.time()
result_custom_unif <- multiScaleAC(
  coords = coords,
  gene1 = gene1,
  gene2 = gene2,
  h = h_vec,
  perms = perms,
  seed = seed,
  kernel = uniform_kernel,
  distance = "euclidean"
)
time_custom_unif <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat("  Time:", round(time_custom_unif, 3), "seconds\n\n")

# Compare built-in vs custom uniform
cat("Uniform Kernel Cross-Validation (built-in vs custom):\n")
obs_diff_unif <- max(abs(result_new_unif$observed - result_custom_unif$observed))
perm_diff_unif <- max(abs(result_new_unif$permutations - result_custom_unif$permutations))

# Compute empirical p-values
pval_new_unif <- (rowSums(result_new_unif$permutations >= result_new_unif$observed) + 1) / (perms + 1)
pval_custom_unif <- (rowSums(result_custom_unif$permutations >= result_custom_unif$observed) + 1) / (perms + 1)
pval_diff_unif <- max(abs(pval_new_unif - pval_custom_unif))

cat("  Observed Moran's I match:     ", ifelse(obs_diff_unif < 1e-10, "YES", "NO"), 
    "(max diff:", format(obs_diff_unif, scientific=TRUE), ")\n")
cat("  Permutation values match:     ", ifelse(perm_diff_unif < 1e-10, "YES", "NO"),
    "(max diff:", format(perm_diff_unif, scientific=TRUE), ")\n")
cat("  Empirical p-values match:     ", ifelse(pval_diff_unif < 1e-10, "YES", "NO"),
    "(max diff:", format(pval_diff_unif, scientific=TRUE), ")\n")

unif_pass <- (obs_diff_unif < 1e-10) && (perm_diff_unif < 1e-10) && (pval_diff_unif < 1e-10)
cat("  UNIFORM TEST:", ifelse(unif_pass, "PASSED", "FAILED"), "\n\n")

# =============================================================================
# Test 3: Custom Kernel with Custom Distance (Minkowski p=3)
# =============================================================================
cat("=", rep("=", 70), "\n", sep="")
cat("TEST 3: Custom Kernel + Custom Distance (Minkowski p=3)\n")
cat("=", rep("=", 70), "\n\n", sep="")

# Define custom distance (Minkowski with p=3)
minkowski_dist <- function(coords, p = 3) {
  as.matrix(dist(coords, method = "minkowski", p = p))
}

# Define custom kernel (exponential decay)
exp_kernel <- function(dist_matrix, h) {
  exp(-dist_matrix / h)
}

cat("Computing multiScaleAC (custom exponential + Minkowski p=3)...\n")
set.seed(seed)
start_time <- Sys.time()
# Note: p parameter is passed to distance function, not kernel
result_custom <- multiScaleAC(
  coords = coords,
  gene1 = gene1,
  gene2 = gene2,
  h = h_vec,
  perms = perms,
  seed = seed,
  kernel = exp_kernel,
  distance = minkowski_dist,
  p = 3  # This goes to distance function
)
time_custom <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat("  Time:", round(time_custom, 3), "seconds\n\n")

cat("Custom Kernel/Distance Test:\n")
cat("  Observed Moran's I computed: ", all(!is.na(result_custom$observed)), "\n")
cat("  Permutations computed:       ", all(dim(result_custom$permutations) == c(length(h_vec), perms)), "\n")
cat("  Observed range:              [", 
    format(min(result_custom$observed), digits=4), ", ",
    format(max(result_custom$observed), digits=4), "]\n")

custom_pass <- all(!is.na(result_custom$observed)) && 
               all(dim(result_custom$permutations) == c(length(h_vec), perms))
cat("  CUSTOM TEST:", ifelse(custom_pass, "PASSED", "FAILED"), "\n\n")

# =============================================================================
# Summary
# =============================================================================
cat("=", rep("=", 70), "\n", sep="")
cat("SUMMARY\n")
cat("=", rep("=", 70), "\n\n", sep="")

all_pass <- gauss_pass && unif_pass && custom_pass

cat("Test Results:\n")
cat("  Gaussian Kernel:  ", ifelse(gauss_pass, "PASSED", "FAILED"), "\n")
cat("  Uniform Kernel:   ", ifelse(unif_pass, "PASSED", "FAILED"), "\n")
cat("  Custom K+D:       ", ifelse(custom_pass, "PASSED", "FAILED"), "\n\n")

cat("OVERALL:", ifelse(all_pass, "ALL TESTS PASSED", "SOME TESTS FAILED"), "\n")

# Exit with appropriate code
quit(status = ifelse(all_pass, 0, 1))
