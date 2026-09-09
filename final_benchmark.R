#!/usr/bin/env Rscript
# Final comprehensive benchmark: multiScaleAC vs ground truth functions
# Validates accuracy and performance

library(KernelCommunication)
library(Matrix)

# Load example data
data(ccRCCExample)
coords <- as.matrix(ccRCCExample[, c("x", "y")])
gene1 <- ccRCCExample$CD274
gene2 <- ccRCCExample$PDCD1

cat("FINAL COMPREHENSIVE BENCHMARK\n")
cat("Dataset: ccRCCExample (", nrow(coords), " cells)\n", sep="")

# Test configuration - use subset for faster single-bandwidth tests
h_vec <- c(10, 20, 40, 80)
perms <- 999
seed <- 333

cat("Bandwidths:", paste(h_vec, collapse=", "), "\n")
cat("Permutations:", perms, "\n")
cat("Seed:", seed, "\n\n")

all_tests_pass <- TRUE

# =============================================================================
# Test 1: multiScaleAC vs soupir_moran_multi (Gaussian)
# =============================================================================
cat("=", rep("=", 70), "\n", sep="")
cat("TEST 1: multiScaleAC vs soupir_moran_multi (Gaussian)\n")
cat("=", rep("=", 70), "\n\n", sep="")

cat("Computing soupir_moran_multi (ground truth)...\n")
set.seed(seed)
start_time <- Sys.time()
result_multi <- soupir_moran_multi(
  coords = coords,
  gene1 = gene1,
  gene2 = gene2,
  h_vec = h_vec,
  perms = perms,
  seed = seed
)
time_multi <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat("  Time:", round(time_multi, 3), "seconds\n")

cat("Computing multiScaleAC (gaussian)...\n")
set.seed(seed)
start_time <- Sys.time()
result_multi_new <- multiScaleAC(
  coords = coords,
  gene1 = gene1,
  gene2 = gene2,
  h = h_vec,
  perms = perms,
  seed = seed,
  kernel = "gaussian",
  distance = "euclidean"
)
time_multi_new <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat("  Time:", round(time_multi_new, 3), "seconds\n\n")

# Compare
obs_diff <- max(abs(result_multi$observed - result_multi_new$observed))
perm_diff <- max(abs(result_multi$permutations - result_multi_new$permutations))

cat("Results:\n")
cat("  Observed match:  ", ifelse(obs_diff < 1e-10, "YES", "NO"), 
    "(max diff:", format(obs_diff, scientific=TRUE), ")\n")
cat("  Permutations match: ", ifelse(perm_diff < 1e-10, "YES", "NO"),
    "(max diff:", format(perm_diff, scientific=TRUE), ")\n")
cat("  Speedup:        ", round(time_multi / time_multi_new, 2), "x\n\n")

test1_pass <- (obs_diff < 1e-10) && (perm_diff < 1e-10)
cat("TEST 1:", ifelse(test1_pass, "PASSED", "FAILED"), "\n\n")
all_tests_pass <- all_tests_pass && test1_pass

# =============================================================================
# Test 2: multiScaleAC vs soupir_moran_bv (single bandwidth validation)
# =============================================================================
cat("=", rep("=", 70), "\n", sep="")
cat("TEST 2: multiScaleAC vs soupir_moran_bv (Single Bandwidth)\n")
cat("=", rep("=", 70), "\n\n", sep="")

# Test each bandwidth individually against soupir_moran_bv
h_test <- 20  # Single bandwidth for this test

cat("Testing at h =", h_test, "...\n\n")

cat("Computing soupir_moran_bv (ground truth)...\n")
set.seed(seed)
start_time <- Sys.time()
result_bv <- soupir_moran_bv(
  h = h_test,
  coords = coords,
  gene1 = gene1,
  gene2 = gene2,
  perms = perms,
  eps = 0
)
time_bv <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat("  Time:", round(time_bv, 3), "seconds\n")

cat("Computing multiScaleAC (single bandwidth)...\n")
set.seed(seed)
start_time <- Sys.time()
result_bv_new <- multiScaleAC(
  coords = coords,
  gene1 = gene1,
  gene2 = gene2,
  h = h_test,
  perms = perms,
  seed = seed,
  kernel = "gaussian",
  distance = "euclidean"
)
time_bv_new <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat("  Time:", round(time_bv_new, 3), "seconds\n\n")

# Compare
obs_diff_bv <- abs(result_bv$observed - result_bv_new$observed)
perm_diff_bv <- max(abs(result_bv$permutations - result_bv_new$permutations))

cat("Results:\n")
cat("  Observed match:  ", ifelse(obs_diff_bv < 1e-10, "YES", "NO"), 
    "(max diff:", format(obs_diff_bv, scientific=TRUE), ")\n")
cat("  Permutations match: ", ifelse(perm_diff_bv < 1e-10, "YES", "NO"),
    "(max diff:", format(perm_diff_bv, scientific=TRUE), ")\n")
cat("  Speedup:        ", round(time_bv / time_bv_new, 2), "x\n\n")

test2_pass <- (obs_diff_bv < 1e-10) && (perm_diff_bv < 1e-10)
cat("TEST 2:", ifelse(test2_pass, "PASSED", "FAILED"), "\n\n")
all_tests_pass <- all_tests_pass && test2_pass

# =============================================================================
# Test 3: Manhattan distance validation
# =============================================================================
cat("=", rep("=", 70), "\n", sep="")
cat("TEST 3: Manhattan Distance\n")
cat("=", rep("=", 70), "\n\n", sep="")

cat("Computing multiScaleAC (Manhattan distance)...\n")
set.seed(seed)
start_time <- Sys.time()
result_manhattan <- multiScaleAC(
  coords = coords,
  gene1 = gene1,
  gene2 = gene2,
  h = h_vec,
  perms = perms,
  seed = seed,
  kernel = "gaussian",
  distance = "manhattan"
)
time_manhattan <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat("  Time:", round(time_manhattan, 3), "seconds\n\n")

cat("Results:\n")
cat("  Observed computed: ", all(!is.na(result_manhattan$observed)), "\n")
cat("  Permutations computed: ", all(dim(result_manhattan$permutations) == c(length(h_vec), perms)), "\n")
cat("  Observed range: [", 
    format(min(result_manhattan$observed), digits=4), ", ",
    format(max(result_manhattan$observed), digits=4), "]\n\n")

test3_pass <- all(!is.na(result_manhattan$observed)) && 
              all(dim(result_manhattan$permutations) == c(length(h_vec), perms))
cat("TEST 3:", ifelse(test3_pass, "PASSED", "FAILED"), "\n\n")
all_tests_pass <- all_tests_pass && test3_pass

# =============================================================================
# Test 4: Uniform kernel validation
# =============================================================================
cat("=", rep("=", 70), "\n", sep="")
cat("TEST 4: Uniform Kernel\n")
cat("=", rep("=", 70), "\n\n", sep="")

cat("Computing multiScaleAC (uniform kernel)...\n")
set.seed(seed)
start_time <- Sys.time()
result_uniform <- multiScaleAC(
  coords = coords,
  gene1 = gene1,
  gene2 = gene2,
  h = h_vec,
  perms = perms,
  seed = seed,
  kernel = "uniform",
  distance = "euclidean"
)
time_uniform <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat("  Time:", round(time_uniform, 3), "seconds\n\n")

cat("Results:\n")
cat("  Observed computed: ", all(!is.na(result_uniform$observed)), "\n")
cat("  Permutations computed: ", all(dim(result_uniform$permutations) == c(length(h_vec), perms)), "\n")
cat("  Observed range: [", 
    format(min(result_uniform$observed), digits=4), ", ",
    format(max(result_uniform$observed), digits=4), "]\n\n")

test4_pass <- all(!is.na(result_uniform$observed)) && 
              all(dim(result_uniform$permutations) == c(length(h_vec), perms))
cat("TEST 4:", ifelse(test4_pass, "PASSED", "FAILED"), "\n\n")
all_tests_pass <- all_tests_pass && test4_pass

# =============================================================================
# Test 5: Custom kernel with custom distance
# =============================================================================
cat("=", rep("=", 70), "\n", sep="")
cat("TEST 5: Custom Kernel + Custom Distance\n")
cat("=", rep("=", 70), "\n\n", sep="")

# Define custom functions
minkowski_dist <- function(coords, p = 3) {
  as.matrix(dist(coords, method = "minkowski", p = p))
}

exp_kernel <- function(dist_matrix, h) {
  exp(-dist_matrix / h)
}

cat("Computing multiScaleAC (custom exponential + Minkowski p=3)...\n")
set.seed(seed)
start_time <- Sys.time()
result_custom <- multiScaleAC(
  coords = coords,
  gene1 = gene1,
  gene2 = gene2,
  h = h_vec,
  perms = perms,
  seed = seed,
  kernel = exp_kernel,
  distance = minkowski_dist,
  p = 3
)
time_custom <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat("  Time:", round(time_custom, 3), "seconds\n\n")

cat("Results:\n")
cat("  Observed computed: ", all(!is.na(result_custom$observed)), "\n")
cat("  Permutations computed: ", all(dim(result_custom$permutations) == c(length(h_vec), perms)), "\n")
cat("  Observed range: [", 
    format(min(result_custom$observed), digits=4), ", ",
    format(max(result_custom$observed), digits=4), "]\n\n")

test5_pass <- all(!is.na(result_custom$observed)) && 
              all(dim(result_custom$permutations) == c(length(h_vec), perms))
cat("TEST 5:", ifelse(test5_pass, "PASSED", "FAILED"), "\n\n")
all_tests_pass <- all_tests_pass && test5_pass

# =============================================================================
# Summary
# =============================================================================
cat("=", rep("=", 70), "\n", sep="")
cat("FINAL SUMMARY\n")
cat("=", rep("=", 70), "\n\n", sep="")

cat("Test Results:\n")
cat("  1. multiScaleAC vs soupir_moran_multi: ", ifelse(test1_pass, "PASSED", "FAILED"), "\n")
cat("  2. multiScaleAC vs soupir_moran_bv:    ", ifelse(test2_pass, "PASSED", "FAILED"), "\n")
cat("  3. Manhattan Distance:                 ", ifelse(test3_pass, "PASSED", "FAILED"), "\n")
cat("  4. Uniform Kernel:                     ", ifelse(test4_pass, "PASSED", "FAILED"), "\n")
cat("  5. Custom Kernel + Distance:           ", ifelse(test5_pass, "PASSED", "FAILED"), "\n\n")

cat("OVERALL:", ifelse(all_tests_pass, "ALL TESTS PASSED ✓", "SOME TESTS FAILED ✗"), "\n\n")

# Performance summary
cat("Performance Summary:\n")
cat("  soupir_moran_multi:  ", round(time_multi, 3), "s\n")
cat("  multiScaleAC:        ", round(time_multi_new, 3), "s (", round(time_multi / time_multi_new, 2), "x speedup)\n")
cat("  soupir_moran_bv:     ", round(time_bv, 3), "s (single h)\n")
cat("  multiScaleAC (1 h):  ", round(time_bv_new, 3), "s (", round(time_bv / time_bv_new, 2), "x speedup)\n\n")

# Exit with appropriate code
quit(status = ifelse(all_tests_pass, 0, 1))
