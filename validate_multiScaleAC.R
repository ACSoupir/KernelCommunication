#!/usr/bin/env Rscript
# Validation script: Compare multiScaleAC vs soupir_moran_bv
# Goal: Verify multiScaleAC produces identical results to soupir_moran_bv (ground truth)

library(KernelCommunication)
library(data.table)

cat("=" , rep("=", 79), "\n", sep="")
cat("VALIDATION REPORT: multiScaleAC vs soupir_moran_bv\n")
cat("=", rep("=", 79), "\n\n", sep="")

# Load example dataset
data(ccRCCExample)
coords <- as.matrix(ccRCCExample[, c("x", "y")])
gene1 <- ccRCCExample$CD274
gene2 <- ccRCCExample$PDCD1

cat("Dataset: ccRCCExample\n")
cat(sprintf("  Number of cells: %d\n", nrow(coords)))
cat(sprintf("  Gene1 range: [%.3f, %.3f]\n", min(gene1), max(gene1)))
cat(sprintf("  Gene2 range: [%.3f, %.3f]\n\n", min(gene2), max(gene2)))

# Test parameters
h_vec <- c(10, 20, 40, 80)
perms <- 999
seed <- 333

cat("Test Configuration:\n")
cat(sprintf("  Bandwidths (h): %s\n", paste(h_vec, collapse=", ")))
cat(sprintf("  Permutations: %d\n", perms))
cat(sprintf("  Random seed: %d\n\n", seed))

# Results storage
results <- list()
summary_table <- data.frame(
  bandwidth = h_vec,
  observed_bv = numeric(length(h_vec)),
  observed_multi = numeric(length(h_vec)),
  difference = numeric(length(h_vec)),
  match = logical(length(h_vec)),
  stringsAsFactors = FALSE
)

cat("Running Tests...\n", rep("-", 80), "\n", sep="")

# Test 1: Compare observed Moran's I values
cat("\n1. OBSERVED MORAN'S I COMPARISON\n")
cat(rep("-", 80), "\n", sep="")

# Run soupir_moran_bv for each bandwidth (ground truth)
cat("\nComputing soupir_moran_bv (ground truth)...\n")
bv_results <- lapply(h_vec, function(h) {
  set.seed(seed)
  soupir_moran_bv(h = h, coords = coords, gene1 = gene1, gene2 = gene2, 
                  perms = perms)
})

# Run multiScaleAC
cat("Computing multiScaleAC...\n")
multi_result <- multiScaleAC(coords = coords, gene1 = gene1, gene2 = gene2, 
                             h = h_vec, perms = perms, seed = seed)

# Compare observed values
cat("\nComparison of observed Moran's I:\n")
cat(sprintf("%-12s %-15s %-15s %-15s %-10s\n", 
            "Bandwidth", "soupir_moran_bv", "multiScaleAC", "Difference", "Match"))
cat(rep("-", 80), "\n", sep="")

all_match <- TRUE
max_diff <- 0

for (i in seq_along(h_vec)) {
  h <- h_vec[i]
  obs_bv <- bv_results[[i]]$observed
  obs_multi <- multi_result$observed[i]
  diff <- abs(obs_bv - obs_multi)
  is_match <- diff < 1e-10  # Tolerance for floating point
  
  summary_table$observed_bv[i] <- obs_bv
  summary_table$observed_multi[i] <- obs_multi
  summary_table$difference[i] <- diff
  summary_table$match[i] <- is_match
  
  cat(sprintf("%-12.1f %-15.10f %-15.10f %-15.2e %-10s\n", 
              h, obs_bv, obs_multi, diff, ifelse(is_match, "YES", "NO")))
  
  if (!is_match) {
    all_match <- FALSE
  }
  max_diff <- max(max_diff, diff)
}

cat(rep("-", 80), "\n", sep="")
cat(sprintf("\nMaximum absolute difference: %.2e\n", max_diff))
cat(sprintf("All observed values match: %s\n\n", ifelse(all_match, "YES ✓", "NO ✗")))

results$observed_match <- all_match
results$max_observed_diff <- max_diff

# Test 2: Compare permutation distributions
cat("\n2. PERMUTATION DISTRIBUTION COMPARISON\n")
cat(rep("-", 80), "\n", sep="")

# Extract permutation matrices
perm_bv_list <- lapply(bv_results, function(x) x$permutations)
perm_multi <- multi_result$permutations

cat("\nComparing permutation distributions (first 10 permutations per bandwidth):\n")
cat(sprintf("%-12s %-10s %-20s %-20s %-10s\n", 
            "Bandwidth", "Perm#", "soupir_moran_bv", "multiScaleAC", "Match"))
cat(rep("-", 80), "\n", sep="")

perm_match <- TRUE
max_perm_diff <- 0

for (i in seq_along(h_vec)) {
  h <- h_vec[i]
  # Compare first 10 permutations
  n_compare <- min(10, perms)
  for (p in 1:n_compare) {
    perm_bv <- perm_bv_list[[i]][p]
    perm_m <- perm_multi[i, p]
    diff <- abs(perm_bv - perm_m)
    is_match <- diff < 1e-10
    
    if (p <= 3 || p == n_compare) {  # Show first 3 and last
      cat(sprintf("%-12.1f %-10d %-20.10f %-20.10f %-10s\n", 
                  h, p, perm_bv, perm_m, ifelse(is_match, "YES", "NO")))
      if (p == 3 && n_compare > 3) {
        cat(sprintf("%-12s %-10s %-20s %-20s %-10s\n", "...", "...", "...", "...", "..."))
      }
    }
    
    if (!is_match) {
      perm_match <- FALSE
    }
    max_perm_diff <- max(max_perm_diff, diff)
  }
}

# Full comparison
full_perm_match <- all(abs(perm_bv_list[[1]] - perm_multi[1,]) < 1e-10)
for (i in 2:length(h_vec)) {
  full_perm_match <- full_perm_match && all(abs(perm_bv_list[[i]] - perm_multi[i,]) < 1e-10)
}

cat(rep("-", 80), "\n", sep="")
cat(sprintf("\nMaximum permutation difference: %.2e\n", max_perm_diff))
cat(sprintf("All permutations match: %s\n\n", ifelse(full_perm_match, "YES ✓", "NO ✗")))

results$permutations_match <- full_perm_match
results$max_perm_diff <- max_perm_diff

# Test 3: Compare p-values (empirical)
cat("\n3. EMPIRICAL P-VALUE COMPARISON\n")
cat(rep("-", 80), "\n", sep="")

cat(sprintf("%-12s %-15s %-15s %-15s %-10s\n", 
            "Bandwidth", "soupir_moran_bv", "multiScaleAC", "Difference", "Match"))
cat(rep("-", 80), "\n", sep="")

pval_match <- TRUE
max_pval_diff <- 0

for (i in seq_along(h_vec)) {
  h <- h_vec[i]
  # Compute empirical p-value from soupir_moran_bv permutations
  obs_bv <- bv_results[[i]]$observed
  perm_bv <- bv_results[[i]]$permutations
  pval_bv <- (sum(perm_bv >= obs_bv) + 1) / (perms + 1)
  
  # Compute empirical p-value from multiScaleAC permutations
  obs <- multi_result$observed[i]
  perm <- multi_result$permutations[i, ]
  pval_multi <- (sum(perm >= obs) + 1) / (perms + 1)
  
  diff <- abs(pval_bv - pval_multi)
  is_match <- diff < 1e-10
  
  cat(sprintf("%-12.1f %-15.10f %-15.10f %-15.2e %-10s\n", 
              h, pval_bv, pval_multi, diff, ifelse(is_match, "YES", "NO")))
  
  if (!is_match) {
    pval_match <- FALSE
  }
  max_pval_diff <- max(max_pval_diff, diff)
}

cat(rep("-", 80), "\n", sep="")
cat(sprintf("\nMaximum p-value difference: %.2e\n", max_pval_diff))
cat(sprintf("All p-values match: %s\n\n", ifelse(pval_match, "YES ✓", "NO ✗")))

results$pvalues_match <- pval_match
results$max_pval_diff <- max_pval_diff

# Test 4: Performance comparison (optional, informative)
cat("\n4. PERFORMANCE COMPARISON\n")
cat(rep("-", 80), "\n", sep="")

cat("\nTiming soupir_moran_bv (sequential, one bandwidth at a time)...\n")
start_bv <- Sys.time()
for (h in h_vec) {
  set.seed(seed)
  soupir_moran_bv(h = h, coords = coords, gene1 = gene1, gene2 = gene2, 
                  perms = perms)
}
time_bv <- as.numeric(difftime(Sys.time(), start_bv, units = "secs"))

cat("Timing multiScaleAC (optimized, all bandwidths)...\n")
start_multi <- Sys.time()
multiScaleAC(coords = coords, gene1 = gene1, gene2 = gene2, 
             h = h_vec, perms = perms, seed = seed)
time_multi <- as.numeric(difftime(Sys.time(), start_multi, units = "secs"))

speedup <- time_bv / time_multi

cat(sprintf("\nsoupir_moran_bv total time: %.3f seconds\n", time_bv))
cat(sprintf("multiScaleAC total time:    %.3f seconds\n", time_multi))
cat(sprintf("Speedup factor:              %.2fx\n\n", speedup))

results$time_bv <- time_bv
results$time_multi <- time_multi
results$speedup <- speedup

# Summary
cat("\n", rep("=", 80), "\n", sep="")
cat("SUMMARY\n")
cat(rep("=", 80), "\n", sep="")

cat(sprintf("\nObserved Moran's I match:     %s (max diff: %.2e)\n", 
            ifelse(results$observed_match, "YES ✓", "NO ✗"), 
            results$max_observed_diff))
cat(sprintf("Permutation values match:     %s (max diff: %.2e)\n", 
            ifelse(results$permutations_match, "YES ✓", "NO ✗"), 
            results$max_perm_diff))
cat(sprintf("Empirical p-values match:     %s (max diff: %.2e)\n", 
            ifelse(results$pvalues_match, "YES ✓", "NO ✗"), 
            results$max_pval_diff))
cat(sprintf("\nPerformance speedup:          %.2fx\n", results$speedup))

overall_pass <- results$observed_match && results$permutations_match && results$pvalues_match

cat("\n", rep("-", 80), "\n", sep="")
if (overall_pass) {
  cat("✓ VALIDATION PASSED: multiScaleAC produces identical results to soupir_moran_bv\n")
} else {
  cat("✗ VALIDATION FAILED: Discrepancies detected between functions\n")
}
cat(rep("-", 80), "\n", sep="")

# Save detailed results
cat("\n\nDetailed Comparison Table:\n")
cat(rep("-", 80), "\n", sep="")
print(summary_table)

# Save results to file
output_file <- "validation_results.txt"
sink(output_file)
cat("VALIDATION RESULTS: multiScaleAC vs soupir_moran_bv\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Dataset: ccRCCExample (", nrow(coords), " cells)\n", sep="")
cat("Bandwidths: ", paste(h_vec, collapse=", "), "\n")
cat("Permutations: ", perms, "\n")
cat("Seed: ", seed, "\n\n")
cat("RESULTS:\n")
cat("--------\n")
cat(sprintf("Observed Moran's I match:     %s (max diff: %.2e)\n", 
            ifelse(results$observed_match, "YES", "NO"), 
            results$max_observed_diff))
cat(sprintf("Permutation values match:     %s (max diff: %.2e)\n", 
            ifelse(results$permutations_match, "YES", "NO"), 
            results$max_perm_diff))
cat(sprintf("Empirical p-values match:     %s (max diff: %.2e)\n", 
            ifelse(results$pvalues_match, "YES", "NO"), 
            results$max_pval_diff))
cat(sprintf("\nPerformance:\n"))
cat(sprintf("  soupir_moran_bv time:  %.3f seconds\n", results$time_bv))
cat(sprintf("  multiScaleAC time:     %.3f seconds\n", results$time_multi))
cat(sprintf("  Speedup:               %.2fx\n", results$speedup))
cat("\n\nDETAILED COMPARISON:\n")
cat("-------------------\n")
print(summary_table)
cat("\n\nOVERALL: ", ifelse(overall_pass, "PASSED", "FAILED"), "\n", sep="")
sink()

cat(sprintf("\nDetailed results saved to: %s\n", output_file))

# Return results invisibly
invisible(list(
  summary = summary_table,
  overall_pass = overall_pass,
  results = results
))
