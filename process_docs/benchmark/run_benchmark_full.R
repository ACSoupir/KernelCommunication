#!/usr/bin/env Rscript

#' Full-Resolution Benchmark Runner
#' 
#' This script runs the complete benchmark with all bandwidths (5-250 by 5)
#' and the full ccRCCExample dataset. WARNING: This may take a long time.

# =============================================================================
# Setup
# =============================================================================

cat("=== KernelCommunication Full-Resolution Benchmark ===\n\n")
cat("WARNING: This benchmark uses the full dataset and all bandwidths.\n")
cat("Expected runtime: 30+ minutes depending on hardware.\n\n")

# Set working directory
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  pkg_dir <- normalizePath(
    file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "..", ".."),
    mustWork = FALSE
  )
  if (dir.exists(pkg_dir)) {
    setwd(pkg_dir)
  }
}

cat("Working directory:", getwd(), "\n\n")

# =============================================================================
# Install/Load Required Packages
# =============================================================================

cat("Step 1: Checking required packages...\n")

required_packages <- c("devtools", "bench", "rmarkdown", "knitr", "ggplot2", "dplyr", "tidyr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("  Installing", pkg, "...\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("  All packages loaded successfully.\n\n")

# =============================================================================
# Load KernelCommunication Package
# =============================================================================

cat("Step 2: Loading KernelCommunication package...\n")

tryCatch({
  devtools::load_all()
  cat("  Package loaded successfully.\n\n")
}, error = function(e) {
  cat("  ERROR loading package:", e$message, "\n")
  quit(status = 1)
})

# =============================================================================
# Verify Dataset
# =============================================================================

cat("Step 3: Verifying example dataset...\n")

data(ccRCCExample, package = "KernelCommunication")
cat("  Dataset loaded:", nrow(ccRCCExample), "cells\n")
cat("  Using FULL dataset (no subsetting)\n\n")

# =============================================================================
# Run Benchmark
# =============================================================================

cat("Step 4: Running full-resolution benchmark...\n")
cat("  Bandwidth range: 5 to 250 (increments of 5)\n")
cat("  Total bandwidths: 50\n")
cat("  Permutations: 100\n")
cat("  Benchmark iterations: 5\n\n")

rmd_file <- "process_docs/benchmark/benchmark_comparison.Rmd"
output_file <- "process_docs/benchmark/benchmark_comparison_full.html"

if (!file.exists(rmd_file)) {
  cat("  ERROR: R Markdown file not found:", rmd_file, "\n")
  quit(status = 1)
}

start_time <- Sys.time()

tryCatch({
  rmarkdown::render(
    input = rmd_file,
    output_file = output_file,
    params = list(
      h_max = 250,
      h_by = 5,          # Full resolution
      perms = 100,       # More permutations for accuracy
      seed = 333,
      bench_iterations = 5,  # More iterations for stable timing
      use_subset = FALSE,    # Use full dataset
      subset_size = 2000
    ),
    envir = new.env()
  )
  cat("  Benchmark completed successfully.\n\n")
}, error = function(e) {
  cat("  ERROR during benchmark:", e$message, "\n")
  quit(status = 1)
})

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

# =============================================================================
# Verify Output
# =============================================================================

cat("Step 5: Verifying output...\n")

if (file.exists(output_file)) {
  cat("  HTML report generated:", normalizePath(output_file), "\n")
  cat("  File size:", file.info(output_file)$size, "bytes\n")
  cat("  Total runtime:", round(runtime, 2), "minutes\n\n")
} else {
  cat("  WARNING: Output file not found\n\n")
}

# =============================================================================
# Summary
# =============================================================================

cat("=== Full-Resolution Benchmark Complete ===\n")
cat("Report location:", normalizePath(output_file), "\n")
cat("Total runtime:", round(runtime, 2), "minutes\n")
