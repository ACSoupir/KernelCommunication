#!/usr/bin/env Rscript

#' Benchmark Runner Script
#' 
#' This script runs the benchmark comparison between soupir_moran_multi() 
#' and multiScaleAC() and generates an HTML report.

# =============================================================================
# Setup
# =============================================================================

cat("=== KernelCommunication Benchmark Runner ===\n\n")

# Set working directory to package root (use current directory if not in RStudio)
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
  cat("  Attempting to install package...\n")
  install.packages(".", repos = NULL, type = "source")
  library(KernelCommunication)
  cat("  Package installed and loaded.\n\n")
})

# =============================================================================
# Verify Dataset
# =============================================================================

cat("Step 3: Verifying example dataset...\n")

data(ccRCCExample, package = "KernelCommunication")
cat("  Dataset loaded:", nrow(ccRCCExample), "cells\n\n")

# =============================================================================
# Run Benchmark
# =============================================================================

cat("Step 4: Running benchmark...\n")

rmd_file <- file.path(getwd(), "process_docs/benchmark/benchmark_comparison.Rmd")
output_file <- file.path(getwd(), "process_docs/benchmark/benchmark_comparison.html")

cat("  Rmd file:", rmd_file, "\n")
cat("  Output file:", output_file, "\n")

if (!file.exists(rmd_file)) {
  cat("  ERROR: R Markdown file not found:", rmd_file, "\n")
  quit(status = 1)
}

tryCatch({
  rmarkdown::render(
    input = rmd_file,
    output_file = output_file,
    params = list(
      h_max = 250,
      h_by = 25,        # Fewer bandwidths for speed (10 instead of 51)
      perms = 50,       # Reduced permutations
      seed = 333,
      bench_iterations = 3,
      use_subset = TRUE,
      subset_size = 2000
    ),
    envir = new.env()
  )
  cat("  Benchmark completed successfully.\n\n")
}, error = function(e) {
  cat("  ERROR during benchmark:", e$message, "\n")
  quit(status = 1)
})

# =============================================================================
# Verify Output
# =============================================================================

cat("Step 5: Verifying output...\n")

if (file.exists(output_file)) {
  cat("  HTML report generated:", normalizePath(output_file), "\n")
  cat("  File size:", file.info(output_file)$size, "bytes\n\n")
} else {
  cat("  WARNING: Output file not found\n\n")
}

# =============================================================================
# Summary
# =============================================================================

cat("=== Benchmark Complete ===\n")
cat("Report location:", normalizePath(output_file), "\n")
