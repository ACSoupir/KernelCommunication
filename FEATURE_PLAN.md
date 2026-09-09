# KernelCommunication Feature Plan

**Last updated**: 2026-07-01

---

## Core Design Principles

### 1. **Modularity First** ⭐
- `multiScaleAC()` returns **raw computation only**: `h`, `observed`, `permutations`
- Separate function `compute_significance()` adds p-values, z-scores, CI
- Users can skip significance if they just want raw results
- Each function has single responsibility → easier to test, maintain

### 2. **S4 Objects for Structure** ⭐
- Enforce data structure via slots with validation
- Prevent users from breaking downstream functions
- Enable formal methods (`plot()`, `summary()`, `export_mxfda()`)
- More professional, consistent with spatial statistics packages

### 3. **Single Object, Smart Routing** ⭐
- Single `kernelcomm_result` S4 class for all cases
- `multiScaleAC()` detects if input is single sample or batch
- Returns single object or list of objects accordingly
- Simpler API, less cognitive load for users

### 4. **Optional Slots (NULL when not computed)** ⭐
- Significance slots are `NULL` until `compute_significance()` is called
- Not filled with NA (which suggests "unfinished" or "error")
- Clear distinction: slot exists + has value = computed; slot exists + NULL = not computed yet

### 5. **Memory-Efficient Parallelization** ⭐
- Use `parallel::mclapply()` with `mc.preschedule = TRUE` for Unix
- Process one sample at a time, clear environment between iterations
- For Windows: `parallel::parLapply()` with explicit cluster management
- Document memory requirements: ~200MB per core for typical dataset

---

## Revised Phase 1: Core Workflow

### 1.1 S4 Class Definition
**File:** `R/classes.R` (new file)

**Design decision**: Single `kernelcomm_result` class for all cases. No separate batch class.

```r
#' KernelCommunication Result Object
#'
#' @slot h Numeric vector of bandwidths
#' @slot observed Numeric vector of observed Moran's I
#' @slot permutations Numeric matrix [n_h × n_perms] of permutation results
#' @slot coords Original coordinates (for plotting)
#' @slot gene1_name Character. Name of gene 1
#' @slot gene2_name Character. Name of gene 2
#' @slot sample_id Character. Sample identifier
#' @slot seed Random seed used
#' @slot n_cells Integer. Number of cells
#' @slot call Function call that created this object
#' @slot p_values Numeric vector of empirical p-values (NULL until computed)
#' @slot z_scores Numeric vector of z-scores (NULL until computed)
#' @slot perm_mean Numeric vector, mean of permutations (NULL until computed)
#' @slot perm_sd Numeric vector, SD of permutations (NULL until computed)
#' @slot ci_lower Numeric vector, lower CI bound (NULL until computed)
#' @slot ci_upper Numeric vector, upper CI bound (NULL until computed)
#' @slot significant Logical vector, p < alpha (NULL until computed)
#' @slot p_adjusted Numeric vector, adjusted p-values (NULL until computed)
#'
#' @export
setClass("kernelcomm_result",
  slots = list(
    # Required slots (always present)
    h = "numeric",
    observed = "numeric",
    permutations = "matrix",
    coords = "matrix",
    gene1_name = "character",
    gene2_name = "character",
    sample_id = "character",
    seed = "numeric",
    n_cells = "integer",
    call = "language",
    
    # Optional slots (NULL until compute_significance() is called)
    p_values = "optional",
    z_scores = "optional",
    perm_mean = "optional",
    perm_sd = "optional",
    ci_lower = "optional",
    ci_upper = "optional",
    significant = "optional",
    p_adjusted = "optional"
  )
)

#' Validity check for kernelcomm_result
setValidity("kernelcomm_result", function(object) {
  # Required slot validation
  if (length(object@h) != length(object@observed))
    return("'h' and 'observed' must have same length")
  if (!all(dim(object@permutations) == c(length(object@h), ncol(object@permutations))))
    return("'permutations' rows must match length of 'h'")
  if (nrow(object@coords) != object@n_cells)
    return("coords rows should match n_cells")
  
  # Optional slot validation: if one significance slot is non-NULL, all must be
  sig_slots <- c("p_values", "z_scores", "perm_mean", "perm_sd", 
                 "ci_lower", "ci_upper", "significant")
  sig_non_null <- sapply(sig_slots, function(s) !is.null(slot(object, s)))
  
  if (any(sig_non_null) && !all(sig_non_null))
    return("All significance slots must be NULL or all must be non-NULL")
  
  # If significance is computed, check lengths match
  if (!is.null(object@p_values) && length(object@p_values) != length(object@h))
    return("Significance slots must match length of h")
  
  TRUE
})
```

**Why single class?**
- Simpler API: users don't need to learn two classes
- `multiScaleAC()` can route internally based on input
- Batch processing returns `list` of `kernelcomm_result` objects
- Easier to maintain, document, and test

---

### 1.2 Refactored multiScaleAC()
**File:** `R/multi_scale_ac.R` (modified)

```r
#' Compute Multi-Scale Bivariate Moran's I
#'
#' Returns raw computation results. Use compute_significance() to add
#' p-values, z-scores, and confidence intervals.
#'
#' @return Object of class \code{kernelcomm_result}
#' @export
multiScaleAC <- function(coords, gene1, gene2, h_vec, perms = 999, 
                          seed = 333, kernel = "gaussian", distance = "euclidean", ...) {
  
  # ... existing computation ...
  
  # Create S4 object with raw results only
  result <- new("kernelcomm_result",
    h = h_vec,
    observed = observed_moran,
    permutations = perm_matrix,  # [n_h × n_perms]
    coords = as.matrix(coords),
    gene1_name = deparse(substitute(gene1)),
    gene2_name = deparse(substitute(gene2)),
    sample_id = "sample1",  # or from metadata
    seed = seed,
    n_cells = nrow(coords),
    call = match.call()
  )
  
  return(result)
}
```

**Key change:** No significance calculations here. Just raw output.

---

### 1.3 compute_significance() Function
**File:** `R/significance.R` (new file)

```r
#' Compute Significance Statistics from Permutation Results
#'
#' Adds empirical p-values, z-scores, and confidence intervals to
#' a kernelcomm_result object.
#'
#' @param result Object of class \code{kernelcomm_result}
#' @param alternative Character. One of "two.sided", "greater", "less"
#' @param alpha Significance level for CI (default: 0.05)
#'
#' @return Same object with added slots:
#'   \item{p_values}{Empirical p-values}
#'   \item{z_scores}{Standardized effect sizes}
#'   \item{perm_mean}{Mean of permutations}
#'   \item{perm_sd}{SD of permutations}
#'   \item{ci_lower}{Lower CI bound}
#'   \item{ci_upper}{Upper CI bound}
#'   \item{significant}{Logical. p < alpha}
#'
#' @export
#'
#' @examples
#' result <- multiScaleAC(coords, gene1, gene2, h = c(10, 20, 40))
#' result <- compute_significance(result)
#' summary(result)  # Now shows p-values
compute_significance <- function(result, 
                                  alternative = "two.sided",
                                  alpha = 0.05) {
  
  # Extract permutation matrix
  perms <- result@permutations
  observed <- result@observed
  
  # Compute permutation statistics
  perm_mean <- rowMeans(perms)
  perm_sd <- apply(perms, 1, sd)
  
  # Empirical p-values (with correction for finite permutations)
  # p = (r + 1) / (n + 1) where r = rank of observed
  p_values <- sapply(1:length(observed), function(i) {
    r <- sum(perms[i, ] >= observed[i])
    (r + 1) / (ncol(perms) + 1)
  })
  
  # Adjust for one-sided tests
  if (alternative == "less") {
    p_values <- 1 - p_values
  }
  
  # Z-scores
  z_scores <- (observed - perm_mean) / perm_sd
  
  # Confidence intervals from permutation distribution
  ci_lower <- apply(perms, 1, quantile, probs = alpha / 2)
  ci_upper <- apply(perms, 1, quantile, probs = 1 - alpha / 2)
  
  # Add slots to object
  result@p_values <- p_values
  result@z_scores <- z_scores
  result@perm_mean <- perm_mean
  result@perm_sd <- perm_sd
  result@ci_lower <- ci_lower
  result@ci_upper <- ci_upper
  result@significant <- p_values < alpha
  
  return(result)
}
```

**Benefits of modularity:**
- Can compute significance with different alpha levels without recomputing
- Can try different alternative hypotheses post-hoc
- Easier to test: test computation and significance separately
- Users can skip significance if doing custom analysis

---

### 1.4 export_mxfda() Function
**File:** `R/export_mxfda.R` (new file)

```r
#' Export to mxfda Format
#'
#' Converts kernelcomm_result to data.frame compatible with
#' mxfda::add_summary_function().
#'
#' @param result Object of class \code{kernelcomm_result} (must have significance computed)
#' @param sample_id Character. Sample identifier (overrides object slot)
#' @param interaction_name Character. Gene pair name (e.g., "CD274_PDCD1")
#'
#' @return data.frame with columns: sample_id, r, sumfun, csr, fundiff, 
#'         p_value, z_score, ci_lower, ci_upper
#' @export
#'
#' @examples
#' result <- multiScaleAC(coords, gene1, gene2, h = c(10, 20, 40))
#' result <- compute_significance(result)
#' mxfda_data <- export_mxfda(result, sample_id = "RCC4_17")
export_mxfda <- function(result, 
                         sample_id = NULL,
                         interaction_name = NULL) {
  
  # Check that significance was computed
  if (!"p_values" %in% slotNames(result) || is.null(result@p_values)) {
    stop("Must call compute_significance() before export_mxfda()")
  }
  
  sample_id <- sample_id %||% result@sample_id
  
  # Create interaction name if not provided
  if (is.null(interaction_name)) {
    interaction_name <- paste0(result@gene1_name, "_", result@gene2_name)
  }
  
  # Build data.frame in mxfda format
  data.frame(
    sample_id = rep(sample_id, length(result@h)),
    interaction_name = rep(interaction_name, length(result@h)),
    r = result@h,
    sumfun = result@observed,
    csr = result@perm_mean,
    fundiff = result@observed - result@perm_mean,
    p_value = result@p_values,
    z_score = result@z_scores,
    ci_lower = result@ci_lower,
    ci_upper = result@ci_upper,
    stringsAsFactors = FALSE
  )
}
```

---

### 1.5 S3 Methods
**Files:** `R/methods.R` (new file)

```r
#' Print method for kernelcomm_result
#' @export
print.kernelcomm_result <- function(x, ...) {
  cat("KernelCommunication Result\n")
  cat("==========================\n")
  cat("Sample:", x@sample_id, "\n")
  cat("Gene pair:", x@gene1_name, "×", x@gene2_name, "\n")
  cat("Bandwidths:", length(x@h), "values from", min(x@h), "to", max(x@h), "\n")
  cat("Permutations:", ncol(x@permutations), "\n")
  cat("Cells:", x@n_cells, "\n")
  
  if (!is.null(x@p_values)) {
    cat("Significant bandwidths (α=0.05):", sum(x@p_values < 0.05), "/", length(x@h), "\n")
  } else {
    cat("Significance: Not computed (use compute_significance())\n")
  }
}

#' Summary method
#' @export
summary.kernelcomm_result <- function(object, alpha = 0.05, ...) {
  if (is.null(object@p_values)) {
    stop("Must call compute_significance() first")
  }
  
  data.frame(
    h = object@h,
    observed = object@observed,
    perm_mean = object@perm_mean,
    perm_sd = object@perm_sd,
    z_score = object@z_scores,
    p_value = object@p_values,
    significant = object@p_values < alpha,
    ci_lower = object@ci_lower,
    ci_upper = object@ci_upper
  )
}

#' Plot method
#' @export
plot.kernelcomm_result <- function(x, y, type = "observed", ...) {
  # Dispatch to appropriate plot function
  plot_scale_curve(x, type = type, ...)
}
```

---

## Additional Features to Strengthen the Package

### 2.1 Bandwidth Selection Helpers ⭐⭐⭐

**Problem:** Users don't know what h values to choose.

**Solution:** Helper functions to guide bandwidth selection.

```r
#' Suggest Bandwidth Values
#'
#' Recommends bandwidth values based on nearest neighbor distances.
#'
#' @param coords Coordinate matrix
#' @param n_bandwidths Number of bandwidth values to suggest
#' @param range Character. One of "local", "medium", "global"
#'
#' @return Numeric vector of suggested h values
#' @export
suggest_bandwidths <- function(coords, n_bandwidths = 5, range = "medium") {
  # Compute nearest neighbor distances
  dists <- as.matrix(dist(coords))
  diag(dists) <- Inf
  nn1 <- apply(dists, 1, min)  # 1st nearest neighbor
  
  # Choose range based on user preference
  if (range == "local") {
    quantiles <- seq(0.1, 0.5, length.out = n_bandwidths)
  } else if (range == "medium") {
    quantiles <- seq(0.25, 0.75, length.out = n_bandwidths)
  } else {  # global
    quantiles <- seq(0.5, 0.9, length.out = n_bandwidths)
  }
  
  quantile(nn1, probs = quantiles)
}

#' Plot Neighborhood Size Distribution
#'
#' Shows how many neighbors each bandwidth captures.
#'
#' @param coords Coordinates
#' @param h_vec Bandwidths to evaluate
#'
#' @return ggplot object
#' @export
plot_neighborhood_stats <- function(coords, h_vec) {
  # For each h, compute: mean neighbors, min, max, % isolated
  # Plot as line + ribbon
}
```

**Why this matters:**
- Prevents arbitrary bandwidth selection
- Makes analysis reproducible
- Helps users understand spatial scale of their data

---

### 2.2 Multiple Testing Correction ⭐⭐⭐

**Problem:** When testing 1100 gene pairs, false positives accumulate.

**Solution:** Built-in FDR control.

```r
#' Adjust P-values for Multiple Testing
#'
#' Applies multiple testing correction across bandwidths and/or gene pairs.
#'
#' @param result kernelcomm_result or kernelcomm_batch object
#' @param method Character. One of "BH", "BY", "bonferroni", "holm", "none"
#' @param across Character. One of "bandwidths", "genes", "both"
#'
#' @return Object with added slot: p_adjusted
#' @export
adjust_pvalues <- function(result, method = "BH", across = "bandwidths") {
  if (across == "bandwidths") {
    # Adjust within this sample (multiple h values)
    result@p_adjusted <- p.adjust(result@p_values, method = method)
  } else if (across == "genes") {
    # Requires batch object
    stop("Use adjust_pvalues_batch() for gene-level correction")
  }
  result
}
```

---

### 2.3 Input Validation & QC ⭐⭐

**File:** `R/qc.R` (new file)

```r
#' Validate Input Data
#'
#' Checks for common issues before computation.
#'
#' @param coords Coordinates
#' @param gene1 Gene 1 expression
#' @param gene2 Gene 2 expression
#'
#' @return List of warnings/errors
#' @export
validate_inputs <- function(coords, gene1, gene2) {
  issues <- list(errors = character(), warnings = character())
  
  # Length checks
  if (nrow(coords) != length(gene1)) {
    issues$errors <- c(issues$errors, "coords and gene1 length mismatch")
  }
  
  # Duplicate coordinates
  if (any(duplicated(coords))) {
    issues$warnings <- c(issues$warnings, 
                         "Duplicate coordinates found. Adding small jitter.")
  }
  
  # Zero inflation
  zero_frac1 <- mean(gene1 == 0)
  zero_frac2 <- mean(gene2 == 0)
  if (zero_frac1 > 0.9) {
    issues$warnings <- c(issues$warnings,
                         paste0("gene1 is ", round(zero_frac1 * 100), 
                                "% zeros. Consider filtering."))
  }
  
  # Expression variance
  if (sd(gene1) == 0) {
    issues$errors <- c(issues$errors, "gene1 has zero variance")
  }
  
  return(issues)
}

#' Compute Neighborhood Quality Metrics
#'
#' @param coords Coordinates
#' @param h Bandwidth
#'
#' @return List with neighborhood statistics
#' @export
qc_neighborhoods <- function(coords, h) {
  dist_mat <- as.matrix(dist(coords))
  W <- (dist_mat <= h) * 1
  diag(W) <- 0
  
  n_neighbors <- rowSums(W)
  
  list(
    mean_neighbors = mean(n_neighbors),
    sd_neighbors = sd(n_neighbors),
    min_neighbors = min(n_neighbors),
    max_neighbors = max(n_neighbors),
    isolated_cells = sum(n_neighbors == 0),
    connected_components <- get_connected_components(W)
  )
}
```

**Integration:** Call `validate_inputs()` at start of `multiScaleAC()`.

---

### 2.4 Diagnostic Plots ⭐⭐

**File:** `R/diagnostics.R` (new file)

```r
#' Plot Permutation Distribution
#'
#' Shows histogram of permutation values with observed statistic marked.
#'
#' @param result kernelcomm_result object
#' @param h_index Which bandwidth to show
#'
#' @return ggplot object
#' @export
plot_permutation_distribution <- function(result, h_index = 1) {
  # Histogram of permutations
  # Vertical line at observed value
  # Show p-value in corner
}

#' QQ Plot of Permutation Distribution
#'
#' @param result kernelcomm_result object
#'
#' @return ggplot object
#' @export
plot_perm_qq <- function(result) {
  # Check if permutations are approximately normal
  # Useful for validating z-score assumption
}

#' Sensitivity Analysis
#'
#' Run with different random seeds to check stability.
#'
#' @param coords Coordinates
#' @param ... Other args to multiScaleAC
#' @param n_seeds Number of seeds to try
#'
#' @return List of results
#' @export
check_permutation_stability <- function(coords, ..., n_seeds = 5) {
  seeds <- sample(10000, n_seeds)
  results <- lapply(seeds, function(s) {
    multiScaleAC(coords, ..., seed = s)
  })
  
  # Compare variance across seeds
  # Should be small if perms converged
}
```

---

### 2.5 Result Manipulation Methods ⭐

```r
#' Subset kernelcomm_result
#' @export
`[.kernelcomm_result` <- function(x, i) {
  # Subset to specific bandwidths
  x@h <- x@h[i]
  x@observed <- x@observed[i]
  x@permutations <- x@permutations[i, , drop = FALSE]
  if (!is.null(x@p_values)) x@p_values <- x@p_values[i]
  # ... etc for all slots
  x
}

#' Combine batch results
#' @export
rbind.kernelcomm_batch <- function(...) {
  # Merge multiple batch results
  # Useful for incremental processing
}

#' Extract specific bandwidth
#' @export
slice_result <- function(result, h_value) {
  # Get results at specific h
  idx <- which.min(abs(result@h - h_value))
  result[idx]
}
```

---

### 2.6 Alternative Significance Methods ⭐

```r
#' Analytical P-value Approximation
#'
#' Fast approximation using normal approximation to permutation distribution.
#' Useful for screening before full permutation test.
#'
#' @param result kernelcomm_result
#'
#' @return Object with p_values_analytical slot
#' @export
compute_analytical_pvalues <- function(result) {
  # Z-test approximation
  z <- (result@observed - result@perm_mean) / result@perm_sd
  p <- 2 * pnorm(-abs(z))  # two-sided
  
  result@p_values_analytical <- p
  result
}

#' Tail Approximation for Extreme P-values
#'
#' When observed is more extreme than all permutations,
#' fit generalized Pareto to tail.
#'
#' @param result kernelcomm_result
#'
#' @return Object with p_values_tail slot
#' @export
compute_tail_pvalues <- function(result) {
  # Use extreme value theory
  # For when p < 1/n_perms
}
```

---

### 2.6 Parallel Processing Strategy ⭐⭐⭐

**Design Decision Document**

#### Option 1: `parallel::mclapply()` (Unix/Mac only)

**Benefits:**
- Simple API: `mclapply(data_list, process_one, mc.cores = 4)`
- Automatic fork-based parallelization
- Good performance, minimal overhead
- Built-in to R (no new dependencies)

**Downsides:**
- **Doesn't work on Windows** (falls back to sequential)
- **Memory duplication**: Each core gets copy of parent environment
- Can't easily control when objects are copied

**Memory concern:** For 160 samples × 1100 gene pairs:
- If parent environment has 500MB data
- Using 8 cores → potentially 4GB additional memory
- **Mitigation**: Load data inside worker function, not in parent

```r
# GOOD: minimal memory duplication
results <- mclapply(sample_ids, function(id) {
  # Load data inside worker (not copied from parent)
  data <- load_sample(id)
  multiScaleAC(data$coords, data$gene1, data$gene2, ...)
}, mc.cores = 4, mc.preschedule = TRUE)

# BAD: large memory duplication
all_data <- lapply(sample_ids, load_sample)  # Loaded in parent
results <- mclapply(all_data, function(d) {  # Copied to each core
  multiScaleAC(d$coords, d$gene1, d$gene2, ...)
}, mc.cores = 4)
```

---

#### Option 2: `future` + `furrr` Package

**Benefits:**
- Cross-platform (Windows, Mac, Linux)
- More control over memory: can use "shared memory" strategies
- Modern, actively maintained
- Easy to switch between sequential/parallel: `plan(multisession)` vs `plan(sequential)`

**Downsides:**
- **New dependency** (future, furrr packages)
- Slightly more complex API
- Can be slower for very short tasks (overhead)

```r
library(furrr)
plan(multisession, workers = 4)

results <- future_map(sample_ids, function(id) {
  data <- load_sample(id)
  multiScaleAC(data$coords, data$gene1, data$gene2, ...)
}, .options = furrr_options(seed = TRUE))
```

---

#### Option 3: Manual Cluster with `parLapply()` (Windows-compatible)

**Benefits:**
- Cross-platform
- Explicit control over what gets exported to workers
- Can minimize memory by only exporting what's needed

**Downsides:**
- More boilerplate code
- Must manually manage cluster lifecycle
- Easy to forget to stop cluster (resource leak)

```r
cl <- parallel::makeCluster(4)

# Only export the function, not the data
parallel::clusterExport(cl, c("multiScaleAC"), envir = environment())

results <- parallel::parLapply(cl, sample_ids, function(id) {
  # Load data inside worker
  data <- load_sample(id)
  multiScaleAC(data$coords, data$gene1, data$gene2, ...)
})

parallel::stopCluster(cl)
```

---

#### **Recommendation for This Application**

**Use `parallel::mclapply()` with careful memory management:**

```r
multiScaleAC_batch <- function(data_list, h_vec, perms = 999, 
                                n_cores = 1, progress = TRUE, ...) {
  
  # Validate inputs
  if (n_cores < 1) n_cores <- 1
  
  # Process function (loads data if needed, computes)
  process_one <- function(i) {
    sample_data <- data_list[[i]]
    
    # Compute
    result <- multiScaleAC(
      coords = sample_data$coords,
      gene1 = sample_data$gene1,
      gene2 = sample_data$gene2,
      h_vec = h_vec,
      perms = perms,
      seed = sample_data$seed %||% (333 + i),
      ...
    )
    
    # Force garbage collection after each sample
    gc(verbose = FALSE)
    
    result
  }
  
  # Progress bar wrapper
  if (progress && n_cores == 1) {
    pb <- txtProgressBar(0, length(data_list), style = 3)
    results <- lapply(seq_along(data_list), function(i) {
      result <- process_one(i)
      setTxtProgressBar(pb, i)
      result
    })
    close(pb)
  } else if (n_cores > 1) {
    if (.Platform$OS.type == "unix") {
      # Unix: use mclapply
      results <- parallel::mclapply(
        seq_along(data_list),
        process_one,
        mc.cores = n_cores,
        mc.preschedule = TRUE,  # Better for uneven workloads
        mc.set.seed = TRUE
      )
    } else {
      # Windows: warn and fall back to sequential
      warning("Parallel processing not supported on Windows. Running sequentially.")
      results <- lapply(seq_along(data_list), process_one)
    }
  } else {
    # Sequential
    results <- lapply(seq_along(data_list), process_one)
  }
  
  # Return list of kernelcomm_result objects
  results
}
```

**Memory management best practices:**

1. **Load data inside workers** (don't pre-load in parent)
2. **Use `mc.preschedule = TRUE`** for uneven workloads
3. **Call `gc()` after each sample** to free memory
4. **Document memory requirements**: ~200-500MB per core for typical dataset
5. **For very large batches**: Process in chunks, save intermediate results

```r
# For 176,000 computations (160 samples × 1100 genes):
# Recommended approach: chunking

chunk_size <- 100  # gene pairs per chunk
chunks <- split(gene_pairs, ceiling(seq_along(gene_pairs) / chunk_size))

for (i in seq_along(chunks)) {
  cat("Processing chunk", i, "of", length(chunks), "\n")
  
  chunk_results <- multiScaleAC_batch(
    data_list = chunks[[i]],
    n_cores = 4
  )
  
  # Save chunk to disk
  saveRDS(chunk_results, sprintf("results_chunk_%03d.rds", i))
  
  # Clear memory
  rm(chunk_results)
  gc()
}

# Combine all chunks
all_results <- list.files(pattern = "results_chunk_.*\\.rds") %>%
  lapply(readRDS) %>%
  do.call(c, .)
```

---

### 2.7 Caching & Resumption ⭐⭐

**File:** `R/caching.R` (new file)

```r
#' Save Result to Disk
#'
#' @param result kernelcomm_result or kernelcomm_batch
#' @param file File path
#' @export
save_result <- function(result, file) {
  saveRDS(result, file)
}

#' Load Result from Disk
#' @param file File path
#' @export
load_result <- function(file) {
  readRDS(file)
}

#' Run Batch with Checkpointing
#'
#' Saves progress after each sample. Can resume if interrupted.
#'
#' @param data_list List of samples
#' @param checkpoint_file Path to save checkpoint
#' @param resume Logical. Resume from checkpoint if exists?
#'
#' @return kernelcomm_batch object
#' @export
multiScaleAC_batch_checkpoint <- function(data_list, 
                                           checkpoint_file = "checkpoint.rds",
                                           resume = TRUE) {
  # Check if checkpoint exists
  if (resume && file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    cat("Resuming from checkpoint. Completed:", length(checkpoint), "\n")
  } else {
    checkpoint <- list()
  }
  
  # Process remaining samples
  for (i in seq_along(data_list)) {
    if (i <= length(checkpoint)) next  # already done
    
    result <- multiScaleAC(...)
    checkpoint[[i]] <- result
    
    # Save checkpoint
    saveRDS(checkpoint, checkpoint_file)
    cat("Completed sample", i, "/", length(data_list), "\n")
  }
  
  # Combine into batch object
  as.kernelcomm_batch(checkpoint)
}
```

**Why this matters:**
- For 176,000 computations, interruptions will happen
- Don't lose 23 hours of work
- Can run overnight safely

---

### 2.8 Export Formats Beyond mxfda ⭐

```r
#' Export to CSV/Parquet
#'
#' @param result kernelcomm_result or kernelcomm_batch
#' @param file Output file path
#' @param format Character. "csv", "parquet", "tsv"
#' @export
export_results <- function(result, file, format = "parquet") {
  df <- as.data.frame(result)
  
  if (format == "parquet") {
    arrow::write_parquet(df, file, compression = "zstd")
  } else if (format == "csv") {
    readr::write_csv(df, file)
  }
}

#' Export to SpatialDM-compatible Format
#'
#' @param result kernelcomm_result
#' @param file Output file
#' @export
export_spatialdm <- function(result, file) {
  # Format matching spatialDM output
  # For comparison/integration
}

#' Export to Giotto Object
#'
#' @param result kernelcomm_result
#' @param giotto_obj Giotto object to add results to
#' @export
export_giotto <- function(result, giotto_obj) {
  # Add as spatial network or gene pair scores
}
```

---

## Revised Class Structure

### Single Sample Result
```r
setClass("kernelcomm_result",
  slots = list(
    # Raw computation (always present)
    h = "numeric",
    observed = "numeric",
    permutations = "matrix",
    coords = "matrix",
    gene1_name = "character",
    gene2_name = "character",
    sample_id = "character",
    seed = "numeric",
    n_cells = "integer",
    call = "language",
    
    # Significance (added by compute_significance())
    p_values = "optional numeric",
    z_scores = "optional numeric",
    perm_mean = "optional numeric",
    perm_sd = "optional numeric",
    ci_lower = "optional numeric",
    ci_upper = "optional numeric",
    significant = "optional logical",
    p_adjusted = "optional numeric"
  )
)
```

### Batch Result
```r
setClass("kernelcomm_batch",
  slots = list(
    results = "list",  # List of kernelcomm_result objects
    metadata = "optional data.frame",  # Sample/gene metadata
    call = "language"
  )
)
```

---

## What NOT to Add (Scope Control) ❌

1. **Too many kernel types** - Gaussian + Uniform cover 95% of use cases
2. **GUI/Shiny apps** - Stay computational
3. **Database connectivity** - Export to files, let others handle storage
4. **Web services** - Not core competency
5. **Machine learning integration** - Stay focused on spatial autocorrelation
6. **Interactive visualization** - ggplot2 is sufficient
7. **Report generation** - Users can make their own reports
8. **Power analysis** - Too speculative for now
9. **Bayesian methods** - Different paradigm, out of scope
10. **Spatiotemporal extensions** - Future work, not now

---

## Design Decisions - Resolved

### ✅ 1. Single S4 Object with Smart Routing

**Decision**: Single `kernelcomm_result` class. No separate batch class.

**Rationale**:
- Simpler API: users don't need to learn two classes
- `multiScaleAC()` can detect if input is batch and route accordingly
- Batch processing returns `list` of `kernelcomm_result` objects
- Easier to maintain, document, and test

**Implementation**:
```r
# Single sample
result <- multiScaleAC(coords, gene1, gene2, h = c(10, 20, 40))
# Returns: kernelcomm_result object

# Batch (list of samples)
results <- multiScaleAC_batch(data_list, h_vec = c(10, 20, 40))
# Returns: list of kernelcomm_result objects
```

---

### ✅ 2. Optional Slots (NULL when not computed)

**Decision**: Significance slots are `NULL` until `compute_significance()` is called.

**Rationale**:
- `NULL` clearly means "not computed yet"
- `NA` suggests "computed but missing/error"
- Validation enforces: all significance slots NULL or all non-NULL

**Implementation**:
```r
setClass("kernelcomm_result",
  slots = list(
    # Required (always present)
    h = "numeric",
    observed = "numeric",
    permutations = "matrix",
    # ...
    
    # Optional (NULL until computed)
    p_values = "optional",
    z_scores = "optional",
    perm_mean = "optional",
    # ...
  )
)

# After multiScaleAC():
result@p_values  # NULL

# After compute_significance():
result@p_values  # numeric vector
```

---

### ✅ 3. Parallel Processing: `parallel::mclapply()` with Memory Management

**Decision**: Use `parallel::mclapply()` (built-in) with careful memory management.

**Rationale**:
- No new dependencies
- Good performance on Mac/Linux
- Windows falls back to sequential with warning
- Memory concerns addressed by loading data inside workers

**Memory management**:
1. Load data inside worker function (not in parent environment)
2. Use `mc.preschedule = TRUE` for uneven workloads
3. Call `gc()` after each sample
4. For very large batches: process in chunks, save intermediates

**Implementation**: See section 2.6 for full code example.

**Estimated memory**: ~200-500MB per core for typical dataset (7000 cells, 999 perms)

---

### ✅ 4. Flexible Heatmap Plotting

**Decision**: Accept multiple input types, user chooses row/column variables.

```r
plot_significance_heatmap <- function(data, 
                                       row_var = "sample_id",
                                       col_var = "h",
                                       value_var = "z_score",
                                       ...)
```

**Accepts**:
- `kernelcomm_result` (single sample)
- `list` of `kernelcomm_result` (batch)
- `data.frame` (from `export_mxfda_batch`)

**User chooses**:
- `row_var`: samples, gene pairs, or interaction of both
- `col_var`: bandwidths
- `value_var`: z_score, p_value, observed, etc.

---

### ✅ 5. Spatial Neighbor Plot

**Decision**: Simple point plot with optional connection lines.

```r
plot_spatial_neighbors <- function(coords, h, 
                                    gene1 = NULL,
                                    max_connections = 10,
                                    connection_alpha = 0.1)
```

**For 10,000 cells**:
- Points: fast with ggplot2 (no subsampling needed)
- Lines: limit to 10 connections per cell for readability
- Faint lines (alpha = 0.1) to avoid clutter

---

## Implementation Priority

### Week 1: Foundation ✅ COMPLETE
- [x] S4 class definitions
- [x] Refactor multiScaleAC() to return raw results
- [x] compute_significance() function
- [x] Basic print/summary methods

### Week 2: Export & Plotting ✅ COMPLETE
- [x] export_mxfda()
- [x] plot_scale_curve() (3 types)
- [x] plot.kernelcomm_result method

### Week 3: Batch Processing 🔄 PARTIAL
- [ ] kernelcomm_batch S4 class *(deferred: using list of kernelcomm_result objects per design decision in section ✅1)*
- [ ] multiScaleAC_batch() *(deferred to Phase 7)*
- [x] export_mxfda_batch()
- [ ] plot_batch_comparison() *(deferred to Phase 7)*

### Week 4: Polish ✅ COMPLETE
- [x] suggest_bandwidths()
- [x] validate_inputs()
- [x] adjust_pvalues()
- [x] Unit tests
- [x] Documentation

**Total:** Weeks 1, 2, 4 complete. Week 3 partial (export_mxfda_batch done, batch processing functions deferred).

---

## Implementation Status Summary

**Completed (Phase 1-6):**
- ✅ S4 class infrastructure with validity checks
- ✅ Modular design: multiScaleAC() → compute_significance() → adjust_pvalues()
- ✅ S3 methods: print(), summary(), plot()
- ✅ mxfda export: export_mxfda(), export_mxfda_batch()
- ✅ Helper functions: validate_inputs(), suggest_bandwidths(), qc_neighborhoods()
- ✅ All tests passing with ccRCCExample data

**Deferred (Phase 7+):**
- ⏳ multiScaleAC_batch() - batch processing orchestration
- ⏳ plot_batch_comparison() - batch visualization
- ⏳ kernelcomm_batch S4 class - intentionally using list instead (simpler API)

---

## Summary of Changes from V1

| Feature | V1 Plan | V2 Plan | Rationale |
|---------|---------|---------|-----------|
| Result structure | List | S4 object | Better validation, methods |
| Significance | In multiScaleAC() | Separate function | Modularity, flexibility |
| Batch output | List | S4 batch object | Structure for downstream |
| Bandwidth helpers | Phase 4 | Phase 1 | Essential for usability |
| Multiple testing | Not mentioned | Phase 1 | Critical for 1100 gene pairs |
| Caching | Not mentioned | Phase 3 | Needed for large batches |
| Diagnostics | Phase 3 | Phase 2 | Build trust in results |

---

**Ready for your feedback!** What should we adjust before starting implementation?
