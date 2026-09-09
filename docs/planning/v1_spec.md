# multiScaleAC Specification v1.0

**Package:** KernelCommunication  
**Date:** 2026-06-25  
**Revised:** 2026-06-26  
**Status:** Approved for Implementation  

---

## Executive Summary

This document specifies `multiScaleAC()`, a new optimized function for computing bivariate Moran's I across multiple kernel bandwidth scales. The implementation uses RcppEigen for core computations, supports both permutation-based and analytical inference methods, and handles datasets from hundreds to hundreds of thousands of cells via sparse matrix operations.

### Key Improvements Over `soupir_moran_multi()`

| Feature | Old Function | New Function |
|---------|--------------|--------------|
| Backend | Pure R | RcppEigen (C++) |
| Matrix storage | Dense (O(N²)) | Sparse (O(N × neighbors)) |
| Max practical N | ~5,000 cells | ~100,000+ cells |
| Inference methods | Permutation only | Permutation + Analytical |
| Analytical assumptions | N/A | Randomization + Normality |
| Progress reporting | None | Console progress bar |
| Memory optimization | None | Batched permutations |
| Large dataset mode | N/A | Auto-activate at N > 10,000 |
| Kernel functions | Gaussian only | Built-in: Gaussian, Uniform; Custom: user-provided |
| Distance metrics | Euclidean only | Built-in: Euclidean, Manhattan; Custom: user-provided |
| Deprecation status | Deprecated | Current |

---

## Architecture Overview

### System Boundaries

```
┌─────────────────────────────────────────────────────────────┐
│                      R Layer (multi_scale_ac.R)             │
│  ┌─────────────────────────────────────────────────────┐    │
│  │ • Input validation (user-friendly errors)           │    │
│  │ • Parameter sanitization                            │    │
│  │ • Permutation matrix generation (R side)            │    │
│  │ • Progress bar orchestration                        │    │
│  │ • Output formatting                                 │    │
│  └─────────────────────────────────────────────────────┘    │
│                            ↓                                  │
│              moran_multi_cpp() call                           │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│                   C++ Layer (moran_multi_cpp.cpp)           │
│  ┌─────────────────────────────────────────────────────┐    │
│  │ • Coordinate centering (numerical stability)        │    │
│  │ • Squared distance matrix (computed once)           │    │
│  │ • Sparse kernel construction per bandwidth          │    │
│  │ • Row normalization during construction             │    │
│  │ • Sparse-dense matrix multiplication                │    │
│  │ • Analytical variance computation                   │    │
│  │ • Progress reporting (Rcpp::Rcout)                  │    │
│  └─────────────────────────────────────────────────────┘    │
└─────────────────────────────────────────────────────────────┘
```

### Data Flow

```
Input Validation (R)
       ↓
Parameter Sanitization (R)
       ↓
Permutation Generation (R) → Y_perm matrix (N × perms)
       ↓
C++ Function Call
       ↓
┌──────────────────────────────────────┐
│ For each bandwidth h in h_vec:       │
│   1. Build sparse kernel W (N × N)   │
│   2. Row-normalize W                 │
│   3. U = W × Y_perm (batched)        │
│   4. I = (g1' × U) / denom           │
│   5. If analytical: compute Var(I)   │
└──────────────────────────────────────┘
       ↓
Output Formatting (R)
       ↓
Return List
```

---

## API Reference

### R Function: `multiScaleAC()`

```r
multiScaleAC(
  coords,
  gene1,
  gene2,
  h,
  perms = 1000,
  seed = 333,
  eps = 0,
  analytical = FALSE,
  assumption = c("randomization", "normality"),
  batch_size = 100,
  verbose = NULL,
  large_mode_n = 10000,
  kernel = c("gaussian", "uniform"),
  dist_method = c("euclidean", "manhattan"),
  kernel_func = NULL,
  dist_func = NULL
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `coords` | numeric matrix | (required) | N × 2 matrix of spatial coordinates (x, y) |
| `gene1` | numeric vector | (required) | Expression values for gene 1 (length N) |
| `gene2` | numeric vector | (required) | Expression values for gene 2 (length N) |
| `h` | numeric vector | (required) | Kernel bandwidths (must be > 0) |
| `perms` | integer | 1000 | Number of permutations (use 0 with `analytical=TRUE`) |
| `seed` | integer | 333 | Random seed for reproducibility |
| `eps` | numeric | 0 | Epsilon threshold for weight pruning (0 = no pruning) |
| `analytical` | logical | FALSE | Use analytical mean/variance instead of permutations |
| `assumption` | character | "randomization" | Distributional assumption: "randomization" or "normality" |
| `batch_size` | integer | 100 | Permutations per batch (memory optimization) |
| `verbose` | logical | NULL | Show progress (NULL = auto for large datasets) |
| `large_mode_n` | integer | 10000 | N threshold for automatic sparse mode |
| `kernel` | character | "gaussian" | Built-in kernel: "gaussian" or "uniform" (ignored if `kernel_func` provided) |
| `dist_method` | character | "euclidean" | Built-in distance: "euclidean" or "manhattan" (ignored if `dist_func` provided) |
| `kernel_func` | function | NULL | Custom kernel function (see Custom Functions section); disables large mode optimization |
| `dist_func` | function | NULL | Custom distance function (see Custom Functions section); disables large mode optimization |

**Returns:** List with components:

| Component | Type | Always Present | Description |
|-----------|------|----------------|-------------|
| `h` | numeric vector | Yes | Bandwidth values used |
| `observed` | numeric vector | Yes | Moran's I for each bandwidth |
| `permutations` | numeric matrix | Yes | n_h × perms matrix (empty if analytical only) |
| `expectation` | numeric vector | If analytical | E[I] per bandwidth |
| `variance` | numeric vector | If analytical | Var(I) per bandwidth |
| `z_score` | numeric vector | If analytical | Z-score per bandwidth |
| `p_value` | numeric vector | If analytical | One-sided p-value per bandwidth |
| `method` | character | Yes | "permutation" or "analytical" |
| `assumption` | character | If analytical | "randomization" or "normality" |

---

### C++ Function: `moran_multi_cpp()`

```cpp
// [[Rcpp::export]]
Rcpp::List moran_multi_cpp(
    const Eigen::MatrixXd& coords,      // N × 2
    const Eigen::VectorXd& gene1,       // N
    const Eigen::VectorXd& gene2,       // N
    const Eigen::VectorXd& h_vec,       // n_h bandwidths
    int perms,                          // total permutations
    int seed,                           // random seed
    double eps,                         // epsilon threshold
    bool analytical,                    // use analytical mean/variance
    std::string assumption,             // "randomization" or "normality"
    int batch_size,                     // permutations per batch
    bool verbose,                       // progress reporting
    int large_mode_n,                   // N threshold for large mode
    std::string kernel,                 // "gaussian" or "uniform"
    std::string dist_method             // "euclidean" or "manhattan"
);
```

**Note:** Custom kernel/distance functions (`kernel_func`, `dist_func`) are handled in the R layer and bypass the C++ function. The R wrapper computes weights/distances in R and passes pre-computed sparse matrices to a separate C++ entry point.

---

## Parameter Specifications

### Input Validation Rules

| Parameter | Validation | Error Message |
|-----------|------------|---------------|
| `coords` | Must be numeric matrix with 2 columns | "coords must be a numeric matrix with 2 columns" |
| `coords` | No NA/Inf values | "Coordinates contain missing or infinite values" |
| `gene1`, `gene2` | Must be numeric vectors | "gene1 and gene2 must be numeric vectors" |
| `gene1`, `gene2` | No NA/Inf values | "Gene expression values contain missing or infinite values" |
| All inputs | Same length N | "coords, gene1, and gene2 must have same number of observations" |
| `h` | Must have at least one value | "h must contain at least one bandwidth value" |
| `h` | Filter negative values | Warning: "Removed X negative bandwidth values" |
| `h` | Filter zero values | Warning: "Removed X zero bandwidth values (undefined kernel)" |
| `h` | At least one valid after filtering | "No valid bandwidth values remaining" |
| `perms` | Integer ≥ 0 | "perms must be a non-negative integer" |
| `perms` | ≥ 1 if !analytical | "perms must be >= 1 when method='permutation'" |
| `seed` | Integer | "seed must be an integer" |
| `eps` | Numeric ≥ 0 | "eps must be a non-negative numeric value" |
| `batch_size` | Integer ≥ 1 | "batch_size must be a positive integer" |
| `large_mode_n` | Integer ≥ 3 | "large_mode_n must be an integer >= 3" |
| `kernel` | Must be "gaussian" or "uniform" (if `kernel_func` is NULL) | "kernel must be 'gaussian' or 'uniform'" |
| `dist_method` | Must be "euclidean" or "manhattan" (if `dist_func` is NULL) | "dist_method must be 'euclidean' or 'manhattan'" |
| `kernel_func` | Must be a function (if provided) | "kernel_func must be a function" |
| `dist_func` | Must be a function (if provided) | "dist_func must be a function" |
| N | ≥ 3 | "At least 3 observations required" |
| `kernel_func` + `large_mode_n` | Warning if N > large_mode_n | "Custom kernel disables large dataset optimization; computation may be slow" |
| `dist_func` + `large_mode_n` | Warning if N > large_mode_n | "Custom distance method disables large dataset optimization; computation may be slow" |

### Derived Parameters

| Parameter | Computation | Purpose |
|-----------|-------------|---------|
| `N` | `nrow(coords)` | Number of observations |
| `n_h` | `length(h)` | Number of bandwidths |
| `large_mode` | `N > large_mode_n && is.null(kernel_func) && is.null(dist_func)` | Enable sparse mode + distance pruning (disabled for custom functions) |
| `verbose` | `if (is.null(verbose)) verbose <- large_mode` | Auto-enable for large datasets |
| `denom` | `sum(g1_s^2)` | Denominator for Moran's I formula |
| `use_custom_kernel` | `!is.null(kernel_func)` | Use R-based custom kernel (vs built-in C++) |
| `use_custom_dist` | `!is.null(dist_func)` | Use R-based custom distance (vs built-in C++) |

---

## Algorithm Details

### 0. Custom Function Handling (R Layer)

**For custom kernel functions:**
```r
if (!is.null(kernel_func)) {
  # Validate function signature
  if (!is.function(kernel_func)) stop("kernel_func must be a function")
  
  # Warn about large mode disable
  if (N > large_mode_n) {
    warning("Custom kernel disables large dataset optimization; computation may be slow")
  }
  
  # For each bandwidth, call kernel_func in R
  # kernel_func must accept: kernel_func(dist_matrix, h)
  # Returns: N × N weight matrix (dense or sparse Matrix::dsCMatrix)
  for (h_idx in seq_along(h_vec)) {
    dist_matrix <- compute_distance_matrix(coords, dist_method, dist_func)
    W <- kernel_func(dist_matrix, h_vec[h_idx])
    
    # If dense, convert to sparse if needed
    if (!inherits(W, "dsCMatrix")) {
      W <- Matrix::Matrix(W, sparse = TRUE)
    }
    
    # Pass W to C++ for permutation multiplication
    # ...
  }
}
```

**For custom distance functions:**
```r
if (!is.null(dist_func)) {
  # Validate function signature
  if (!is.function(dist_func)) stop("dist_func must be a function")
  
  # Warn about large mode disable
  if (N > large_mode_n) {
    warning("Custom distance method disables large dataset optimization; computation may be slow")
  }
  
  # Compute distance matrix once in R
  # dist_func must accept: dist_func(coords)
  # Returns: N × N distance matrix (or squared distance matrix)
  dist_matrix <- dist_func(coords)
}
```

### 1. Squared Distance Matrix (Built-in Methods, Computed Once)

**Euclidean Distance (squared):**
```cpp
// Center coordinates for numerical stability
Eigen::VectorXd coords_mean = coords.colwise().mean();
coords.rowwise() -= coords_mean.transpose();

// Compute squared Euclidean distances
Eigen::MatrixXd dist_sq(N, N);
for (int i = 0; i < N; i++) {
  for (int j = i; j < N; j++) {
    double dx = coords(i, 0) - coords(j, 0);
    double dy = coords(i, 1) - coords(j, 1);
    double d2 = dx * dx + dy * dy;
    dist_sq(i, j) = d2;
    dist_sq(j, i) = d2;  // Symmetric
  }
}
```

**Manhattan Distance:**
```cpp
// Compute Manhattan distances (not squared - used directly)
Eigen::MatrixXd dist(N, N);
for (int i = 0; i < N; i++) {
  for (int j = i; j < N; j++) {
    double dx = std::abs(coords(i, 0) - coords(j, 0));
    double dy = std::abs(coords(i, 1) - coords(j, 1));
    double d = dx + dy;
    dist(i, j) = d;
    dist(j, i) = d;  // Symmetric
  }
}
```

**Note:** For Manhattan distance with Gaussian kernel, use `exp(-d²/(2h²))` where `d` is Manhattan distance (not squared Euclidean).

### 2. Sparse Kernel Construction (Per Bandwidth, Built-in Kernels)

**Gaussian Kernel:**
```cpp
for each bandwidth h in h_vec:
  std::vector<Triplet<double>> triplets;
  Eigen::VectorXd rowsums(N);
  rowsums.setZero();
  
  const double inv2h2 = 1.0 / (2.0 * h * h);
  
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i == j) continue;  // Skip diagonal
      
      double w = std::exp(-dist_sq(i, j) * inv2h2);
      
      if (w > eps) {  // Epsilon threshold (0 = include all)
        rowsums(i) += w;
        triplets.push_back(Triplet<double>(i, j, w));
      }
    }
  }
  
  // Build sparse matrix
  Eigen::SparseMatrix<double> W(N, N);
  W.setFromTriplets(triplets.begin(), triplets.end());
  
  // Row-normalize: W(i,j) = W(i,j) / rowsum(i)
  for (int i = 0; i < N; i++) {
    if (rowsums(i) > 0) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(W, i); it; ++it) {
        it.valueRef() /= rowsums(i);
      }
    }
  }
```

**Uniform Kernel:**
```cpp
for each bandwidth h in h_vec:
  std::vector<Triplet<double>> triplets;
  Eigen::VectorXd rowsums(N);
  rowsums.setZero();
  
  const double h_sq = h * h;
  
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i == j) continue;  // Skip diagonal
      
      // Uniform: w = 1 if d <= h, else 0
      if (dist_sq(i, j) <= h_sq) {
        rowsums(i) += 1.0;
        triplets.push_back(Triplet<double>(i, j, 1.0));
      }
    }
  }
  
  // Build and row-normalize (same as Gaussian)
  // ...
```

**Kernel Formulas:**

| Kernel | Formula | Support | Notes |
|--------|---------|---------|-------|
| Gaussian | `exp(-d²/(2h²))` | Infinite | Smooth decay, standard choice |
| Uniform | `1` if `d ≤ h`, else `0` | Compact | Binary neighbors, sharp cutoff |

**Note:** For uniform kernel, all non-zero weights are equal, so row-normalization gives `W(i,j) = 1/n_neighbors(i)`.

### 3. Batched Permutation Handling

```cpp
int n_batches = (perms + batch_size - 1) / batch_size;
Eigen::MatrixXd all_perm_vals(n_h, perms);

for (int b = 0; b < n_batches; b++) {
  int batch_start = b * batch_size;
  int batch_end = std::min(batch_start + batch_size, perms);
  int current_batch_size = batch_end - batch_start;
  
  // Y_batch: N × current_batch_size (passed from R)
  Eigen::MatrixXd U_batch = W * Y_batch;  // Sparse-dense multiply
  
  // Compute dot products: (g1' × U_batch) / denom
  Eigen::RowVectorXd perm_vals = (g1_s.transpose() * U_batch) / denom;
  
  // Store results
  all_perm_vals.block(h_idx, batch_start, 1, current_batch_size) = perm_vals;
}
```

**Memory Savings:**

| Scenario | Unbatched | Batched (100) | Savings |
|----------|-----------|---------------|---------|
| N=20,000, perms=1,000 | 160 MB | 16 MB | 10× |
| N=100,000, perms=10,000 | 8 GB | 80 MB | 100× |

### 4. Analytical Mean and Variance

**Expectation (both assumptions):**

```
E[I] = -1 / (N - 1)
```

**Variance under Randomization:**

```
Var(I) = [N·S₄ - S₃·S₅] / [(N-1)(N-2)(N-3)·S₀²] - E[I]²

where:
  S₀ = Σᵢ Σⱼ wᵢⱼ              (sum of all weights)
  S₁ = ½ Σᵢ Σⱼ (wᵢⱼ + wⱼᵢ)²
  S₂ = Σᵢ (Σⱼ wᵢⱼ + Σⱼ wⱼᵢ)²
  S₃ = [N⁻¹ Σᵢ(xᵢ-x̄)⁴] / [N⁻¹ Σᵢ(xᵢ-x̄)²]²   (kurtosis of gene1)
  S₄ = (N² - 3N + 3)·S₁ - N·S₂ + 3·S₀²
  S₅ = (N² - N)·S₁ - 2N·S₂ + 6·S₀²
```

**Variance under Normality:**

```
Var(I) = [N²·S₁ - N·S₂ + 3·S₀²] / [(N-1)(N+1)·S₀²] - E[I]²
```

**For Row-Normalized Weights:**
- `S₀ = N` (each row sums to 1)
- Simplifies variance formulas

**Z-score and P-value:**

```
z = (I - E[I]) / sqrt(Var(I))
p = 1 - Φ(z)  // one-sided, upper tail (pnorm in R)
```

### 5. Progress Reporting

```cpp
if (verbose && n_h > 1) {
  int progress_interval = std::max(1, n_h / 20);  // Report every 5%
  
  for (int h_idx = 0; h_idx < n_h; h_idx++) {
    // ... compute bandwidth h_idx ...
    
    if ((h_idx + 1) % progress_interval == 0 || h_idx == n_h - 1) {
      int pct = (h_idx + 1) * 100 / n_h;
      Rcpp::Rcout << "\r  Computing bandwidths: " << pct 
                  << "% (" << (h_idx + 1) << "/" << n_h << ")" 
                  << std::flush;
    }
  }
  Rcpp::Rcout << "\n" << std::flush;
}
```

---

## Custom Functions

### Custom Kernel Functions

**Interface:**
```r
kernel_func <- function(dist_matrix, h) {
  # dist_matrix: N × N matrix of pairwise distances
  # h: scalar bandwidth for current iteration
  # Returns: N × N weight matrix (dense matrix or sparse Matrix::dsCMatrix)
}
```

**Example - Exponential Kernel:**
```r
exponential_kernel <- function(dist_matrix, h) {
  exp(-dist_matrix / h)  # Note: uses distance, not squared distance
}

result <- multiScaleAC(
  coords = coords,
  gene1 = g1,
  gene2 = g2,
  h = c(5, 10, 20),
  kernel_func = exponential_kernel
)
```

**Example - Bisquare Kernel:**
```r
bisquare_kernel <- function(dist_matrix, h) {
  u <- dist_matrix / h
  w <- ifelse(u <= 1, (1 - u^2)^2, 0)
  return(w)
}
```

**Performance Note:** Custom kernel functions are called once per bandwidth in R. For large datasets (N > 10,000), this bypasses the C++ sparse optimization and may be significantly slower. Consider using built-in kernels for large datasets.

---

### Custom Distance Functions

**Interface:**
```r
dist_func <- function(coords) {
  # coords: N × 2 matrix of coordinates
  # Returns: N × N distance matrix (or squared distance matrix)
}
```

**Example - Minkowski Distance (p=3):**
```r
minkowski_dist <- function(coords, p = 3) {
  N <- nrow(coords)
  dist_matrix <- matrix(0, N, N)
  for (i in 1:N) {
    for (j in i:N) {
      d <- sum(abs(coords[i, ] - coords[j, ])^p)^(1/p)
      dist_matrix[i, j] <- d
      dist_matrix[j, i] <- d
    }
  }
  return(dist_matrix)
}

result <- multiScaleAC(
  coords = coords,
  gene1 = g1,
  gene2 = g2,
  h = c(5, 10, 20),
  dist_func = function(coords) minkowski_dist(coords, p = 3)
)
```

**Example - Fixed Radius Distance:**
```r
fixed_radius_dist <- function(coords) {
  # All pairs within radius R have distance 1, others 0
  # Useful for defining neighbors by threshold
  R <- 10
  dist_matrix <- as.matrix(dist(coords))
  dist_matrix[dist_matrix <= R] <- 1
  dist_matrix[dist_matrix > R] <- 0
  return(dist_matrix)
}
```

**Performance Note:** Custom distance functions are called once (before bandwidth loop). The returned matrix is used for all bandwidths. Memory usage is O(N²) for dense distance matrices.

---

### Combining Custom Kernel and Distance

```r
result <- multiScaleAC(
  coords = coords,
  gene1 = g1,
  gene2 = g2,
  h = c(5, 10, 20),
  kernel_func = exponential_kernel,
  dist_func = minkowski_dist
)
```

**Warning:** Using both custom functions disables all C++ optimizations. Recommended only for small datasets (N < 5,000) or when custom methods are essential.

---

## Error Handling Matrix

| Scenario | Severity | Message | Recovery |
|----------|----------|---------|----------|
| NA in coords | Fatal | "Coordinates contain missing or infinite values" | None |
| Inf in coords | Fatal | "Coordinates contain missing or infinite values" | None |
| NA in gene1/gene2 | Fatal | "Gene expression values contain missing or infinite values" | None |
| Length mismatch | Fatal | "coords, gene1, and gene2 must have same number of observations" | None |
| N < 3 | Fatal | "At least 3 observations required" | None |
| No valid h values | Fatal | "No valid bandwidth values remaining" | None |
| perms < 1 (when !analytical) | Fatal | "perms must be >= 1 when method='permutation'" | None |
| Negative h values | Warning | "Removed X negative bandwidth values" | Continue with filtered h |
| Zero h values | Warning | "Removed X zero bandwidth values (undefined kernel)" | Continue with filtered h |
| Isolated points (zero row sum) | Warning (once) | "X cells have no neighbors within bandwidth; lag values set to 0" | Set lag = 0 for isolated cells |
| All weights pruned (eps too high) | Warning | "All weights pruned for bandwidth h=X; consider increasing eps or h" | Return NA for that bandwidth |
| Invalid `kernel` string | Fatal | "kernel must be 'gaussian' or 'uniform'" | None |
| Invalid `dist_method` string | Fatal | "dist_method must be 'euclidean' or 'manhattan'" | None |
| `kernel_func` not a function | Fatal | "kernel_func must be a function" | None |
| `dist_func` not a function | Fatal | "dist_func must be a function" | None |
| Custom function with large N | Warning | "Custom kernel/distance disables large dataset optimization; computation may be slow" | Continue with dense computation |
| Custom kernel returns invalid matrix | Fatal | "kernel_func must return a square matrix with dimensions N × N" | None |
| Custom distance returns invalid matrix | Fatal | "dist_func must return a square matrix with dimensions N × N" | None |

**Isolated Point Handling:**

```cpp
// Check for zero row sums after kernel construction
for (int i = 0; i < N; i++) {
  if (rowsums(i) == 0) {
    isolated_count++;
    // Lag value will be 0 (no neighbors contribute)
  }
}

if (isolated_count > 0 && verbose) {
  Rcpp::Rcout << "Warning: " << isolated_count 
              << " cells have no neighbors within bandwidth\n";
}
```

---

## Module Breakdown

### File Structure

```
KernelCommunication/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── multi_scale_ac.R        # NEW: multiScaleAC function
│   ├── soupir_moran_multi.R    # MODIFY: Add deprecation warning
│   ├── RcppExports.R           # AUTO: Add moran_multi_cpp binding
│   └── ...
├── src/
│   ├── moran_multi_cpp.cpp     # NEW: C++ implementation
│   ├── RcppExports.cpp         # AUTO: Add moran_multi_cpp binding
│   └── ...
├── tests/
│   └── testthat/
│       ├── test-multi_scale_ac.R    # NEW: Comprehensive tests
│       └── test-deprecated.R        # NEW: Backward compatibility tests
├── benchmarks/
│   └── benchmark_multi_scale.R      # NEW: Performance benchmarks
└── docs/
    └── planning/
        └── v1_spec.md          # THIS DOCUMENT
```

### File-by-File Implementation

#### 1. `DESCRIPTION` (Modify)

```
Package: KernelCommunication
Version: 1.1.0
Imports:
    Rcpp (>= 1.0.13),
    Matrix,
    progress,
    lifecycle
LinkingTo: Rcpp, RcppEigen
Suggests:
    testthat (>= 3.0.0)
```

#### 2. `NAMESPACE` (Modify)

```
export(multiScaleAC)
export(soupir_moran_multi)
importFrom(lifecycle, deprecate_warn)
useDynLib(KernelCommunication, .registration = TRUE)
```

#### 3. `R/multi_scale_ac.R` (Create)

**Responsibilities:**
- Input validation with user-friendly error messages
- Parameter sanitization (filter negative/zero h)
- Kernel/distance method validation
- Custom function validation and handling
- Permutation matrix generation (single seed)
- C++ function call (built-in kernels/distances)
- R-based computation (custom kernels/distances)
- Progress bar integration (optional)
- Output formatting
- Warning aggregation
- Large mode auto-disable for custom functions

**Key Functions:**
```r
multiScaleAC <- function(...) { ... }
.validate_inputs <- function(...) { ... }
.sanitize_h <- function(h) { ... }
.validate_kernel <- function(kernel, kernel_func, N, large_mode_n) { ... }
.validate_distance <- function(dist_method, dist_func, N, large_mode_n) { ... }
.compute_distance_matrix <- function(coords, dist_method, dist_func) { ... }
.compute_kernel_weights <- function(dist_matrix, h, kernel, kernel_func) { ... }
.generate_permutations <- function(gene2, perms, seed) { ... }
.format_output <- function(result, analytical, assumption) { ... }
```

#### 4. `src/moran_multi_cpp.cpp` (Create)

**Responsibilities:**
- Coordinate centering
- Distance matrix computation (built-in methods: euclidean, manhattan)
- Sparse kernel construction (built-in kernels: gaussian, uniform)
- Row normalization
- Batched permutation multiplication
- Analytical variance computation
- Progress reporting

**Key Functions:**
```cpp
Rcpp::List moran_multi_cpp(...) { ... }
Eigen::MatrixXd compute_distance_matrix_euclidean(...) { ... }
Eigen::MatrixXd compute_distance_matrix_manhattan(...) { ... }
Eigen::SparseMatrix<double> build_gaussian_kernel(...) { ... }
Eigen::SparseMatrix<double> build_uniform_kernel(...) { ... }
void compute_analytical_variance(...) { ... }
```

**Note:** Custom kernel/distance functions bypass this C++ module and are computed in R.

#### 5. `R/soupir_moran_multi.R` (Modify)

**Changes:**
- Add `@keywords internal` tag
- Add deprecation warning via `lifecycle::deprecate_warn()`
- Wrap `multiScaleAC()` call with compatible defaults
- Preserve original output format

#### 6. `tests/testthat/test-multi_scale_ac.R` (Create)

**Test Categories:**
1. Basic functionality
2. Backward compatibility
3. Analytical method
4. Analytical vs permutation convergence
5. Randomization vs normality
6. Error handling
7. Large mode
8. Batched permutations
9. Isolated points
10. Deprecation warning
11. Built-in kernels (gaussian, uniform)
12. Built-in distance methods (euclidean, manhattan)
13. Custom kernel functions
14. Custom distance functions
15. Custom function large mode warnings

#### 7. `tests/testthat/test-deprecated.R` (Create)

**Test Categories:**
1. Old function still produces correct output
2. Deprecation warning is shown

#### 8. `benchmarks/benchmark_multi_scale.R` (Create)

**Benchmark Configurations:**
- N = [5000, 10000, 20000, 50000, 100000]
- n_h = [26, 100, 500]
- perms = [100, 1000, 10000]
- analytical = [TRUE, FALSE]
- kernel = ["gaussian", "uniform"]
- dist_method = ["euclidean", "manhattan"]
- Custom function overhead (N = [1000, 5000] only due to R callback cost)

---

## Testing Strategy

### Test Framework: `testthat`

**Test File Organization:**
```
tests/
└── testthat/
    ├── test-multi_scale_ac.R    # Main test suite
    ├── test-deprecated.R        # Backward compatibility
    └── setup.R                  # Test fixtures (optional)
```

### Test Coverage Goals

| Component | Target Coverage |
|-----------|-----------------|
| Input validation | 100% of error paths |
| Parameter sanitization | 100% of edge cases |
| Core algorithm | 90%+ (via output comparison) |
| Error handling | 100% of documented scenarios |
| Output formatting | 100% of return types |

### Test Data

**Primary Dataset:** `ccRCCExample` (7,122 cells)

**Synthetic Datasets:**
```r
# Small (quick tests)
N = 100

# Medium (integration tests)
N = 500

# Large (performance tests)
N = 20000

# Edge cases
N = 3 (minimum)
N = 10 with isolated points
```

### Test Execution

```bash
# Run all tests
R -e "devtools::test()"

# Run specific test file
R -e "testthat::test_file('tests/testthat/test-multi_scale_ac.R')"

# Run with coverage
R -e "covr::package_coverage()"
```

---

## Performance Targets

### Benchmark Configurations

| Test ID | N | n_h | perms | Mode | Target Time |
|---------|---|-----|-------|------|-------------|
| B1 | 5,000 | 26 | 1,000 | Permutation | < 1 sec |
| B2 | 10,000 | 100 | 1,000 | Permutation | < 5 sec |
| B3 | 20,000 | 500 | 1,000 | Permutation | < 30 sec |
| B4 | 50,000 | 100 | 1,000 | Permutation | < 60 sec |
| B5 | 100,000 | 50 | 1,000 | Permutation | < 120 sec |
| B6 | 20,000 | 500 | 0 | Analytical | < 20 sec |
| B7 | 50,000 | 100 | 0 | Analytical | < 30 sec |

### Memory Targets

| N | Max Memory (RSS) |
|---|------------------|
| 5,000 | < 500 MB |
| 10,000 | < 1 GB |
| 20,000 | < 2 GB |
| 50,000 | < 4 GB |
| 100,000 | < 8 GB |

### Comparison: Old vs New

| Metric | `soupir_moran_multi()` | `multiScaleAC()` | Improvement |
|--------|------------------------|------------------|-------------|
| N=5,000, h=26, perms=100 | ~2 sec | < 0.5 sec | 4× |
| N=10,000, h=100, perms=100 | OOM | < 3 sec | N/A |
| N=20,000 | Impossible | < 30 sec | N/A |

---

## Backward Compatibility

### Deprecation Workflow

**Phase 1 (v1.1.0):** Soft Deprecation
- `soupir_moran_multi()` shows warning via `lifecycle::deprecate_warn()`
- Function remains fully functional
- Documentation updated with deprecation notice

**Phase 2 (v1.3.0):** Superseded
- Warning message strengthened
- Vignette migration guide added

**Phase 3 (v2.0.0):** Defunct (Optional)
- Function removed or made to error with migration message
- Only if user feedback supports removal

### Deprecation Warning Example

```r
lifecycle::deprecate_warn(
  when = "1.1.0",
  what = "soupir_moran_multi()",
  with = "multiScaleAC()",
  details = "The new function supports analytical methods, large datasets, and progress reporting."
)
```

**User sees:**
```
Warning message:
`soupir_moran_multi()` was deprecated in KernelCommunication 1.1.0.
ℹ Please use `multiScaleAC()` instead.
ℹ The new function supports analytical methods, large datasets, and progress reporting.
```

### Wrapper Implementation

```r
soupir_moran_multi <- function(coords, gene1, gene2, h_vec, perms, seed = 333) {
  lifecycle::deprecate_warn(
    when = "1.1.0",
    what = "soupir_moran_multi()",
    with = "multiScaleAC()"
  )
  
  result <- multiScaleAC(
    coords = coords,
    gene1 = gene1,
    gene2 = gene2,
    h = h_vec,
    perms = perms,
    seed = seed,
    eps = 0,
    analytical = FALSE,
    verbose = FALSE
  )
  
  # Return format compatible with original
  list(
    h = result$h,
    observed = result$observed,
    permutations = result$permutations
  )
}
```

---

## Build and Test Commands

### Development Workflow

```bash
# 1. Create conda environment
conda create -n rcpp_dev -c conda-forge r-base r-devtools r-matrix r-rcpp r-rcppeigen
conda activate rcpp_dev

# 2. Install additional R packages
R -e "install.packages(c('progress', 'lifecycle', 'testthat', 'microbenchmark'), repos='https://cloud.r-project.org')"

# 3. Load package (no installation needed)
R -e "devtools::load_all()"

# 4. Regenerate Rcpp exports (if C++ changes)
R -e "Rcpp::compileAttributes()"

# 5. Run tests
R -e "devtools::test()"

# 6. Run benchmarks
R -e "source('benchmarks/benchmark_multi_scale.R')"

# 7. Check package
R -e "devtools::check()"

# 8. Build documentation
R -e "roxygen2::roxygenise()"
```

### CRAN Submission Checklist

- [ ] All tests pass (`R CMD check` clean)
- [ ] No NOTEs/WARNINGs (except expected deprecation)
- [ ] Vignettes build correctly
- [ ] Examples run without errors
- [ ] Memory usage within CRAN limits (< 2 GB typical)
- [ ] Cross-platform tested (macOS, Linux, Windows)
- [ ] NEWS file updated
- [ ] Man pages complete and accurate

---

## Appendix A: Mathematical Formulas

### Moran's I Formula

For bivariate Moran's I with row-normalized weights:

```
I = (g1_s' × W × g2_s) / (g1_s' × g1_s)

where:
  g1_s = scale(gene1)  # centered and scaled
  g2_s = scale(gene2)
  W = row-normalized spatial weight matrix
  W(i,j) = exp(-d(i,j)² / (2h²)) / Σ_k exp(-d(i,k)² / (2h²))
```

### Gaussian Kernel

```
w(i,j) = exp(-d(i,j)² / (2h²))

where:
  d(i,j)² = squared distance (method-dependent)
  h = bandwidth parameter
```

### Uniform Kernel

```
w(i,j) = 1  if d(i,j) ≤ h
       = 0  otherwise

where:
  d(i,j) = distance (method-dependent)
  h = bandwidth (radius threshold)
```

### Distance Methods

**Euclidean Distance:**
```
d(i,j) = √[(x_i - x_j)² + (y_i - y_j)²]

Squared Euclidean (used in Gaussian):
d(i,j)² = (x_i - x_j)² + (y_i - y_j)²
```

**Manhattan Distance:**
```
d(i,j) = |x_i - x_j| + |y_i - y_j|

Note: For Gaussian kernel with Manhattan distance:
w(i,j) = exp(-d_manhattan(i,j)² / (2h²))
```

### Weight Normalization

Row-normalization ensures each row sums to 1:

```
W_norm(i,j) = w(i,j) / Σ_k w(i,k)
```

This makes `S₀ = N` for row-normalized weights, simplifying variance formulas.

### Kernel Comparison

| Kernel | Formula | Support | Continuity | Use Case |
|--------|---------|---------|------------|----------|
| Gaussian | `exp(-d²/(2h²))` | Infinite | Smooth (C∞) | Standard, smooth spatial decay |
| Uniform | `1{d ≤ h}` | Compact | Discontinuous | Binary neighbor definition |

**Support:**
- Infinite: All points contribute (weights decay with distance)
- Compact: Only points within radius `h` contribute (natural sparsity)

---

## Appendix B: Glossary

| Term | Definition |
|------|------------|
| **Moran's I** | Measure of spatial autocorrelation; ranges from -1 (dispersed) to +1 (clustered) |
| **Bivariate** | Involves two variables (gene1 and gene2); measures spatial correlation between them |
| **Bandwidth (h)** | Kernel parameter controlling spatial scale; larger h = broader influence |
| **Row-normalized** | Each row of weight matrix sums to 1; ensures comparability across bandwidths |
| **Sparse matrix** | Matrix with mostly zero elements; stored efficiently as (row, col, value) triplets |
| **Epsilon threshold** | Minimum weight value; weights below this are set to 0 (pruned) |
| **Randomization** | Null hypothesis: gene values randomly assigned to fixed locations |
| **Normality** | Null hypothesis: gene values drawn from normal distribution |
| **Large mode** | Automatic sparse mode for N > 10,000; includes distance-based pruning |
| **Kernel** | Function that converts distances to weights; determines spatial influence pattern |
| **Gaussian kernel** | `exp(-d²/(2h²))`; smooth decay with distance, infinite support |
| **Uniform kernel** | `1{d ≤ h}`; binary weights, compact support (sharp cutoff at h) |
| **Euclidean distance** | Straight-line distance: `√(Δx² + Δy²)` |
| **Manhattan distance** | Grid-based distance: `|Δx| + |Δy|` |
| **Compact support** | Kernel is zero beyond finite distance (e.g., uniform); enables natural sparsity |
| **Infinite support** | Kernel is non-zero everywhere (e.g., gaussian); relies on epsilon pruning for sparsity |

---

## Appendix C: Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-06-25 | Alex Soupir | Initial specification |
| 1.1 | 2026-06-26 | Alex Soupir | Added kernel/distance flexibility: built-in gaussian/uniform kernels, euclidean/manhattan distances, custom function support via `kernel_func` and `dist_func` parameters |

---

## Appendix D: Open Questions

| Question | Status | Decision |
|----------|--------|----------|
| Progress bar: R or C++? | Resolved | C++ with `Rcpp::Rcout` |
| Permutation batching? | Resolved | Yes, default batch_size=100 |
| Analytical return format? | Resolved | Empty matrix (n_h × 0) for compatibility |
| Testing framework? | Resolved | `testthat` |
| Backward compatibility? | Resolved | Deprecation via `lifecycle` package |
| OpenMP parallelization? | Deferred | Single-threaded for v1.0 |
| Custom kernel interface? | Resolved | `kernel_func(dist_matrix, h)` → N×N matrix |
| Custom distance interface? | Resolved | `dist_func(coords)` → N×N matrix |
| Large mode with custom functions? | Resolved | Auto-disable with warning |
| Which built-in kernels? | Resolved | gaussian, uniform (extensible in future) |
| Which built-in distances? | Resolved | euclidean, manhattan (extensible in future) |

---

**END OF SPECIFICATION**
