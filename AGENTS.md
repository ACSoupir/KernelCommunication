# AGENTS.md

## Package Overview

R package with Rcpp/C++ backend for computing bivariate Moran's I (spatial autocorrelation) over kernel bandwidth scales. Uses RcppEigen for sparse matrix operations.

## Key Commands

```bash
# Build/install package
R CMD INSTALL .
# or from R session
devtools::install()

# Regenerate Rcpp exports (if C++ changes)
Rcpp::compileAttributes()

# Rebuild documentation
roxygen2::roxygenise()

# Load package in R session
devtools::load_all()
```

## Structure

- `R/` - R source files
  - `soupir_moran_bv.R` - single-bandwidth bivariate Moran's I
  - `soupir_moran_multi.R` - multi-bandwidth curve computation
  - `RcppExports.R` - auto-generated C++ bindings (do not edit)
  - `utils.R` - helper functions (sparse matrix conversion, lag operations)
  - `data.R` - documentation for `ccRCCExample` dataset
- `src/` - C++ source (RcppEigen)
  - `gaussian_kernel_weights.cpp` - sparse weight matrix from Gaussian kernel
  - `sparse_dense_mult.cpp` - sparse-dense matrix multiplication
  - `RcppExports.cpp` - auto-generated (do not edit)
- `data/` - packaged example dataset (`ccRCCExample.rda`)
- `data-raw/` - dataset preparation scripts

## Important Notes

- `RcppExports.R` and `RcppExports.cpp` are generated; run `Rcpp::compileAttributes()` after C++ changes
- Documentation is roxygen2-based; run `roxygen2::roxygenise()` to regenerate `.Rd` files and `NAMESPACE`
- Package requires R >= 3.5, Rcpp >= 1.0.13, RcppEigen (linking), Matrix

## Main Functions

- `soupir_moran_bv(h, coords, gene1, gene2, perms, eps)` - single bandwidth
- `soupir_moran_multi(coords, gene1, gene2, h_vec, perms, seed)` - multiple bandwidths
