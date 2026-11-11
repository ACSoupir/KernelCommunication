# AGENTS.md

## Package Overview

R package with an Rcpp/C++ backend computing bivariate Moran's I (spatial cross-correlation) as a
function of Gaussian kernel bandwidth, producing a functional curve `I(h)` rather than a single
value at one chosen scale. Implements the method in the companion preprint
*"multiScaleAC: Cell-Cell interaction with Moran's I as a function of kernel bandwidth"*
(bioRxiv 2026.07.14.738306).

**The R package and the GitHub repository are both named `multiScaleAC`
(https://github.com/ACSoupir/multiScaleAC).** The local working-copy directory on this machine is
still named `KernelCommunication` from before the repository rename; that is just a folder name on
disk, unrelated to git, and does not need to match.

## Non-negotiable rules

1. **`R/soupir_moran_multi.R` is the reference implementation and must never be modified.** It is the
   regression oracle: every fixture in `tests/testthat/fixtures/` is generated from it, and a test
   hashes it to detect edits. `soupir_moran_multi()` also stays exported permanently for backwards
   compatibility — it is not deprecated.
2. **`RcppExports.R`, `RcppExports.cpp`, `NAMESPACE` and `man/*.Rd` are generated. Never hand-edit
   them.** Run `Rcpp::compileAttributes()` after any C++ change, then `roxygen2::roxygenise()`.
3. **Work is local only. Never push to a remote.**

## Key Commands

```bash
R CMD INSTALL .                     # build and install
R CMD check --no-manual .           # run during development
R CMD check --as-cran .             # release gate
```

```r
Rcpp::compileAttributes()           # after any C++ change (reads DESCRIPTION -- rename first)
roxygen2::roxygenise()              # regenerate .Rd and NAMESPACE
testthat::test_local()              # run the test suite
```

Note: `devtools` is not installed in this environment (its `curl` dependency needs system libcurl
headers). Use the base `R CMD` tooling and `testthat::test_local()` directly rather than
`devtools::install()` / `devtools::test()`.

## Structure

- `R/`
  - `soupir_moran_multi.R` — **reference implementation, never modify**
  - `soupir_moran_bv.R` — single-bandwidth sparse variant
  - `multi_scale_ac.R` — `multiScaleAC()`, the main entry point
  - `classes.R` — S4 `kernelcomm_result` class and validity
  - `significance.R` — permutation p-values, multiple-testing correction
  - `moments.R`, `approx.R` — exact analytical permutation moments and tail approximations
  - `methods.R`, `plot.R` — `show`/`print`/`summary`/`plot` and the `plot_*` functions
  - `helpers.R` — input validation, bandwidth suggestion, neighbourhood QC
  - `export_mxfda.R` — bridge to `mxfda::add_summary_function()`
  - `utils.R` — sparse matrix conversion and lag operations
  - `data.R` — documentation for the `ccRCCExample` dataset
  - `RcppExports.R` — generated, do not edit
- `src/` — C++ (RcppEigen)
  - `moran_kernel.cpp` — tiled two-pass multi-bandwidth core
  - `gauss_kernel_weights.cpp` — sparse weight matrix from a Gaussian kernel
  - `sparse_dense_mult.cpp` — sparse-dense multiplication
  - `Makevars`, `Makevars.win` — OpenMP flags via `$(SHLIB_OPENMP_CXXFLAGS)`
  - `RcppExports.cpp` — generated, do not edit
- `tests/testthat/` — suite pinned to the reference implementation, with fixtures
- `data/`, `data-raw/` — packaged `ccRCCExample` dataset and its preparation script
- `process_docs/` — development notes from earlier design rounds. Treat them as historical: they
  predate the current compute core and inference code, and some describe behaviour that no longer
  exists. Build-ignored, so they do not ship in the package.
- `inst/validation/` — benchmarks and the approximation calibration harness
- `plans/` — implementation plans

## Bandwidth units — the most common error in this package

`ccRCCExample` coordinates are in **millimetres**: they span 1.02 × 1.02 mm with a median
nearest-neighbour distance of 0.0073 mm. The manuscript's bandwidth range of 10–250 is in
**micrometres**, i.e. `h = 0.01` to `0.25` in this dataset's units.

Using `h = 10` or `h = 25` here makes every off-diagonal Gaussian weight ≥ 0.9999 and the mean row sum
over 7119 of 7122 cells — the kernel becomes a constant matrix and Moran's I collapses toward zero.
Earlier validation in this repository did exactly that, which is why it appeared to pass while testing
nothing. Always sanity-check `mean(rowSums(K))` against `N` before trusting a result, and use
`suggest_bandwidths()` rather than hard-coding numbers.

## Dependencies

R >= 3.5, `Rcpp` (>= 1.0.13), `RcppEigen` (LinkingTo), `Matrix`, `methods`, `ggplot2`, `stats`,
`utils`, `parallel`. Suggested: `testthat`, `knitr`, `rmarkdown`, `mxfda`, `PearsonDS`.

## Main Functions

- `multiScaleAC(coords, gene1, gene2, h, perms, ...)` — main entry point; returns an S4
  `kernelcomm_result`
- `multiScaleAC_batch(data_list, ...)` — many samples, optionally in parallel, with checkpointing
- `compute_significance(result, ...)` — p-values, z-scores, null envelope, supremum test
- `adjust_pvalues(result, method)` — multiple-testing correction
- `export_mxfda(result)` — data frame for `mxfda::add_summary_function()`
- `suggest_bandwidths(coords, ...)`, `qc_neighborhoods(coords, h)` — bandwidth selection and QC
- `soupir_moran_multi(coords, gene1, gene2, h_vec, perms, seed)` — reference implementation
- `soupir_moran_bv(h, coords, gene1, gene2, perms, eps)` — single bandwidth, sparse
