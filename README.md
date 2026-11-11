# multiScaleAC

Bivariate Moran's I as a function of Gaussian kernel bandwidth: a functional curve
`I(h)` over spatial scales, rather than a single value at one arbitrarily chosen
scale.

Intended for ligand–receptor interaction analysis in spatial transcriptomics.
Weights come directly from the kernel, row standardised, with **no minimum weight
threshold**, so no cells or spots are dropped from the weight matrix. That is the
main methodological difference from approaches that filter low kernel density
estimates before computing the statistic.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ACSoupir/multiScaleAC")
```

Requires a C++ toolchain (`Rcpp`, `RcppEigen`). OpenMP is used when available;
builds without it, including stock Apple clang, work single-threaded.

## Quick start

```r
library(multiScaleAC)
data(ccRCCExample)

coords <- as.matrix(ccRCCExample[, c("x", "y")])

# Bandwidths are in COORDINATE UNITS. This dataset is in millimetres, so these
# come back in millimetres. Do not guess -- ask.
h <- suggest_bandwidths(coords, n_bandwidths = 8)

# Sanity-check before interpreting anything. `saturation` near 1 means the kernel
# is effectively a constant matrix and no result at that bandwidth is meaningful.
qc_neighborhoods(coords, h)

res <- multiScaleAC(
  coords, ccRCCExample$CD274, ccRCCExample$PDCD1,
  h = h, perms = 999
)
res <- compute_significance(res)

summary(res)
plot(res)

# Whole-curve significance in one number, needing no multiplicity correction
res@sup_p_value
```

## The thing that goes wrong most often

**Bandwidths are in the units of your coordinates.** Get the scale wrong and the
output is meaningless in a way that looks like a real answer.

`ccRCCExample` spans about 1.02 mm, with a median nearest-neighbour distance of
0.0073 mm:

```r
qc_neighborhoods(coords, h = c(0.02, 25))
#>       h mean_neighbors  ...  saturation
#> 1  0.02          20.31  ...      0.0029
#> 2 25.00        7119.29  ...      0.9996
```

At `h = 25` every pairwise weight is ≥ 0.9999 and each cell is effectively
everyone's neighbour. The kernel is constant, so `I(h)` collapses toward zero for
every gene pair — about `2e-06` here, against `5e-02` at a sensible bandwidth.
Earlier validation in this repository ran entirely inside that regime and reported
agreement to `1e-18`, which is agreement on an arithmetic identity rather than
evidence of anything.

The manuscript's range of 10–250 is in **micrometres**, i.e. `0.01`–`0.25` here.
Use `suggest_bandwidths()`, then check `qc_neighborhoods()`.

## Inference

Two routes, both available on the same result.

**Permutations** (the default). Only the second gene is permuted. `alternative`
defaults to `"greater"`, the manuscript's one-tailed upper test; `"less"` and a
genuine `"two.sided"` are available.

**Exact analytical moments.** Because the weight matrix is row standardised and
both genes are z-scored, the statistic is a *linear permutation statistic*, so the
first four moments of its null and the full cross-bandwidth covariance have closed
forms. They cost `O(N)` per bandwidth and are computed by default.

```r
# No permutations at all
res <- multiScaleAC(coords, g1, g2, h = h, permute = "none")
res <- compute_significance(res, method = "moment")

# Or permute only where the diagnostics say an approximation cannot be trusted
res <- multiScaleAC(coords, g1, g2, h = h, perms = 999, permute = "flagged")
res <- compute_significance(res, method = "auto")
res@p_method   # records which method produced each p-value
```

**What this buys is resolution, not speed.** With 999 permutations the smallest
reportable p-value is 1/1000, so `p = 0.001` and `p = 1e-6` are
indistinguishable, and resolving `1e-5` needs roughly 100,000 permutations. The
moments give it for free, and they are correct there: they agree with a
499,999-permutation run to under 1%. For a single gene pair the kernel sweep
dominates, so dropping permutations saves about 1.7× on memory and only ~10% of
wall clock.

Do **not** run `pnorm()` on the `z_scores` slot. It is an effect size. The
permutation null is skewed whenever the expression data is, and on this package's
own example data a normal approximation is 100× anti-conservative in the tail —
the direction that manufactures false positives. The measurements, and the
reasoning behind the default approximation, are in
[`inst/validation/results/README.md`](inst/validation/results/README.md).

## Downstream: functional data analysis

Each gene pair yields a curve, which is what makes functional data analysis
natural. `export_mxfda()` produces a table for
`mxfda::add_summary_function()`, as the manuscript uses for FPCA. Note that the
bandwidth occupies mxfda's radius column; `?export_mxfda` explains why those are
not the same quantity.

## Many samples

```r
fovs <- split(ccRCCExample, ccRCCExample$fov)
results <- multiScaleAC_batch(
  fovs, h_vec = h, perms = 999,
  n_cores = 4, checkpoint_dir = "ckpt"    # resumable
)
plot_batch_comparison(results)
```

The manuscript applies FDR across *gene pairs*. To reproduce that, collect one
p-value per pair — `sup_p_value` is the natural choice — and correct across the
collection. `adjust_pvalues()` corrects across bandwidths *within* one result,
which are strongly dependent and rarely what you want.

## Reference implementation

`soupir_moran_multi()` is the reference implementation. It is permanently
supported, is not deprecated, and is the oracle the test suite is pinned to: every
fixture is generated from it, and a test detects any edit to it.

## Documentation

- `vignette("multiScaleAC")` — choosing bandwidths, reading the curve, inference
- [`inst/validation/results/README.md`](inst/validation/results/README.md) — what
  was measured, and what failed

## Citation

Soupir AC, Hayes MT, Manley BJ, Wang X, Peres LC, Wrobel J, Fridley BL.
*multiScaleAC: Cell-Cell interaction with Moran's I as a function of kernel
bandwidth.* bioRxiv 2026.07.14.738306.
<https://doi.org/10.64898/2026.07.14.738306>

## License

GPL (>= 2)
