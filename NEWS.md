# multiScaleAC 0.2.0

First release of the `multiScaleAC()` interface. The package is renamed from
`KernelCommunication` to `multiScaleAC`, matching the companion manuscript; the
repository keeps its old name, because that is the URL the paper prints.

Earlier releases contained only `soupir_moran_multi()` and `soupir_moran_bv()`.
Both are unchanged, remain exported, and are not deprecated — `soupir_moran_multi()`
is the reference implementation the whole test suite is pinned to. So results
computed with previous versions are unaffected by this release.

## New: the multiScaleAC interface

- `multiScaleAC()` computes the bivariate Moran's I curve over a bandwidth vector
  and returns an S4 `kernelcomm_result`. Working memory is `O(N × length(h))`, with
  no N × N matrix at any point: 256 MB at N = 7122 against 2225 MB for the
  reference, and N = 20000 runs in 440 MB. Results match the reference exactly.
- `multiScaleAC_batch()` runs many samples, optionally forked across cores, with
  per-sample checkpoint files so an interrupted batch resumes where it stopped.
- `compute_significance()` and `adjust_pvalues()` for inference, kept separate from
  computation so one curve can be tested several ways.
- `suggest_bandwidths()` and `qc_neighborhoods()` for choosing and checking
  bandwidths, and `validate_inputs()` for pre-flight checks.
- `plot()`, `summary()`, `print()`, `show()`, plus `plot_scale_curve()`,
  `plot_permutation_distribution()`, `plot_batch_comparison()` and
  `plot_significance_summary()`.
- `export_mxfda()` and `export_mxfda_batch()` for handing curves to functional data
  analysis, as the manuscript does with FPCA.

## Inference

- Permutation p-values, permuting only the second gene. `alternative` defaults to
  `"greater"`, the manuscript's one-tailed upper test; `"less"` and `"two.sided"`
  are available. `estimator` offers `(r+1)/(n+1)` (default, cannot return zero) or
  the manuscript's `r/n`.
- **The whole-curve supremum test** the manuscript defines, in `sup_p_value`. It
  needs no multiplicity correction across bandwidths, which matters because
  bandwidths of one gene pair are strongly dependent.
- **Exact analytical moments of the permutation null.** The statistic is a linear
  permutation statistic, so the first four moments of its null and the full
  cross-bandwidth covariance have closed forms. They cost `O(N)` per bandwidth and
  are computed by default. `compute_significance(method = "moment")` uses them;
  `method = "auto"` uses them only where diagnostics say a continuous approximation
  is trustworthy; `multiScaleAC(permute = "flagged")` draws permutations only where
  they are needed.

  The gain is resolution rather than speed. With 999 permutations the smallest
  reportable p-value is 1/1000, so `p = 0.001` and `p = 1e-6` are
  indistinguishable, and resolving `1e-5` needs about 100,000 permutations. The
  moments have no such floor and agree with a 499,999-permutation run to under 1%.

  The moments are exact; mapping them to a tail probability is not. The default
  mapping was chosen by measurement against 200,000 permutations rather than by
  theory, and `inst/validation/results/README.md` records the comparison, including
  two outcomes that contradict what theory would predict.
- `null_lower` and `null_upper` give the null envelope: quantiles of the
  permutation distribution, showing where the statistic falls under complete
  spatial randomness. They are not a confidence interval for the observed value,
  and the envelope level is set separately from the test threshold.
- `z_scores` is an effect size only. Do not convert it with `pnorm()`: the
  permutation null is skewed whenever the expression data is, and on the bundled
  example data a normal approximation is anti-conservative by two orders of
  magnitude in the tail.

## Bandwidth units

Bandwidths are in the units of `coords`, and this is the easiest way to get a
meaningless result. `ccRCCExample` is in **millimetres**, spanning 1.02 mm with a
median nearest-neighbour distance of 0.0073 mm, so useful bandwidths are about 0.01
to 0.25 — which is the manuscript's 10–250 **micrometres** expressed in this
dataset's units. Passing `h = 25` drives every pairwise weight above 0.9999, making
the kernel effectively constant and collapsing `I(h)` toward zero for every gene
pair.

`qc_neighborhoods()` reports a `saturation` column for exactly this, and a test
fails if any fixture uses a bandwidth in that regime.

## Kernels and distances

`"gaussian"` (default, and what the manuscript uses), `"uniform"` and
`"exponential"` kernels; `"euclidean"`, `"manhattan"` and `"minkowski"` distances.
All run inside the same tiled core, so none can bypass the `eps` predicate or the
zero-row-sum policy. User-supplied kernel and distance functions are still
supported and take an R fallback path.

`eps` defaults to 0, meaning no thresholding, which is the point of the method:
computing the statistic without discarding low weights is what distinguishes it
from approaches that filter low kernel density estimates first. Raising it is
supported but changes the estimand.

`nan_policy` makes the zero-row-sum behaviour explicit. The default reproduces the
reference implementation, which reports 0 for an entire affected bandwidth.

## Also in this release

- `lag.sparse()` is renamed `lag_sparse()`. The old name was only ever registered
  as an S3 method for a class that `Matrix` objects do not have, so it was never
  callable and no shim is provided.
- Compiled objects and shared libraries are no longer tracked in git, and the
  source tarball no longer ships development notes, the manuscript PDF or
  build artefacts.
- `man/rcpp_hello_world.Rd`, which documented a function that does not exist, is
  removed.
- README, this file, and `vignette("multiScaleAC")`.

## Testing

451 tests, pinned to `soupir_moran_multi()`, which is never modified and is
hash-guarded by a test so it cannot drift. Fixtures are generated from it alone,
at unit-correct bandwidths, and cover the degenerate cases — isolated cells,
zero-variance genes, saturated and empty kernels — as well as the ordinary ones.
`R CMD check --as-cran` runs the examples, the tests and the vignette.
