# multiScaleAC 0.2.0

The package is renamed from `KernelCommunication` to `multiScaleAC`, matching the
companion manuscript. The repository keeps its old name, because that is the URL
the paper prints.

This release corrects results. Read the first section before reusing any earlier
output.

## Results that were wrong

**`kernel = "gaussian"` with `distance = "manhattan"` returned garbage.** For any
sample below 10,000 cells, that combination computed the kernel from uninitialised
memory: three identical calls returned `4.5e-05`, then `-5.3e-04`, then `NaN`,
against a correct value of `2.6e-02`, and one of them returned the same number for
every bandwidth. Any result from that configuration must be recomputed. The default
`distance = "euclidean"` was never affected.

**`alternative = "two.sided"` was one-sided.** It was the default, and its
two-sided branch was empty, so it returned a one-sided upper-tail p-value under a
two-sided label. Negative spatial cross-correlation could not be detected at any
setting. The default is now `"greater"`, which is honest and matches the
manuscript; `"two.sided"` now tests both tails.

**`alternative = "less"` was not a p-value.** It computed `1 - p`, which can be
exactly 0 and double-counts ties. It now counts the lower tail.

**Bandwidths in the documentation and examples were wrong by three orders of
magnitude.** `ccRCCExample` coordinates are in millimetres, spanning 1.02 mm, and
the shipped examples used `h = 10`–`250`. At those bandwidths every pairwise weight
exceeds 0.9999, the kernel is effectively constant, and Moran's I collapses toward
zero. The manuscript's 10–250 is in micrometres, i.e. `0.01`–`0.25` here. Use
`suggest_bandwidths()` and check `qc_neighborhoods()`, whose new `saturation`
column makes this visible.

**`eps` was silently ignored** in the code path taken for samples under 10,000
cells, so the sparse and dense paths disagreed on the row sums and therefore on
`I`. There is now one path.

## New

- **Exact analytical moments of the permutation null.** The statistic is a linear
  permutation statistic, so its first four moments and its full cross-bandwidth
  covariance have closed forms; they cost `O(N)` per bandwidth and are computed by
  default. `compute_significance(method = "moment")` uses them, and
  `method = "auto"` uses them only where diagnostics say a continuous
  approximation is trustworthy. `multiScaleAC(permute = "flagged")` draws
  permutations only where they are needed. The value is resolution rather than
  speed: they resolve below the `1/(perms + 1)` floor and agree with a
  499,999-permutation run to under 1%. The choice of approximation was made by
  measurement; see `inst/validation/results/README.md`, including two results that
  contradict what theory suggests.
- **The whole-curve supremum test** the manuscript defines, in `sup_p_value`. It
  needs no multiplicity correction across bandwidths.
- `estimator = c("add_one", "paper")` on `compute_significance()`, so the
  manuscript's `r/n` is available alongside the default `(r+1)/(n+1)`.
- `nan_policy` makes the zero-row-sum behaviour explicit. The default reproduces
  the reference implementation, which reports 0 for an entire affected bandwidth.
- `perms = 0` is now allowed.
- A `show` method, so printing a result no longer dumps the whole permutation
  matrix to the console.
- README, this file, and `vignette("multiScaleAC")`.
- An `"exponential"` kernel and a `"minkowski"` distance, both inside the same
  tiled core as the existing options, so neither can bypass the `eps` predicate or
  the zero-row-sum policy. `p = 1` and `p = 2` are routed to the Manhattan and
  Euclidean paths, which compute the same thing without calling `pow()`.
  `qc_neighborhoods()` accepts the new options, so a diagnostic always describes
  the model actually being fitted.

## Changed

- **Memory.** The compute core no longer allocates any N × N matrix. Working memory
  is `O(N × length(h))`: measured at N = 7122, peak RSS fell from 2230 MB to
  208 MB, and N = 20000 now runs in 217 MB where the previous code would have
  needed 11.9 GB. Slightly faster too, single-threaded.
- `batch_size` bounds the permutation buffer at `N × batch_size` instead of
  `N × perms`. It previously had no effect.
- `ci_lower` and `ci_upper` are renamed `null_lower` and `null_upper`. They are
  quantiles of the permutation distribution — an acceptance region for the null —
  and were named and plotted as a confidence interval for the observed statistic.
  The envelope level is now separate from the test threshold.
- Gene names come from the caller's expressions. They were always the literals
  `"gene1"` and `"gene2"`, in every plot title and every mxfda interaction name.
- `suggest_bandwidths()` returns genuinely different scales for `"local"`,
  `"medium"` and `"global"`. All three previously returned quantiles of the same
  nearest-neighbour distribution.
- `qc_neighborhoods()` accepts a bandwidth vector, reports the kernel row sum
  actually in use rather than a hard-radius count, and adds `saturation`. A vector
  previously recycled against the distance matrix and returned meaningless numbers.
- Neither helper builds an N × N distance matrix any more; both were unusable at
  the sample sizes this package targets.
- A custom `distance` with a built-in `kernel` now works. It always errored, which
  made the Minkowski example in the documentation throw.
- Custom kernels have no self-weight, matching the built-in path. The diagonal was
  not zeroed, so custom results were silently on a different scale.
- Batch checkpointing works under parallel processing, and resuming a completed
  batch no longer reprocesses it in reverse. Columns are validated before any
  computation starts.
- `plot_batch_comparison()` no longer errors on its own default argument.
- `show_permutations`, `cluster_rows` and `sort_by = "sample_id"` are implemented;
  they were documented and inert. `alpha` in `plot_scale_curve()` is renamed
  `ribbon_alpha`, and `ci_level` is removed, having had no effect.
- The z-score plot no longer draws ±1.96 lines. Its envelope comes from the
  permutation distribution, because a normal band contradicts the p-values in the
  same object.
- `export_mxfda()` always emits `p_adjusted`, so batch exports bind. It also now
  documents that the bandwidth occupies mxfda's radius column and that these are
  different quantities.
- `large_mode_n` is deprecated and warns. It described a threshold hard-coded in
  C++ and was never read.

## Removed

- The separate sparse compute path. At the default `eps = 0` it was slower and used
  more memory than the dense path while disagreeing with it numerically.
- `lag.sparse` is gone, with no shim: it was only ever registered as an S3 method
  for a class that `Matrix` objects do not have, so it was unreachable.

## Testing and documentation integrity

The suite went from zero tests to 407, pinned to `soupir_moran_multi()`, which is
never modified and is hash-guarded by a test.

Six documents in this repository asserted results that do not hold, including one
that declared the reference implementation "mathematically invalid" with a "73%
relative error" and advised noting the correction in publications. It was wrong.
Each is marked retracted in place, with the reasoning, in `docs/RETRACTIONS.md`.
