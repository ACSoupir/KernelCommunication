# Analytical permutation moments: what was measured

Regenerate with `Rscript inst/validation/calibrate_moments.R`. Raw table in `calibration.csv`.

## What is exact and what is approximate

**Exact.** Because the row-normalised weight matrix satisfies `W %*% 1 = 1` and both genes are
z-scored, the statistic is a *linear permutation statistic*: `I_pi(h) = sum_i a_i b_pi(i) / S`, with
`a = v_h - mean(v_h)` and `b = y`. Its first four moments and its full cross-bandwidth covariance
have closed forms in the power sums of `a` and `b`. Nothing is estimated.

Checked against 20,000 permutations on three datasets:

| quantity | agreement |
|---|---|
| variance | ratio 0.983 – 1.010 |
| skewness | 1.217 exact vs 1.210 empirical; 0.661 vs 0.636; 0.184 vs 0.174 |
| excess kurtosis | ratio ~1.0 wherever skewness is large enough to resolve it |
| mean | empirical mean 6.3e-05, below its own Monte Carlo error of 1.3e-04, consistent with exactly 0 |
| cross-bandwidth covariance | within 0.8% |

**Approximate.** Turning those moments into a tail probability. Four methods are provided; the
default was chosen by the measurement below, not by theory.

## Calibration of the tail approximations

200,000 permutations per configuration; four configurations from benign to adversarial; scored at
target upper-tail probabilities down to 1e-4, because an FDR screen over the ~1100
ligand-receptor pairs per sample in the companion manuscript operates there. The figure reported is
approximate / empirical. **Below 1 is anti-conservative**, which manufactures false positives.

Over the bandwidths the diagnostics do *not* flag, i.e. those that would actually use moment
inference:

| method | median | range | anti-conservative (< 0.5) |
|---|---|---|---|
| **gamma** (default) | 1.05 | [0.86, 3.11] | **0 of 24** |
| edgeworth | 1.01 | [0.78, 3.11] | 0 of 24 |
| saddlepoint | 1.00 | [0.28, 1.53] | 2 of 24 |
| normal | 1.01 | [0.12, 5.47] | 2 of 24 |

**The three-moment shifted gamma is the default.** It is the only method with no anti-conservative
failures and a bounded worst case; where it errs, it errs conservatively.

### Two results that contradict what theory suggests

**1. The saddlepoint measured worse than the gamma.** It is the best-founded method here — a double
saddlepoint on the conditional representation of Robinson (1980) and Booth & Butler (1990) — and it
is excellent under strong skew, where it stayed within 0.91–1.02 against the gamma's 0.83–1.20. But
on the real `ccRCCExample` data at h = 0.08, a near-symmetric case, it returned 0.59, 0.38 and 0.28
at targets of 1e-2, 1e-3 and 1e-4, understating the tail. That discrepancy is **unexplained**;
candidates are the 128-atom representation of a continuous gene and the conditioning step. It ships
as an option, marked experimental, and is not the default.

**2. The normal approximation fails on real data, not just contrived examples.** On
`ccRCCExample` at h = 0.02, where the exact skewness is 1.91, it returned ratios of 0.01 at a target
of 1e-2 and **0.00** at 1e-3 and 1e-4. Anyone converting the `z_scores` slot with `pnorm()` on data
like this would generate false positives at exactly the threshold they were screening at. This is
why `compute_significance()` documents that the slot is an effect size only.

Edgeworth is included for comparison, not recommended: it is not monotone, and on one configuration
it returned 32.52 and then 0.00 at adjacent targets.

## The diagnostic and the fallback

`leverage` is the largest single cell's share of the null variance. It falls monotonically with
bandwidth — 0.156 at h = 0.01 down to 0.005 at h = 0.25 on the bundled data — because small
bandwidths concentrate the statistic on few cells, which is exactly the region the manuscript reports
as unstable. Bandwidths with leverage above 0.10, or excess kurtosis above 20, are flagged;
`permute = "flagged"` draws permutations only there, and `method = "auto"` uses them only there.

The 0.10 threshold is empirical. At 0.05 it flagged a configuration where every method was accurate
to within 60%; 0.10 separates the configurations where the gamma stayed within 6% of the truth from
those where it did not.

## What this actually buys

Measured on the full `ccRCCExample` (N = 7122), six bandwidths, single-threaded:

| mode | wall clock | peak RSS | permutations stored |
|---|---|---|---|
| `permute = "all"`, 999 perms | 1.99 s | 282 MB | 999 |
| `permute = "flagged"`, 999 perms | 1.79 s | 166 MB | 0 (none flagged at this N) |
| `permute = "none"` | 1.78 s | 175 MB | 0 |

**The gain is resolution, not speed.** For a single gene pair the kernel sweep dominates, so removing
permutations saves 1.7x on memory and only ~10% on time. What it does buy is a p-value below the
permutation floor. With 999 permutations the smallest reportable p-value is 1/1000, so p = 0.001 and
p = 1e-6 are indistinguishable; reaching 1e-5 resolution needs ~100,000 permutations. The moments
give it for nothing, and they are correct there:

| bandwidth | moment | 999 perms | 9,999 | 99,999 | 499,999 |
|---|---|---|---|---|---|
| 0.04 | 0.4692 | 0.458 | 0.4757 | 0.4701 | 0.4701 |
| 0.08 | 0.2097 | 0.190 | 0.2133 | 0.2103 | 0.2103 |
| 0.15 | 0.02473 | 0.026 | 0.0267 | 0.02534 | 0.0249 |

The moment values agree with a half-million-permutation run to under 1%.

Note also that at N = 7122 no bandwidth was flagged: leverage falls as cells are added, so the
approximation gets *better* on larger samples, which is the opposite of how permutation cost scales.

## Honest limitations

- The gamma matches three moments. The exact fourth moment is computed and reported but not used by
  the default method; a four-moment fit was not shown to beat it.
- The simulation-based whole-curve supremum test that the exact cross-bandwidth covariance makes
  possible is **not implemented**. It would need a Gaussian copula for the joint null, and that
  assumption has not been validated here. The permutation supremum test in
  `compute_significance()` is used instead, and it requires every bandwidth to be permuted.
- Calibration used four configurations. It is enough to rule out the normal and Edgeworth methods
  decisively and to establish the gamma as safe on these data; it is not a broad survey.
- Permutations remain the default for `compute_significance()`, because they require no
  approximation at all.
