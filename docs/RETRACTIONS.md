# Retracted validation and status documents

Several documents in this repository asserted results that do not hold. They have been removed rather
than corrected in place, because each was written as a settled finding and would otherwise be
rediscovered as fact. This file records what each claimed, why it is wrong, and what is actually true,
so the same conclusions are not reached a second time.

The retracted files were, at commit `9b14cca`:

- `docs/validation/issue_001_legacy_normalization_bug.md`
- `docs/validation/validation_report.md`
- `docs/validation/issues.md` (issue #004 only; the rest is sound)
- `docs/planning/analytical_variance_test_report.md`
- `IMPLEMENTATION_COMPLETE.md`
- `validation_results.txt`

---

## 1. `soupir_moran_multi()` was never wrong

**Claimed** (`issue_001_legacy_normalization_bug.md`, severity "Critical"): the reference
implementation "uses an incorrect row-normalization formula that produces mathematically invalid
Moran's I values", with a "73% relative error". It advised users to "Trust the new values", "Re-run
analyses", and "Note in publications that corrected formula was used".

**Why it is wrong.** The document's own algebra disproves it. It writes the two expressions as

```
legacy  = Σ_i Σ_j g2[j] · K[i,j] · g1[i] / r[i]        (its line 76)
correct = Σ_i Σ_j g1[i] · K[i,j] · g2[j] / r[i]        (its line 82)
```

These are the same expression — the summand is an identical product of identical factors, merely
written in a different order. No symmetry argument is even required. The document then applies a
relabelling of the bound dummy indices to one side only, compares the relabelled form against the
original, and concludes they are "NOT equivalent unless `r_sums[i] = r_sums[j]`". Relabelling bound
indices cannot change the value of a sum, so that conclusion does not follow.

The "Empirical Demonstration" that produced the 73% figure uses

```r
W_correct <- t(t(K) / r_sums)
```

In R, `t(K) / r_sums` divides *row* `i` of `t(K)` by `r_sums[i]`, giving `K[j,i] / r_sums[i]`;
transposing yields `W_correct[i,j] = K[i,j] / r_sums[j]`. That is **column**-normalized, not
row-normalized. The 73% is the gap between row- and column-normalization, not an error in the
reference implementation. The snippet also references `g1_s` and `g2_s`, which it never defines
(it defines `g1` and `g2`), so it cannot have run as written.

**What is actually true.** Because `W` is row-normalized (`W %*% 1 = 1`) and `K` is symmetric,

```
v_h = K %*% (x / r) = t(W) %*% x        so        I = <v_h, y> / S = x' W y / S
```

Verified numerically: `max(abs(v_h - t(W) %*% x)) = 2.2e-16` and
`abs(I_vh - x'Wy/S) = 2.4e-17`. The reference formulation is the standard row-standardized
bivariate Moran's I. It is exact, and it is also the cheaper of the two forms, since it needs only a
matrix–vector product and never materializes `W`.

Two further errors in the same document: its "Verification" section asserts that under the null
`E[I]` should equal `-1/(N-1)`. That is the *univariate* Moran's I result. Under the permutation
scheme actually used here — permuting `y` only, with both genes centred — `E[I] = 0` exactly, because
`mean(v_h) = 0` follows from `W %*% 1 = 1` and `sum(x) = 0`.

## 2. The two validation documents contradict each other

`validation_report.md` states "Observed values match exactly (within 1e-14)" and
"confirming mathematical equivalence". `issue_001` states the same comparison shows a 73% error and
that the reference is "mathematically invalid". Both cannot hold. The match claim is the correct one.

`validation_report.md` additionally claims "Analytical inference methods validated". No analytical
inference code was active: the only such function was commented out at
`R/multi_scale_ac.R:467-484`.

## 3. `issues.md` #004 describes code that does not exist

**Claimed:** permutation values differ between the reference and the new implementation because of
"Different random number generators and shuffle algorithms — Legacy: R's Mersenne Twister +
`sample()`; New: C++ `std::rand()` + Fisher-Yates shuffle", and this is "EXPECTED BEHAVIOR".

**Why it is wrong.** The C++ code contains no random number generator at all. The permutation matrix
is built in R (`R/multi_scale_ac.R:187-188`) with the same `set.seed` / `replicate` / `sample`
sequence as the reference, and the one matrix is reused across all bandwidths in both. The
permutations therefore agree exactly — as `validation_results.txt` itself reported
("Permutation values match: YES"), contradicting this issue in the same repository.

Issues #001, #002 and #003 in that file are sound and are preserved in `docs/validation/issues.md`.

## 4. The analytical variance formula is right; the claims made about it were not

`analytical_variance_test_report.md` gives `E[I] = 0` and `Var(I) = ||W'x||² / (N-1)²`. **This
formula is correct**, and it is retained and generalized in the current implementation — it is the
special case of the exact second moment

```
Var(I) = A2 · B2 / ((n-1) · S²),   A2 = Σ(v_h - mean(v_h))²,   B2 = Σy²
```

when both genes are z-scored, so `B2 = S = n-1`. Nothing here retracts the mathematics.

What is retracted is the surrounding evidence. The report is titled a test report, marks every case
"✓ PASSED", and concludes "READY FOR PRODUCTION" — for code that is commented out and therefore
never executed. Its numbers cannot have been produced by the state of the repository it describes.
Its accuracy tables are reported at bandwidths 25–250 on `ccRCCExample`, which is a degenerate
regime (see §6). And its framing is inverted: it attributes "~2% Monte Carlo error" to the analytical
value, when the analytical moments are exact and it is the empirical estimate that carries Monte
Carlo error.

## 5. There were no tests

`IMPLEMENTATION_COMPLETE.md` claimed "Comprehensive Test Suite: 17/17 PASSING ✓" and
"- [x] Unit tests (17/17 passing)", and enumerated seventeen named tests. No `tests/` directory
existed anywhere in the repository or in its history, and no test file existed under any name. The
document also declared "Status: ✅ COMPLETE - Ready for Production Use" for a state containing an
uninitialized-memory read in the default dense compute path.

A real suite now lives in `tests/testthat/`.

## 6. `validation_results.txt` measured a degenerate configuration

**Claimed:** "OVERALL: PASSED", observed Moran's I matching to `1.39e-18`, permutations matching to
`1.34e-18`, and a `71.65x` speedup, on `ccRCCExample` (N = 7122) at `h = 10, 20, 40, 80`.

**Why it establishes nothing.** The coordinates in `ccRCCExample` are in **millimetres**: they span
1.02 × 1.02 mm, with a median nearest-neighbour distance of 0.0073 mm. At `h = 10` — let alone
`h = 80` — every off-diagonal Gaussian weight is ≥ 0.9999 and the mean row sum is over 7119 of the
7122 cells. The kernel is a constant matrix, so Moran's I collapses toward zero; the reported
observed values are ~2e-6. Agreement to 1e-18 between two implementations of a constant kernel is
agreement on an arithmetic identity, and would hold whether or not either implementation was correct.

The same unit error affects every benchmark and documentation example that shipped in that state
(`h = 10..80` and `h = 25..250`). The manuscript's bandwidth range of 10–250 is in **micrometres**,
which in this dataset's units is `h = 0.01` to `0.25`. Current fixtures, examples and the vignette use
that corrected range, and a test now fails if any fixture uses a bandwidth where the mean row sum
exceeds most of the sample.

The `71.65x` speedup is also not a like-for-like comparison: it timed four separate
`soupir_moran_bv()` calls, each rebuilding the kernel, against one `multiScaleAC()` call that builds
the distance matrix once and reuses it across all four bandwidths. That is an amortization gain, not
a C++-versus-R gain.
