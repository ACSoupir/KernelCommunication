# multiScaleAC — repository remediation and roadmap

## Context

This repo implements the method in the companion preprint *"multiScaleAC: Cell-Cell interaction with
Moran's I as a function of kernel bandwidth"* (Soupir et al., bioRxiv 2026.07.14.738306). A previous
round of development built a v2 API on top of the original ground-truth function. That work has genuine value — an S4
result class, plotting, an mxfda bridge, a batch driver — but the review found it is not trustworthy as
it stands, and **none of it is committed to git**.

Verified defects (each confirmed by reading the code, not inferred):

- **Silent wrong answers.** `src/moran_multi_cpp.cpp:94` reads `dist_sq` unconditionally in the dense
  Gaussian branch, but `dist_sq` is only filled when `distance == "euclidean"` (`:43-53`). So
  `kernel="gaussian", distance="manhattan"` at N < 10000 computes the kernel from uninitialized heap.
- **Wrong p-values.** `R/significance.R:52-65`: the default `alternative="two.sided"` computes a
  one-sided upper-tail p — the two-sided branch is empty. `alternative="less"` uses `1 - p`, which is
  not a p-value (it can be exactly 0, and it double-counts ties).
- **The validation evidence base is vacuous.** `ccRCCExample` coordinates are in **millimetres**,
  spanning 1.02 × 1.02 mm with a median nearest-neighbour distance of **0.0073 mm**. Every benchmark,
  validation script and doc example in the repo uses `h = 10, 20, 40, 80` or `h = 25..250`. At those
  bandwidths every off-diagonal weight is ≥ 0.9999 and the mean row sum is ~7120 of 7122 cells — the
  kernel is a constant matrix and Moran's I collapses. That is why `validation_results.txt` reports
  observed I ≈ 2e-6 and agreement to 1e-18: it is agreement on a degenerate all-ones kernel.
  **The "OVERALL: PASSED" result tests nothing.** The paper's h = 10–250 is in µm, i.e. 0.01–0.25 mm here.
- **No tests exist at all** (no `tests/` directory), while `IMPLEMENTATION_COMPLETE.md` claims
  "17/17 PASSING". All three root benchmark scripts call `$` on what is now an S4 object and cannot run.
- **`docs/validation/issue_001_legacy_normalization_bug.md` is provably wrong.** It calls the
  ground-truth formula "mathematically invalid" with "73% relative error". Verified numerically: because
  `W` is row-normalized, `v_h = K %*% (x/r)` equals `t(W) %*% x` to 2.2e-16, so `I = x'Wy/S` exactly.
  The doc's "correct" reference `t(t(K)/r_sums)` is *column*-normalized — that is the fabricated 73%.
- Plus: ~1.6 GB peak memory at N=7122 (4 unconditional N×N allocations) rising to >4.4 GB in the
  "large" path; `eps` silently ignored in dense mode; three different NaN policies; gene names never
  captured; `plot_batch_comparison()` fails on its own default; `parallel` used but undeclared;
  4.53 MB of `.o`/`.so` tracked in git; ~6.5 MB of cruft shipping in the tarball.

**The opportunity.** Because `W` is row-normalized (`W1 = 1`) and both genes are z-scored, the statistic
is exactly a **linear permutation statistic** — which has closed-form moments. This was derived and
numerically verified against 2×10⁶ permutations in a strongly-skewed regime (ratios 1.002 / 1.006 /
1.008 for the 2nd/3rd/4th moments; cross-bandwidth covariance 0.999). It is precisely the "future
opportunity ... identification of analytical mean and variance of the distributions, eliminating the
requirement of permutations" that the manuscript's Limitations section names.

**Outcome sought:** a well-documented, user-friendly, correctly-named `multiScaleAC` package whose
numbers are pinned to the ground-truth function by an actual test suite, which scales to real CosMx
sample sizes at `eps = 0`, and which can produce calibrated significance with far fewer permutations.

## Decisions taken

| Decision | Choice |
|---|---|
| Package name | R package → `multiScaleAC`. Local directory and GitHub repo stay `KernelCommunication` so the manuscript's printed URL resolves. |
| Prior work | Keep the API surface; rewrite the compute core; fix every confirmed bug; delete the false docs; add real tests. |
| `eps` | **Kept** as an opt-in parameter, default `0`, honored identically everywhere. |
| Execution | All phases, end to end, one commit per phase. |

Two constraints hold throughout: **`R/soupir_moran_multi.R` is never modified** (it is the regression
oracle, and `soupir_moran_multi()` stays exported permanently), and **nothing is pushed** — all work is
local commits. No file, comment, doc or commit message mentions Claude, AI, agents or Anthropic.

### One refinement to the `eps` decision, flagged explicitly

You chose "keep `eps` as opt-in, fix the divergence". I will keep `eps` exactly as you asked, but I
recommend fixing the divergence **by construction rather than by discipline**: collapse the two
multi-bandwidth compute paths (dense for N<10000, sparse for N≥10000) into a single tiled path that
applies the `w > eps` predicate once. Reasons: at the default `eps = 0` the "sparse" path is strictly
worse on both axes (~3.2 GB of triplets *plus* the two dense distance matrices), and as long as two
paths exist, "eps behaves identically in both" is a property maintained by vigilance — which is exactly
how the current divergence arose. One path makes the bug inexpressible. The standalone sparse
capability users can already reach (`gaussian_kernel_weights()`, `sparse_dense_mult()`,
`soupir_moran_bv()`) is untouched. If you would rather keep two paths, say so and I will instead add a
cross-path `identical()` test over an `eps` grid.

---

## Phase 1 — Checkpoint and hygiene

Nothing of the current package is in git; a single `rm -rf` loses it. Commit first, clean second.

- **Commit A (safety):** `git add` the untracked v2 sources exactly as they are — `R/{classes,helpers,
  significance,methods,plot,export_mxfda,multi_scale_ac,KernelCommunication-package}.R`,
  `src/moran_multi_cpp.cpp`, `src/Makevars`, the new `man/*.Rd`, plus the staged modifications. This is
  a known-imperfect baseline, and that is the point: every later phase is diffable against it.
- **Commit B (hygiene):**
  - `git rm --cached` the four tracked binaries (`src/KernelCommunication.so`, `RcppExports.o`,
    `gauss_kernel_weights.o`, `sparse_dense_mult.o`); add `*.o`, `*.so`, `*.dll`, `src/symbols.rds`,
    `.DS_Store`, `*.bak` to `.gitignore`. Untrack `KernelCommunication.Rproj` (already gitignored).
  - Delete `R/plot.R.bak` (a stale 424-line duplicate of `plot.R`, and an invalid file in `R/`).
  - `.Rbuildignore`: add the three root benchmark scripts, `validation_results.txt`, `AGENTS.md`,
    `plans/`, `*.pdf`, and `^docs$` — ~6.5 MB that currently ships.
  - **Retract, don't preserve,** the four documents whose central claims are false:
    `docs/validation/issue_001_legacy_normalization_bug.md`, `docs/validation/validation_report.md`
    (asserts the opposite of issue_001), `docs/planning/analytical_variance_test_report.md` (documents
    code that is commented out at `multi_scale_ac.R:467-484`), `IMPLEMENTATION_COMPLETE.md` ("17/17
    PASSING" with no tests). Replace with one short `docs/RETRACTIONS.md` recording what each claimed,
    why it was wrong, and the correct result — so the error cannot be rediscovered as fact. Delete
    `validation_results.txt` and the three unrunnable benchmark scripts; their role is taken over by
    `tests/` and a real benchmark in Phase 12.
  - Update `AGENTS.md`: it predates the refactor, misnames `src/gauss_kernel_weights.cpp`, and never
    mentions `multiScaleAC()`.

*Verify:* `git status` clean but for intended deletions; `git ls-files | grep -E '\.(o|so)$'` empty.
*Commits:* `checkpoint work-in-progress modular sources` / `remove build artifacts and retract incorrect validation notes`

## Phase 2 — Toolchain and the golden reference

`RcppEigen`, `roxygen2`, `devtools` and `mxfda` are not installed, so **no C++ in this package can
currently compile**. Bootstrap them (R 4.6.1; `PearsonDS` too, for Phase 8), then confirm
`R CMD INSTALL .` succeeds on the Phase 1 baseline — establishing that later breakage is ours.

Then build the regression oracle, in `tests/testthat/fixtures/`:

- `make_golden.R` — generates fixtures **using only `soupir_moran_multi()`**, which is pure base R and
  therefore needs no build and cannot be invalidated by the rename or the core rewrite.
- Cases: a small synthetic set (n = 200/300, several seeds) plus a subsample of `ccRCCExample`, each at
  **unit-correct bandwidths**. Define one canonical vector and reuse it in fixtures, examples and the
  vignette so the three cannot drift: `H_MM <- c(0.01, 0.02, 0.04, 0.08, 0.15, 0.25)` — the paper's
  10–250 µm expressed in this dataset's millimetres. At `h = 0.02` the mean row sum is ~20 effective
  neighbours; at `h = 0.25` it is a substantial fraction of the core. Include a degenerate case
  (isolated cell under a uniform kernel) to pin the NaN semantics, and a zero-variance gene.
- Store `h`, `observed`, and the full `permutations` matrix per case.

*Verify:* fixtures regenerate byte-identically on a second run (RNG is seeded); baseline install succeeds.
*Commit:* `add golden reference fixtures generated from the reference implementation`

## Phase 3 — Rename to `multiScaleAC`

Placed **before** the core rewrite deliberately: the rename forces `Rcpp::compileAttributes()` to
regenerate every `_KernelCommunication_*` C entry point, and the core rewrite forces the same. Doing the
rename first means that churn happens once, against a baseline already pinned by Phase 2 fixtures,
instead of twice with a rewritten core confounding the diff.

Ordered mechanics (order matters — `compileAttributes()` reads `DESCRIPTION`):

1. `DESCRIPTION`: `Package: multiScaleAC`. While here: `Version: 0.2.0` (1.0 overstates it),
   `Date` → current, drop the redundant `Author:`/`Maintainer:` beside `Authors@R`, remove the
   `ORCID = "XX"` placeholder, add `URL`/`BugReports` pointing at the unchanged repo, add
   `Suggests: testthat (>= 3.0.0), knitr, rmarkdown, mxfda, PearsonDS`, and add the missing
   `Imports: parallel` (currently an `R CMD check` ERROR at `multi_scale_ac.R:720`).
2. `git mv R/KernelCommunication-package.R R/multiScaleAC-package.R`; update `@keywords internal` doc.
3. Run `Rcpp::compileAttributes()` then `roxygen2::roxygenise()` — this regenerates `NAMESPACE`
   (`useDynLib(multiScaleAC, .registration = TRUE)`), `src/RcppExports.cpp`, `R/RcppExports.R` and the
   package `.Rd`. Never hand-edit these (per `AGENTS.md`).
4. Delete the stale `man/KernelCommunication-package.Rd`, `man/rcpp_hello_world.Rd` (documents a
   function that no longer exists) and `man/lag.sparse.Rd`; delete stale `src/symbols.rds`.
5. Rename `KernelCommunication.Rproj` → `multiScaleAC.Rproj`.
6. `grep -rn KernelCommunication` and reconcile every remaining hit: keep it in URLs and in a README
   note explaining the repo/package name difference; remove it everywhere else.

Also restore a `lag.sparse` shim: `lag_sparse` was renamed with no deprecation path, silently breaking
any existing caller.

*Verify:* `R CMD INSTALL .`; `library(multiScaleAC)`; every Phase 2 fixture reproduces unchanged.
*Commit:* `rename package to multiScaleAC to match the companion manuscript`

## Phase 4 — Compute core rewrite

Replace `src/moran_multi_cpp.cpp` with `src/moran_kernel.cpp`. One path, two streaming passes over
tiles of the pair matrix, with the pair distance computed **once** and reused across all bandwidths:

- **Pass 1** accumulates `rowsums[h]` for every bandwidth simultaneously.
- **Pass 2** accumulates `v_h = K(h) %*% (x / rowsums[h])` for every bandwidth.

Working memory becomes **O(N · n_h)** instead of O(N²) — ~1.1 MB rather than ~1.6 GB at N = 7122, and
tractable at N = 100,000 where the current code needs ~80 GB. The distance is evaluated once per pair
rather than once per pair per bandwidth. Both passes parallelize over rows with OpenMP.

Fixes folded in:

- The uninitialized-`dist_sq` read: there is no separate `dist_sq`/`dist` buffer to get out of sync —
  the tile computes the requested metric inline, so `gaussian` + `manhattan` is correct by construction.
- `eps`: a single `w > eps` predicate (strict, matching today's sparse path and
  `gaussian_kernel_weights`) applied **before** the row sum accumulates, so thresholded weights are
  excluded from the denominator in every configuration. At `eps = 0` this reduces exactly to ground truth.
- **NaN policy**, unified and made explicit via a new `nan_policy = c("zero", "drop", "na")` argument,
  default `"zero"`. Note the subtlety: ground truth divides by a zero row sum, producing `NaN` in `v_h`,
  which propagates through the dot product, so `out$observed[is.nan(...)] = 0` zeroes the **entire
  bandwidth**, not just the offending cell. `"zero"` reproduces that exactly (required — it is the
  oracle) and warns; `"drop"` excludes isolated cells from both numerator and denominator; `"na"`
  returns `NA` so it cannot pass unnoticed. Documented in one place.
- Dead weight removed: the unread `large_mode_n` and `seed` C++ arguments, `colsums`, the shadowed
  `SparseMatrixRM W`, and the four `perms == 0` guards for a state `.validate_inputs` forbids — except
  that `perms = 0` becomes **legal** in Phase 8, so the guard is kept deliberately and tested.
- `src/Makevars` currently sets only `CXX11*` with no `CXX_STD`, so it is inert and `R CMD check`
  flags it. Replace with `PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)` / `PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS)`
  and add `Makevars.win`.

**Floating-point tolerance, decided rather than assumed.** Ground truth computes
`as.matrix(dist(coords))^2` — a square root followed by a square — whereas the tile computes
`dx*dx + dy*dy`. These differ by ~1 ulp before `exp()` amplifies them, so bit-identity is not available
and the fixtures must be compared at a tolerance. Eigen's vectorized `exp` adds a second source of
divergence. **First action after the build succeeds is to measure both** (tile row sums vs
`rowSums(K_ref)`) and record the observed maximum in the test file, asserting at `1e-12` relative. If
the measurement exceeds ~1e-13, the fix is scalar `std::exp` in the kernel (bit-comparable, 4–8×
slower) — not a loosened tolerance.

*Verify:* every Phase 2 fixture matches to 1e-12 on both `observed` and `permutations`; peak RSS
measured at N = 7122 and at N = 20000 (previously impossible); a timing comparison recorded.
*Commit:* `rewrite compute core with tiled two-pass kernel and constant working memory`

## Phase 5 — R wiring and confirmed bug fixes

- `multi_scale_ac.R:502-503` — `deparse(substitute(gene1))` resolves in the wrong frame, so gene names
  are *always* the literals `"gene1"`/`"gene2"`, polluting every plot title and every mxfda
  `interaction_name`. Capture from `match.call()` in the caller's frame instead.
- `multi_scale_ac.R:403` — the custom path errors whenever `kernel` is built-in, so a custom
  `distance` with the default kernel always fails; the documented example at `:106-117` throws. Route
  built-in kernels with custom distances through the shared assembler.
- `:442` `sprintf("%d%%", pct)` on a double; `:323` `.sanitize_h` crashing on `NA` (`if (NA > 0)`);
  `h` never checked for `NA`/`Inf`/duplicates.
- Fold the never-called exported `validate_inputs()` (`helpers.R:26`, a divergent duplicate of the
  internal `.validate_inputs`) into one validator, fixing its ordering bugs (`:34` uses `nrow()` before
  the matrix check; `:57`/`:73` run before the `NA` check at `:82`, so `NA` input dies with an obscure
  error instead of the friendly message it was written to give).
- `suggest_bandwidths()` and `qc_neighborhoods()` each materialize a full N×N `dist()`, so both are
  unusable at the N this package targets. Rewrite over tiles. `suggest_bandwidths()` currently returns
  quantiles of the *same* 1-NN distribution for all three `range` presets, so `"global"` is not global —
  base `"local"` on nearest-neighbour distances and `"global"` on the tissue extent, and return values
  in the data's own units (the mm/µm trap above). `qc_neighborhoods()` reports isolation using a hard
  `dist <= h` cutoff — the *uniform* kernel — while the default kernel is Gaussian with infinite
  support, so under the default analysis it can never report an isolated cell; report effective row
  sums for the kernel actually in use, and accept a bandwidth vector.
- Extract a shared `.assemble()` that owns the NaN policy and result construction. It must keep the
  `V` matrix (`N × n_h`) alive until the permutation step completes, so Phase 8's selective-permutation
  mode does not have to recompute the kernel.
- Add `setMethod("show", ...)`. There is none, so bare auto-printing at the REPL dumps the whole
  permutation and coordinate matrices; the S3 `print()` only fires on an explicit `print(x)` call.

*Verify:* fixtures still match; gene names appear correctly in `print`/`plot`/`export_mxfda`.
*Commit:* `fix input validation, gene name capture, and helper scalability`

## Phase 6 — Test suite

`tests/testthat/`, with `Suggests: testthat (>= 3.0.0)`. Beyond fixture agreement:

- **A test per fixed bug**, each written to fail against the Phase 1 baseline. The one that matters
  most: `kernel="gaussian", distance="manhattan"` compared against an explicit R Manhattan
  implementation — the uninitialized-memory bug produces plausible-looking non-`NA` values, so the
  existing `!is.na()` style of check reports PASSED on garbage.
- A guard asserting `R/soupir_moran_multi.R` is unmodified (hash the `deparse()` of the function), so
  the oracle cannot be silently edited.
- p-value correctness: `"greater"`/`"less"`/`"two.sided"` against hand-computed values; `p > 0` always.
- Bandwidth sanity: a test that fails if an example or fixture uses an `h` where the mean row sum
  exceeds ~80% of N — the degeneracy that made the old validation vacuous.
- Edge cases: single bandwidth, `perms = 1`, isolated cells under each `nan_policy`, zero-variance gene,
  duplicate coordinates, `N < 5`.
- Because Phase 8 may change the default inference method, tests must read defaults via `formals()`
  rather than hard-coding them.

*Verify:* `devtools::test()` green; `R CMD check --no-manual` clean.
*Commit:* `add test suite pinned to the reference implementation`

## Phase 7 — Significance rewrite

- **Honest tails.** `alternative` default becomes `"greater"` — the paper's one-tailed upper test, which
  is what the code already computed while claiming `"two.sided"`. Implement a real two-sided
  (`min(p_lo, p_hi)` doubled, or `|I| ≥ |I_obs|`) and a correct `"less"`
  (`(sum(perms <= obs) + 1)/(n + 1)`), replacing the current `1 - p`.
- **Estimator made explicit.** `estimator = c("add_one", "paper")`: `(r+1)/(n+1)` (default; the
  defensible choice, and it cannot return 0) versus the manuscript's `r/n`. Both available, the
  difference documented rather than argued.
- **Envelope, not CI.** `ci_lower`/`ci_upper` are quantiles of the *null* distribution — an acceptance
  envelope — but are named and plotted as a confidence interval around `observed`. Rename to
  `null_lower`/`null_upper` (keeping the old slots as deprecated aliases), and separate the envelope
  level from the test threshold, which currently share one `alpha`.
- Clear stale `p_adjusted` on recompute (`significance.R:81-87` leaves it, so `summary()` can show
  adjusted p-values that do not correspond to the current `p_values`).
- Note in the docs that `adjust_pvalues()` across bandwidths of one gene pair corrects across strongly
  positively dependent tests (nested kernels on the same data); the manuscript applies FDR across
  *gene pairs*, which is what `multiScaleAC_batch()` output supports.
- **Add the whole-curve supremum test** that the paper defines and the repo never implemented: for each
  permutation, the maximum absolute deviation from the mean curve, compared against the observed
  maximum. This is a pure O(n_h × perms) reduction of the matrix already in hand — no analytical
  machinery — so it lands here, cheaply, ahead of Phase 8. New slots `sup_statistic`, `sup_p_value`.

*Verify:* new tests; one-tailed p-values reproduce the Phase 1 values (confirming the old default was
one-sided); supremum p checked against a direct loop.
*Commit:* `correct one and two-sided p-values and add whole-curve supremum test`

## Phase 8 — Analytical moments (the permutation-reduction work)

The core result, verified numerically before planning. With `a = v_h - mean(v_h)`, `b = y`,
`A_k = sum(a^k)`, `B_k = sum(b^k)`, `S = sum(x^2)`, and `U = S · I_perm`:

```
E[U]   = 0
E[U^2] = A2·B2 / (n-1)
E[U^3] = n·A3·B3 / ((n-1)(n-2))
E[U^4] = A4·B4·(1/n + 4/(n(n-1))) + 3(A2²-A4)(B2²-B4)/(n(n-1))
       + 6(A2²-2A4)(B2²-2B4)/(n(n-1)(n-2)) + 9(A2²-2A4)(B2²-2B4)/(n(n-1)(n-2)(n-3))
Cov(I(h), I(h')) = <a_h, a_h'>·B2 / ((n-1)·S²)
```

`v_h` is already computed by the existing algorithm, `B2..B4` are bandwidth-independent, so the moments
cost **O(n) per bandwidth** — replacing the O(n_h · n · perms) permutation product and eliminating the
`N × perms` matrix entirely (57 MB at n = 7122, perms = 999). `eps = 0` is preserved exactly; no
thresholding or gene filtering is introduced. `mean(v_h)` is computed empirically rather than assumed
zero, so isolated rows (which break `W1 = 1`) stay correct.

**The moments→p-value step is an open question and will be settled by measurement, not assertion.**
Measured ratios of approximate to empirical upper-tail p in an adversarial skewed case (n = 180, 14
non-zero cells, skew 2.37, excess kurtosis 7.25):

| z | normal | 3-moment gamma | Edgeworth |
|---|---|---|---|
| 1.645 | 0.72× | **1.01×** | 0.60× |
| 2.326 | 0.24× | 0.89× | 1.29× |
| 3.090 | **0.04×** | 0.73× | 1.62× |
| 3.720 | **0.01×** | 0.64× | 0.89× |

So a normal/z-score p-value is catastrophically anti-conservative under skew — which matters because
the existing `z_scores` slot invites exactly that mistake — and Edgeworth is erratic and non-monotone
and can leave [0,1]. A shifted gamma is much better but not tail-exact. A double saddlepoint
(Robinson 1980; Booth & Butler 1990 — condition `sum(a_i Y_i)` on `sum(Y_i) = 0` with `Y_i` drawn from
the empirical distribution of `b`) reached 1.27× at z = 1.645 in a prototype before overflowing; it
needs log-sum-exp stabilization and becomes O(n·m) per Newton step if `b` is binned to `m` atoms.

Deliverables:

1. `R/moments.R` — exact moments and the cross-bandwidth covariance; C++ helper for the power sums.
2. `R/approx.R` — `"saddlepoint"`, `"gamma"`, `"edgeworth"`, `"normal"`, each monotone and clamped to
   (0, 1), with the saddlepoint stabilized.
3. **Bake-off harness** (`inst/validation/`, `.Rbuildignore`d output): each method against 10⁵–10⁶
   permutation ground truth, across bandwidths spanning the degenerate-to-saturated range, on
   `ccRCCExample` and on deliberately zero-inflated synthetic genes. Report tail ratios at
   p = 10⁻¹…10⁻⁵, since FDR over ~1100 gene pairs lives at 10⁻⁴–10⁻⁵. The winner becomes the default.
4. **Trustworthiness diagnostic + automatic fallback.** Exact excess kurtosis and max leverage
   `max(a_i²)/A2` per bandwidth; when either exceeds a calibrated threshold the bandwidth falls back to
   permutations. This matters most at small `h` — exactly where the manuscript says the curve is
   unstable. `permute = c("none", "flagged", "all")` reuses the retained `V` so a fallback never
   recomputes the kernel. This is the design that actually delivers the goal: permutations become a
   targeted safety net rather than the bulk cost.
5. **Simulation-based supremum** using the closed-form covariance — draw curves from the n_h × n_h
   covariance with moment-matched marginals, O(M · n_h²), never touching `n`. Held to a higher
   evidentiary bar than the marginal p-values, since a Gaussian copula is not the true joint
   permutation law; ships only if the bake-off supports it, and is reported alongside the Phase 7
   permutation supremum rather than replacing it.

API: `multiScaleAC(..., perms = 1000)` gains `moments = TRUE` (cheap, always on) and permits
`perms = 0`; `compute_significance(result, method = c("permutation", "moment", "auto"), moment_dist = ...)`.
New slots `moment_mean/var/skew/kurt`, `moment_cov`, `p_method` (per-bandwidth provenance — which method
produced each p), `leverage`. The all-or-nothing significance-slot validity rule at `classes.R:77-86`
must be revised to per-group consistency (permutation-derived and moment-derived groups checked
independently), or a moments-only result cannot be valid; `permutations` must be allowed zero columns;
and the tautological second element of the `dim()` check at `:67` fixed.

*Verify:* moments match large-permutation empirics to <1%; bake-off table committed with the default it
justified; a permutation-free run on `ccRCCExample` compared end-to-end against a 10⁵-permutation run
with peak memory and wall-clock for both.
*Commit:* `add exact analytical permutation moments with approximation bake-off and fallback`

## Phase 9 — Batch, plotting, mxfda

- `multiScaleAC_batch()`: `checkpoint[[i]] <<- result` (`:690`) is discarded inside `mclapply` forks and
  every child races on the same `saveRDS` path — replace with per-worker files reduced by the parent.
  `for (i in start_idx:n_samples)` (`:751`, and `:712`) counts **down** when a checkpoint already covers
  every sample, silently reprocessing everything in reverse — use `seq_len`/`seq_along`. Validate
  columns before subsetting (`:654` subsets first, making three carefully-written error messages
  unreachable). Guard `parallel` for Windows.
- `plot.R:337-338`: the significance-asterisk layer passes `data = sig_df` with `inherit.aes = TRUE`
  while the plot maps `fill = value`, so `plot_batch_comparison()` errors on its **default**
  `value_var = "z_score"` as soon as anything is significant. Also: `show_permutations` and `ci_level`
  are documented but never used; `cluster_rows` does no clustering; `sort_by = "sample_id"` is
  documented but unimplemented; `alpha` in `plot_significance_summary` only reaches the title. Implement
  or remove each — no argument may be documented and inert.
- `export_mxfda()`: `export_mxfda_batch()` `rbind`-errors when some results carry `p_adjusted` and
  others do not (a likely state after `lapply(results, compute_significance)`) — bind on a fixed column
  set. Add `mxfda` to `Suggests` and a **skip-if-not-installed** test that actually calls
  `mxfda::add_summary_function()`, since the "mxfda-compatible" claim has never been executed. Document
  explicitly that bandwidth `h` occupies mxfda's radius column `r`: they are different objects (a
  Gaussian scale versus a hard radius) and the substitution needs stating, not assuming.

*Verify:* batch run over `ccRCCExample` split by FOV, sequential and parallel, results identical;
checkpoint interruption and resume exercised; every plot function called on a real result.
*Commit:* `fix batch checkpointing, plot layers, and mxfda export column handling`

## Phase 10 — Documentation and user-friendliness

- **README.md** is currently 21 bytes. Write: what the method does, install, a 20-line quickstart on
  `ccRCCExample` at unit-correct bandwidths, the citation, and a short note that the package is
  `multiScaleAC` while the repository remains `KernelCommunication` (with the manuscript URL).
- **Vignette** `vignettes/multiScaleAC.Rmd`: choosing bandwidths (with the mm/µm trap made explicit and
  the degeneracy demonstrated — this is the single most likely user error), interpreting the curve,
  significance and the permutation/moment trade-off, and the mxfda/FPCA hand-off the paper performs.
- **Runnable examples.** Almost every example is wrapped in `\dontrun{}`, so the flagship function has
  no executed example. Make them run on small synthetic data within CRAN time limits, at bandwidths from
  the `H_MM` constant.
- `NEWS.md` recording the rename, the corrected p-values, the memory reduction and the moments feature,
  including an explicit note that results from the affected `gaussian`+`manhattan` configuration and
  from `alternative = "less"` were wrong in the prior state.
- pkgdown: **yes**, but `docs/` already holds hand-written notes, so set `destination: reference` in
  `_pkgdown.yml` or first move the notes to `inst/notes/` — do not let pkgdown overwrite them.
- Deprecate gently: `soupir_moran_multi()` and `soupir_moran_bv()` stay exported and undeprecated (they
  are the oracle and the backwards-compatibility contract). `multiScaleAC()` docs currently call
  `soupir_moran_multi()` "deprecated" — remove that, since it is not.

*Verify:* `R CMD check --as-cran`; vignette builds; every example runs; pkgdown site builds.
*Commit:* `add README, vignette, runnable examples, and release notes`

## Phase 11 — Kernels and distances (lowest priority)

Explicitly secondary. Minimal correct version: `gaussian`, `uniform`, `exponential` kernels and
`euclidean`, `manhattan`, `minkowski(p)` distances, all inside the single tiled core so none can bypass
the `eps` predicate or the NaN policy. The current custom-kernel path (`multi_scale_ac.R:392-418`) never
zeroes the diagonal, so custom results are on a different scale from built-in ones and not comparable —
fix or document loudly. Keep the R escape hatch for a user-supplied function, but state its cost and
that it is exact rather than tiled. Kernels must be normalizable (row sums > 0) for the moment
identities to hold; validate and refuse otherwise.

*Verify:* each built-in kernel/distance combination against an explicit R implementation.
*Commit:* `support additional kernels and distance measures in the tiled core`

## Phase 12 — Final gate

`R CMD check --as-cran` clean (0 errors, 0 warnings, notes explained); full test suite; a benchmark
recording time and peak memory before/after at N = 7122 and N = 20000, permutation versus moment-based;
`grep -ri` sweep confirming no forbidden attribution anywhere including commit messages
(`git log --format=%B | grep -i`); confirm `git remote` was never pushed to and the working tree is clean.

*Commit:* `finalize checks, benchmarks, and documentation for release`

---

## Verification strategy, end to end

The oracle is `soupir_moran_multi()` — pure base R, never modified, hash-guarded by a test. Every phase
runs the Phase 2 fixtures and must reproduce `observed` and `permutations` to 1e-12 relative; the one
tolerance decision (sqrt-then-square versus `dx*dx+dy*dy`, plus vectorized `exp`) is **measured in
Phase 4 and the observed value recorded in the test**, not assumed. `R CMD check --no-manual` runs every
phase so check findings do not pile up; `--as-cran` gates Phase 10 and 12.

The analytical-moments claim is verified twice over: the moments against large-permutation empirics
(<1%, already demonstrated), and the resulting p-values against 10⁵–10⁶ permutation ground truth in the
far tail where FDR screening actually operates. If a method cannot be shown calibrated at p ≈ 10⁻⁴, it
does not become the default — permutations remain the default and the moments ship as an opt-in
accelerator with a diagnostic. I will report what the bake-off actually shows, including if it shows
that no approximation is trustworthy enough to replace permutations.

## Risks

- **The bake-off may not produce a usable default.** That is a real possible outcome. The exact moments
  are still valuable on their own (a free, exact variance and skewness at every bandwidth, plus the
  supremum covariance), and the fallback design degrades to today's behaviour.
- **OpenMP on macOS** is not available under stock Apple clang. The `$(SHLIB_OPENMP_CXXFLAGS)` form
  degrades to single-threaded rather than failing to build, which is the correct trade.
- **Correcting `alternative = "less"` and the manhattan bug changes published-adjacent numbers.** Any
  prior result from those configurations was wrong; `NEWS.md` will say so plainly.
- **Phase 8 depends on Phase 5 retaining `V`.** Noted as a comment in `.assemble()` so it is not lost.
