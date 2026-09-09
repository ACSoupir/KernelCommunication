> # ⛔ RETRACTED — THIS DOCUMENT IS WRONG. DO NOT ACT ON IT.
>
> The conclusion below is false. `soupir_moran_multi()` is **correct** and is the reference
> implementation for this package. Nothing in it needs fixing, and no published analysis needs
> re-running on account of it.
>
> Two specific errors:
>
> 1. **The algebra.** Lines 76 and 82 below write the "legacy" and "correct" sums as
>    `Σ_i Σ_j g2[j]·K[i,j]·g1[i]/r[i]` and `Σ_i Σ_j g1[i]·K[i,j]·g2[j]/r[i]`. These are the *same
>    expression* — an identical product written in a different factor order. The document then
>    relabels the bound dummy indices on one side only and compares the relabelled form to the
>    original, concluding non-equivalence. Relabelling bound indices cannot change a sum's value.
> 2. **The 73% figure.** The demonstration uses `W_correct <- t(t(K) / r_sums)`, which evaluates to
>    `K[i,j] / r_sums[j]` — **column**-normalized, not row-normalized. The 73% is the gap between
>    row- and column-normalization, not an error in the reference. That snippet also uses `g1_s`/`g2_s`
>    without ever defining them, so it cannot have run as written.
>
> Verified numerically: `v_h = K %*% (x/r)` equals `t(W) %*% x` to `2.2e-16`, so `I = x'Wy/S` exactly.
> Also note the "Verification" section's claim that `E[I] = -1/(N-1)`: that is the *univariate* result.
> Under this package's permutation scheme (permuting `y` only, both genes centred), `E[I] = 0` exactly.
>
> Full analysis: [`docs/RETRACTIONS.md`](../RETRACTIONS.md) §1. Retained only as a record of the claim.

---

# Issue #1: Legacy Function Row-Normalization Bug

**Severity:** ~~Critical (affects mathematical correctness)~~ — retracted, see banner above

**Date:** 2026-06-26

**Status:** RETRACTED — the reference implementation was never buggy

---

## Description

The legacy `soupir_moran_multi()` function uses an incorrect row-normalization formula that produces mathematically invalid Moran's I values.

---

## Root Cause Analysis

### Legacy Implementation (BUGGY)

```r
# From soupir_moran_multi.R lines 49-54
x_weighted <- g1_s / r_sums
v_h <- as.vector(K %*% x_weighted)
perm_vals <- (v_h %*% Y_perm) / denom
obs_val   <- (v_h %*% g2_s) / denom
```

This computes:
```
I_legacy = ((K × (g1_s ⊘ r_sums))' × g2_s) / denom
```

Where `⊘` denotes element-wise division.

### New Implementation (CORRECT)

```cpp
// From moran_multi_cpp.cpp lines 137-141
// Row-normalize: W(i,j) = W(i,j) / rowsum(i)
for (int i = 0; i < N; i++) {
  if (rowsums(i) > 0) {
    for (SparseMatrix<double>::InnerIterator it(W, i); it; ++it) {
      it.valueRef() /= rowsums(i);
    }
  }
}

// Compute observed Moran's I
Eigen::VectorXd U = W * gene2;
observed(h_idx) = (gene1.dot(U)) / gene1.squaredNorm();
```

This computes:
```
I_new = (g1_s' × (row_norm(K)) × g2_s) / denom
```

Where `row_norm(K)[i,j] = K[i,j] / r_sums[i]`.

---

## Mathematical Proof of Difference

The standard bivariate Moran's I formula is:

```
I = (g1' × W × g2) / (g1' × g1)
```

Where W should be **row-normalized**: `W[i,j] = K[i,j] / Σ_k K[i,k]`

### Legacy approach expands to:
```
I_legacy = Σ_j g2_s[j] × Σ_i K[i,j] × (g1_s[i] / r_sums[i])
         = Σ_j Σ_i g2_s[j] × K[i,j] × g1_s[i] / r_sums[i]
```

### Correct approach expands to:
```
I_correct = Σ_i g1_s[i] × Σ_j (K[i,j] / r_sums[i]) × g2_s[j]
          = Σ_i Σ_j g1_s[i] × K[i,j] × g2_s[j] / r_sums[i]
```

### Are they equivalent?

For them to be equal, we need:
```
Σ_j Σ_i g2_s[j] × K[i,j] × g1_s[i] / r_sums[i] = Σ_i Σ_j g1_s[i] × K[i,j] × g2_s[j] / r_sums[i]
```

Since K is symmetric (K[i,j] = K[j,i]), we can swap indices on the left:
```
Left = Σ_i Σ_j g2_s[i] × K[j,i] × g1_s[j] / r_sums[j]
     = Σ_i Σ_j g2_s[i] × K[i,j] × g1_s[j] / r_sums[j]
```

For equality we need:
```
Σ_i Σ_j g2_s[i] × K[i,j] × g1_s[j] / r_sums[j] = Σ_i Σ_j g1_s[i] × K[i,j] × g2_s[j] / r_sums[i]
```

These are **NOT equivalent** unless `r_sums[i] = r_sums[j]` for all i,j (which is generally false).

---

## Empirical Demonstration

```r
set.seed(1)
N <- 100
coords <- matrix(runif(2*N, -50, 50), ncol = 2)
g1 <- rnorm(N)
g2 <- rnorm(N)
h <- 10

# Legacy (buggy)
dist_sq <- as.matrix(dist(coords))^2
K <- exp(-dist_sq / (2 * h^2))
diag(K) <- 0
r_sums <- rowSums(K)
x_weighted <- g1_s / r_sums
v_h <- as.vector(K %*% x_weighted)
I_legacy <- sum(v_h * g2_s) / sum(g1_s^2)
# Result: 0.01800614

# Correct row-normalization
W_correct <- t(t(K) / r_sums)
v_h_correct <- as.vector(W_correct %*% g2_s)
I_correct <- sum(g1_s * v_h_correct) / sum(g1_s^2)
# Result: 0.01037501

# Difference: 0.007631126 (73% relative error!)
```

---

## Impact

The legacy function produces **incorrect Moran's I values** that do not correspond to the standard mathematical definition. The error magnitude depends on the heterogeneity of row sums, which varies with:
- Spatial point distribution
- Bandwidth parameter h
- Kernel function

For the test case above, the relative error is **73%**, which is substantial.

---

## Resolution

**Decision:** The new `multiScaleAC()` function implements the **correct** formula and should NOT be modified to match the buggy legacy behavior.

**Actions:**
1. ✅ Keep new implementation as-is (mathematically correct)
2. ✅ Document the legacy bug in validation report
3. ✅ Update legacy function with deprecation warning noting the bug fix
4. ✅ Add note to migration guide explaining the difference

---

## Verification

To verify the new implementation is correct:

1. **Compare to standard spatial econometrics packages:**
   - Compare to `spdep::moran()` with appropriate weight matrix
   - Compare to `spatialreg` package functions

2. **Test with known spatial patterns:**
   - Perfect positive autocorrelation should give I ≈ 1
   - Perfect negative autocorrelation should give I ≈ -1
   - Random pattern should give I ≈ E[I] = -1/(N-1)

3. **Check expectation:**
   - Under null hypothesis of no spatial autocorrelation, E[I] should equal -1/(N-1)
   - Verify this holds for the new implementation

---

## Recommendation for Users

Users migrating from `soupir_moran_multi()` to `multiScaleAC()` should:
1. **Expect different absolute values** of Moran's I
2. **Trust the new values** as mathematically correct
3. **Re-run analyses** with the new function
4. **Note in publications** that corrected formula was used

---

## References

1. Moran, P. A. P. (1950). Notes on continuous stochastic phenomena. *Biometrika*, 37(1/2), 17-23.
2. Cliff, A. D., & Ord, J. K. (1981). *Spatial processes: Models & applications*. Pion.
3. Bivand, R. S., Pebesma, E., & Gómez-Rubio, V. (2013). *Applied spatial data analysis with R*. Springer.
