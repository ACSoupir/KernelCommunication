# Issue Log: multiScaleAC() Validation

## Issue #001: Legacy Function Missing Import

**Date:** 2026-06-26  
**Severity:** Low  
**Status:** RESOLVED  

**Description:** `soupir_moran_multi()` used `%>%` pipe operator without importing magrittr.

**Resolution:** Replaced pipe operations with base R equivalents.

**Files Modified:** `R/soupir_moran_multi.R`

---

## Issue #002: C++ Coordinate Centering Numerical Error

**Date:** 2026-06-26  
**Severity:** Medium  
**Status:** RESOLVED  

**Description:** C++ code centered coordinates before computing distances, introducing floating-point errors that caused distance mismatches with R's `dist()` function.

**Root Cause:**
```cpp
Eigen::VectorXd coords_mean = coords.colwise().mean();
coords_centered.rowwise() -= coords_mean.transpose();
```

**Resolution:** Removed coordinate centering since Euclidean distances are translation-invariant.

**Files Modified:** `src/moran_multi_cpp.cpp`

---

## Issue #003: C++ Sparse Matrix Row-Normalization Bug

**Date:** 2026-06-26  
**Severity:** Critical  
**Status:** RESOLVED  

**Description:** C++ implementation used column-major sparse matrix with InnerIterator for row-normalization, causing incorrect normalization.

**Root Cause:**
```cpp
// Column-major storage, InnerIterator iterates columns
Eigen::SparseMatrix<double> W(N, N);
for (int i = 0; i < N; i++) {
  for (SparseMatrix<double>::InnerIterator it(W, i); it; ++it) {
    it.valueRef() /= rowsums(i);  // Wrong iteration order
  }
}
```

**Resolution:** Changed to RowMajor storage so InnerIterator correctly iterates over columns within each row:
```cpp
using SparseMatrixRM = Eigen::SparseMatrix<double, Eigen::RowMajor>;
SparseMatrixRM W(N, N);
for (int i = 0; i < N; i++) {
  for (SparseMatrixRM::InnerIterator it(W, i); it; ++it) {
    it.valueRef() /= rowsums(i);  // Correct
  }
}
```

**Files Modified:** `src/moran_multi_cpp.cpp`

---

## Issue #004: Permutation Value Mismatch — RETRACTED

This entry claimed that permutation values differ between the R and C++ paths because the C++ code
used `std::rand()` with a Fisher-Yates shuffle, and recorded that as expected behaviour.

The C++ code contains no random number generator. The permutation matrix is built in R with the same
`set.seed` / `replicate` / `sample` sequence as the reference implementation, and the same matrix is
reused across all bandwidths in both, so the permutations agree exactly. See `docs/RETRACTIONS.md` §3.

---

## Summary

| Issue | Severity | Status |
|-------|----------|--------|
| #001: Missing import | Low | RESOLVED |
| #002: Coordinate centering | Medium | RESOLVED |
| #003: Row-normalization | Critical | RESOLVED |
| #004: Permutation mismatch | — | RETRACTED (see `docs/RETRACTIONS.md`) |

Issues #001-#003 are sound. For claims made elsewhere in this repository that did not hold, see
`docs/RETRACTIONS.md`.
