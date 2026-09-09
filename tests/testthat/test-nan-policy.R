# Zero-row-sum semantics.
#
# When a cell has no neighbours carrying weight, its row sum is 0. The reference divides gene1 by
# that row sum, so v_h picks up Inf/NaN, the NaN propagates through the dot product, and
# `observed[is.nan(observed)] <- 0` then zeroes THE ENTIRE BANDWIDTH -- not merely the offending
# cell. The same applies to a zero-variance gene, where scale() yields NaN everywhere.
#
# That behaviour is surprising, and arguably not what a user wants, but it is what the reference
# does, so it is what the default must reproduce. Improving on it silently would break the oracle,
# which is the most valuable asset in this repository. Alternatives are offered via nan_policy.
#
# MEASURED DIVERGENCE at the point these tests were written (commit 6aba802): the compute core does
# not apply the NaN sweep at all, giving max|observed difference| = 6.0e-02 on isolated_cell and NaN
# on zero_variance, against a reference that returns exactly 0. Three different policies existed
# across the three code paths (none in dense C++, a guard in sparse C++, row_sums[row_sums==0] <- 1
# in the R custom path). These tests are expected to FAIL until the compute-core rewrite unifies
# them, and they are committed failing on purpose so the gap is visible rather than described.

DEGENERATE <- c("isolated_cell", "zero_variance")

test_that("zero row sums reproduce the reference, which zeroes the whole bandwidth", {
  for (nm in DEGENERATE) {
    cs <- golden_cases()[[nm]]
    res <- suppressWarnings(
      multiScaleAC(cs$coords, cs$gene1, cs$gene2,
                   h = cs$h, perms = cs$perms, seed = cs$seed, verbose = FALSE)
    )
    expect_equal(res@observed, cs$observed,
                 tolerance = REFERENCE_TOLERANCE, info = paste("case:", nm))
    expect_equal(res@permutations, cs$permutations,
                 tolerance = REFERENCE_TOLERANCE, ignore_attr = TRUE,
                 info = paste("case:", nm))
  }
})

test_that("no result ever contains a non-finite value under the default policy", {
  # The default policy must never let Inf or NaN reach the result object, because they propagate
  # silently into perm_sd, quantile() and the p-values downstream.
  for (nm in DEGENERATE) {
    cs <- golden_cases()[[nm]]
    res <- suppressWarnings(
      multiScaleAC(cs$coords, cs$gene1, cs$gene2,
                   h = cs$h, perms = cs$perms, seed = cs$seed, verbose = FALSE)
    )
    expect_true(all(is.finite(res@observed)), info = paste("observed, case:", nm))
    expect_true(all(is.finite(res@permutations)), info = paste("permutations, case:", nm))
  }
})

test_that("an isolated cell is reported to the user, not silently absorbed", {
  # Returning 0 for an entire bandwidth is a big enough surprise that it must be announced.
  cs <- golden_cases()$isolated_cell
  expect_warning(
    multiScaleAC(cs$coords, cs$gene1, cs$gene2,
                 h = cs$h, perms = 19, seed = 1, verbose = FALSE),
    regexp = "isolat|no neighbours|no neighbors|row sum",
    ignore.case = TRUE
  )
})
