# Guards against the failure mode that made this repository's earlier validation worthless.
#
# ccRCCExample coordinates are in millimetres, spanning 1.02 x 1.02 mm. Earlier benchmarks and
# examples used h = 10..80 and h = 25..250, which drives every off-diagonal Gaussian weight to
# >= 0.9999 and the mean row sum to over 7119 of 7122 cells. The kernel becomes a constant matrix,
# Moran's I collapses to ~2e-6, and two implementations agree to 1e-18 on what is essentially an
# arithmetic identity. "OVERALL: PASSED" was reported. Nothing had been tested.
#
# See docs/RETRACTIONS.md section 6.

test_that("no fixture bandwidth is in the saturated regime", {
  for (cs in golden_cases()) {
    n   <- nrow(cs$coords)
    sat <- cs$diagnostics[, "mean_rowsum"] / n
    expect_true(
      all(sat < MAX_SATURATION),
      info = sprintf(
        "case %s: bandwidths %s give mean row sums %s of N=%d (saturation %s). A saturated kernel is constant, so agreement there proves nothing.",
        cs$name,
        paste(cs$h[sat >= MAX_SATURATION], collapse = ", "),
        paste(round(cs$diagnostics[sat >= MAX_SATURATION, "mean_rowsum"], 1), collapse = ", "),
        n,
        paste(round(sat[sat >= MAX_SATURATION], 3), collapse = ", ")
      )
    )
  }
})

test_that("fixtures span a wide range of neighbourhood sizes", {
  # A suite that only ever exercises one effective neighbourhood size cannot detect scale-dependent
  # errors. Across the whole fixture set we want both the near-empty and the well-populated regime.
  all_rowsums <- unlist(lapply(golden_cases(), function(cs) cs$diagnostics[, "mean_rowsum"]))
  expect_lt(min(all_rowsums), 1)      # near-empty neighbourhoods, where instability lives
  expect_gt(max(all_rowsums), 100)    # well-populated neighbourhoods
})

test_that("at least one fixture has a strongly non-zero statistic", {
  # If every fixture sat near I = 0, a sign error or a dropped normalisation could hide. The
  # synth_correlated case is built to have real spatial cross-correlation.
  peak <- max(vapply(golden_cases(), function(cs) max(abs(cs$observed)), numeric(1)))
  expect_gt(peak, 0.1)
})

test_that("at least one fixture has a skewed permutation null", {
  # The regime where a normal approximation to the permutation null fails. Required for the
  # analytical-moment calibration to be meaningful.
  skew <- unlist(lapply(golden_cases(), function(cs) {
    apply(cs$permutations, 1, function(r) {
      s <- stats::sd(r)
      if (s == 0) return(0)
      mean((r - mean(r))^3) / s^3
    })
  }))
  expect_gt(max(abs(skew)), 1)
})

test_that("degenerate fixtures pin the zero-rowsum semantics", {
  # The reference divides gene1 by a zero row sum, producing NaN in v_h, which propagates through
  # the dot product; `observed[is.nan(observed)] <- 0` then zeroes the ENTIRE bandwidth, not just the
  # offending cell. Any reimplementation must reproduce that. These two cases capture it.
  iso <- golden_cases()$isolated_cell
  expect_true(all(iso$observed == 0))
  expect_true(all(iso$diagnostics[, "n_isolated"] >= 1))

  zv <- golden_cases()$zero_variance
  expect_true(all(zv$observed == 0))
})

test_that("fixture metadata is self-consistent", {
  g <- golden()
  expect_identical(g$generated_by,
                   "soupir_moran_multi() -- reference implementation, base R only")
  for (cs in g$cases) {
    expect_length(cs$observed, length(cs$h))
    expect_identical(nrow(cs$permutations), length(cs$h))
    expect_identical(ncol(cs$permutations), as.integer(cs$perms))
    expect_identical(nrow(cs$coords), length(cs$gene1))
    expect_identical(length(cs$gene1), length(cs$gene2))
    expect_identical(nrow(cs$diagnostics), length(cs$h))
  }
})
