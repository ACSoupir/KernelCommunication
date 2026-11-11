# The core invariant of this package: multiScaleAC() must reproduce soupir_moran_multi().
#
# This must stay green through every refactor. It is the test that gives the compute-core rewrite,
# the package rename, and the analytical-moment work their licence to proceed.
#
# Cases with a zero row sum are handled separately in test-nan-policy.R, because they exercise a
# semantics question rather than plain numerical agreement.

WELL_BEHAVED <- c("synth_uniform", "synth_correlated", "synth_zeroinflated", "ccrcc_subsample")

test_that("multiScaleAC reproduces the reference observed curve", {
  for (nm in WELL_BEHAVED) {
    cs <- golden_cases()[[nm]]
    res <- multiScaleAC(cs$coords, cs$gene1, cs$gene2,
                        h = cs$h, perms = cs$perms, seed = cs$seed, verbose = FALSE)
    expect_equal(res@observed, cs$observed,
                 tolerance = REFERENCE_TOLERANCE, info = paste("case:", nm))
  }
})

test_that("multiScaleAC reproduces the reference permutation matrix", {
  # Permutations must agree exactly, not just distributionally: both paths draw the same permutation
  # set from the same seed and reuse it across all bandwidths. An earlier note in this repository
  # claimed they differ because of a C++ RNG; there is no RNG in the C++ at all.
  for (nm in WELL_BEHAVED) {
    cs <- golden_cases()[[nm]]
    res <- multiScaleAC(cs$coords, cs$gene1, cs$gene2,
                        h = cs$h, perms = cs$perms, seed = cs$seed, verbose = FALSE)
    expect_equal(dim(res@permutations), dim(cs$permutations), info = paste("case:", nm))
    # ignore_attr: the reference builds its matrix via t(sapply(...)), which carries "h_<value>"
    # rownames; the compute core returns a bare matrix. Only the values are being compared here.
    # Row naming is asserted separately below.
    expect_equal(res@permutations, cs$permutations,
                 tolerance = REFERENCE_TOLERANCE, ignore_attr = TRUE,
                 info = paste("case:", nm))
  }
})

test_that("agreement is far tighter than the asserted tolerance", {
  # Bit-identity is unavailable (the reference squares a square root; the core forms dx*dx + dy*dy),
  # so this records the margin actually available rather than leaving the tolerance unexamined. If
  # this starts failing, the compute core has drifted even though the assertions above still pass.
  worst <- 0
  for (nm in WELL_BEHAVED) {
    cs <- golden_cases()[[nm]]
    res <- multiScaleAC(cs$coords, cs$gene1, cs$gene2,
                        h = cs$h, perms = cs$perms, seed = cs$seed, verbose = FALSE)
    scale_obs <- max(abs(cs$observed), .Machine$double.eps)
    worst <- max(worst,
                 max(abs(res@observed - cs$observed)) / scale_obs,
                 max(abs(res@permutations - cs$permutations)) / max(abs(cs$permutations)))
  }
  # Measured on R 4.6.1 / aarch64-apple-darwin: ~1e-15 relative. Asserted at 1e-13 to leave room for
  # a different BLAS or a vectorised exp() while still being 100x below REFERENCE_TOLERANCE.
  expect_lt(worst, 1e-13)
})

test_that("a single bandwidth gives the same answer as that element of a vector", {
  # Guards against any bandwidth-loop state leaking between iterations.
  cs <- golden_cases()$synth_correlated
  full <- multiScaleAC(cs$coords, cs$gene1, cs$gene2,
                       h = cs$h, perms = cs$perms, seed = cs$seed, verbose = FALSE)
  for (i in seq_along(cs$h)) {
    one <- multiScaleAC(cs$coords, cs$gene1, cs$gene2,
                        h = cs$h[i], perms = cs$perms, seed = cs$seed, verbose = FALSE)
    expect_equal(one@observed, full@observed[i],
                 tolerance = REFERENCE_TOLERANCE, info = paste("h =", cs$h[i]))
  }
})

test_that("bandwidth order does not affect the result", {
  cs <- golden_cases()$synth_uniform
  fwd <- multiScaleAC(cs$coords, cs$gene1, cs$gene2,
                      h = cs$h, perms = 19, seed = 7, verbose = FALSE)
  rev <- multiScaleAC(cs$coords, cs$gene1, cs$gene2,
                      h = rev(cs$h), perms = 19, seed = 7, verbose = FALSE)
  expect_equal(rev@observed, base::rev(fwd@observed), tolerance = REFERENCE_TOLERANCE)
})

test_that("results are deterministic for a fixed seed", {
  # Not a formality: reading uninitialised memory produces plausible values that differ between
  # runs, which is exactly how the gaussian+manhattan defect stayed hidden. Any non-determinism
  # here is a bug.
  cs <- golden_cases()$synth_uniform
  a <- multiScaleAC(cs$coords, cs$gene1, cs$gene2, h = cs$h, perms = 29, seed = 11, verbose = FALSE)
  b <- multiScaleAC(cs$coords, cs$gene1, cs$gene2, h = cs$h, perms = 29, seed = 11, verbose = FALSE)
  expect_identical(a@observed, b@observed)
  expect_identical(a@permutations, b@permutations)
})
