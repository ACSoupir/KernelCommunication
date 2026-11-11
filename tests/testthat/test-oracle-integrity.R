# soupir_moran_multi() is the reference implementation for this package: every golden fixture is
# generated from it, and every claim of correctness is relative to it. It must not be modified.
#
# This test exists because the oracle is only useful if it is stable. If it is edited -- even
# "harmlessly" -- the fixtures silently start describing a different function, and the entire test
# suite continues to pass while measuring nothing.
#
# If this test fails, do NOT update the checksum to make it pass. Establish why the function changed.
# The only legitimate reason to update the constant below is a deliberate, reviewed decision to move
# the reference, which also requires regenerating every fixture.

test_that("the reference implementation is unmodified", {
  fn <- multiScaleAC::soupir_moran_multi

  # Hash the deparsed function rather than the file, so comments and roxygen edits are permitted
  # while any change to the code itself is caught.
  tf <- tempfile()
  on.exit(unlink(tf), add = TRUE)
  writeLines(deparse(fn), tf)

  expect_identical(
    unname(tools::md5sum(tf)),
    "6b05a6a77206c95a57e14cd2c79c55e0",
    info = paste("R/soupir_moran_multi.R appears to have been modified.",
                 "It is the regression oracle and must not change.")
  )
})

test_that("the reference implementation keeps its signature", {
  # Argument order matters: coords comes first, unlike soupir_moran_bv(). Callers and fixtures
  # depend on this, and an earlier version of the documented example had the order wrong.
  expect_identical(
    names(formals(multiScaleAC::soupir_moran_multi)),
    c("coords", "gene1", "gene2", "h_vec", "perms", "seed")
  )
  expect_identical(formals(multiScaleAC::soupir_moran_multi)$seed, 333)
})

test_that("the reference implementation stays exported for backwards compatibility", {
  # soupir_moran_multi() is a permanent part of the public API. It is not deprecated.
  exports <- getNamespaceExports("multiScaleAC")
  expect_true("soupir_moran_multi" %in% exports)
  expect_true("soupir_moran_bv" %in% exports)
})

test_that("the reference implementation still reproduces its own fixtures", {
  # Guards against an environment change (RNG stream, sample() algorithm, BLAS) silently altering
  # what the oracle produces. If this fails, the fixtures are stale relative to this machine, not
  # necessarily wrong.
  for (cs in golden_cases()) {
    res <- multiScaleAC::soupir_moran_multi(
      cs$coords, cs$gene1, cs$gene2,
      h_vec = cs$h, perms = cs$perms, seed = cs$seed
    )
    expect_equal(res$observed, cs$observed, tolerance = 1e-14,
                 info = paste("case:", cs$name))
    expect_equal(res$permutations, cs$permutations, tolerance = 1e-14,
                 info = paste("case:", cs$name))
  }
})
