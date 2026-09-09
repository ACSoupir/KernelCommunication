# Tests for the input-checking and bandwidth-selection helpers.

test_that("validate_inputs reports structural problems before dependent ones", {
  # Previously nrow(coords) was called before the matrix check had short-circuited, so a plain
  # vector produced "argument is of length zero" rather than the message written for it.
  res <- validate_inputs(1:10, rnorm(10), rnorm(10))
  expect_false(res$valid)
  expect_match(paste(res$errors, collapse = " "), "matrix")

  set.seed(1); co <- matrix(runif(20), ncol = 2)
  res <- validate_inputs(co, rnorm(5), rnorm(10))
  expect_false(res$valid)
  expect_match(paste(res$errors, collapse = " "), "length")
})

test_that("validate_inputs reports missing values rather than failing on them", {
  # mean(gene1 == 0) and sd(gene1) used to run before the NA check, so NA input died with
  # "missing value where TRUE/FALSE needed" instead of producing the intended message.
  set.seed(2); co <- matrix(runif(40), ncol = 2)
  g <- rnorm(20); g[3] <- NA
  res <- validate_inputs(co, g, rnorm(20))
  expect_false(res$valid)
  expect_match(paste(res$errors, collapse = " "), "missing")

  co2 <- co; co2[1, 1] <- NA
  expect_match(paste(validate_inputs(co2, rnorm(20), rnorm(20))$errors, collapse = " "), "missing")

  co3 <- co; co3[1, 1] <- Inf
  expect_match(paste(validate_inputs(co3, rnorm(20), rnorm(20))$errors, collapse = " "), "infinite")
})

test_that("validate_inputs flags zero variance and heavy zero inflation", {
  set.seed(3); co <- matrix(runif(200), ncol = 2)
  expect_match(paste(validate_inputs(co, rep(1, 100), rnorm(100))$errors, collapse = " "),
               "zero variance")

  sparse_gene <- c(rep(0, 95), rnorm(5))
  wrn <- validate_inputs(co, sparse_gene, rnorm(100))$warnings
  expect_match(paste(wrn, collapse = " "), "zeros")
  # The warning should say why it matters, not merely that it is true.
  expect_match(paste(wrn, collapse = " "), "skew")
})

test_that("multiScaleAC and validate_inputs agree on what is invalid", {
  # These were two separate implementations with different rules, and the exported one was never
  # called. They now share one implementation, so anything validate_inputs() calls an error (other
  # than zero variance, which is a degenerate case rather than a failure) must also stop the
  # computation.
  set.seed(4); co <- matrix(runif(40), ncol = 2)
  g <- rnorm(20); g[2] <- NA
  expect_false(validate_inputs(co, g, rnorm(20))$valid)
  expect_error(multiScaleAC(co, g, rnorm(20), h = 0.1, perms = 5, verbose = FALSE))

  expect_false(validate_inputs(co, rnorm(5), rnorm(20))$valid)
  expect_error(multiScaleAC(co, rnorm(5), rnorm(20), h = 0.1, perms = 5, verbose = FALSE))

  # Zero variance is the deliberate exception: the reference returns an all-zero curve rather than
  # failing, so multiScaleAC warns and reproduces that.
  expect_false(validate_inputs(co, rep(1, 20), rnorm(20))$valid)
  expect_warning(res <- multiScaleAC(co, rep(1, 20), rnorm(20), h = 0.1, perms = 5,
                                     verbose = FALSE),
                 regexp = "zero variance", ignore.case = TRUE)
  expect_identical(res@observed, 0)
})

test_that("suggest_bandwidths returns genuinely different scales per range", {
  # Every range setting previously returned quantiles of the SAME nearest-neighbour distribution, so
  # "global" was not global; all three gave near-nearest-neighbour bandwidths.
  set.seed(5); co <- matrix(runif(2 * 400), ncol = 2)
  loc <- suggest_bandwidths(co, 5, range = "local")
  med <- suggest_bandwidths(co, 5, range = "medium")
  glb <- suggest_bandwidths(co, 5, range = "global")

  expect_length(loc, 5L); expect_length(med, 5L); expect_length(glb, 5L)
  expect_true(all(diff(loc) > 0)); expect_true(all(diff(glb) > 0))

  # "global" must reach a scale that "local" never does -- by a wide margin, not marginally.
  expect_gt(max(glb) / max(loc), 5)
  expect_gt(max(med), max(loc))
})

test_that("suggest_bandwidths produces bandwidths that are actually usable", {
  # The real test of a bandwidth suggestion: feeding it back in must give sensible neighbourhood
  # sizes, neither empty nor saturated.
  set.seed(6); co <- matrix(runif(2 * 500), ncol = 2)
  h <- suggest_bandwidths(co, 6, range = "medium")
  qc <- qc_neighborhoods(co, h)
  expect_true(all(qc$saturation < MAX_SATURATION))
  expect_gt(max(qc$mean_neighbors), 5)      # the wide end sees a real neighbourhood
  expect_lt(min(qc$mean_neighbors), 20)     # the narrow end is genuinely local
})

test_that("suggest_bandwidths reports the scale it derived from", {
  set.seed(7); co <- matrix(runif(2 * 300), ncol = 2)
  h <- suggest_bandwidths(co, 4)
  ref <- attr(h, "scale_reference")
  expect_true(all(c("median_nn_distance", "extent") %in% names(ref)))
  # Sanity: on a unit square with 300 points the median nearest-neighbour distance is a few
  # hundredths, and the diagonal extent is near sqrt(2).
  expect_lt(ref[["median_nn_distance"]], 0.15)
  expect_gt(ref[["extent"]], 1.0)
})

test_that("suggest_bandwidths gives millimetre bandwidths for the millimetre example data", {
  # The unit trap in concrete form. ccRCCExample is in mm, so suggestions must come back on that
  # scale -- not the h = 10..250 that earlier code in this repository used.
  skip_if_not(exists("ccRCCExample"), "example data unavailable")
  data(ccRCCExample, envir = environment())
  co <- as.matrix(ccRCCExample[, c("x", "y")])
  h <- suggest_bandwidths(co, 6, range = "medium")
  expect_true(all(h < 1))
  expect_true(all(h > 1e-4))
})

test_that("nn_distance matches an explicit calculation", {
  set.seed(8); co <- matrix(runif(2 * 200), ncol = 2)
  for (dm in c("euclidean", "manhattan")) {
    d <- as.matrix(stats::dist(co, method = dm)); diag(d) <- Inf
    expect_equal(as.vector(multiScaleAC:::nn_distance(co, dm)),
                 unname(apply(d, 1, min)), tolerance = 1e-12, info = dm)
  }
})

test_that("qc_neighborhoods accepts a bandwidth vector", {
  # A vector previously recycled against the distance matrix, emitting only a recycling warning and
  # returning meaningless numbers.
  set.seed(9); co <- matrix(runif(2 * 200), ncol = 2)
  qc <- qc_neighborhoods(co, h = c(0.05, 0.1, 0.2))
  expect_s3_class(qc, "data.frame")
  expect_identical(nrow(qc), 3L)
  expect_identical(qc$h, c(0.05, 0.1, 0.2))
  # Neighbourhood size must increase with bandwidth.
  expect_true(all(diff(qc$mean_neighbors) > 0))
})

test_that("qc_neighborhoods reports effective weight for gaussian, counts for uniform", {
  # The old version applied a hard distance <= h cutoff regardless of kernel, so under the default
  # Gaussian it could never report an isolated cell and described a model that was not being fitted.
  set.seed(10); co <- matrix(runif(2 * 150), ncol = 2)

  # For the uniform kernel the row sum IS a neighbour count, so it must be a whole number.
  qu <- qc_neighborhoods(co, h = 0.1, kernel = "uniform")
  expect_equal(qu$mean_neighbors * 150, round(qu$mean_neighbors * 150), tolerance = 1e-8)

  # For the Gaussian it is a weight, generally fractional, and strictly larger than the hard count
  # at the same h because distant pairs still contribute.
  qg <- qc_neighborhoods(co, h = 0.1, kernel = "gaussian")
  expect_gt(qg$mean_neighbors, qu$mean_neighbors)
})

test_that("qc_neighborhoods exposes saturation, the failure mode that broke earlier validation", {
  # This is the diagnostic that makes the degenerate regime visible. At h = 25 on unit-square
  # coordinates every pairwise weight is ~1, the kernel is constant, and Moran's I collapses.
  set.seed(11); co <- matrix(runif(2 * 200), ncol = 2)
  qc <- qc_neighborhoods(co, h = c(0.05, 25))
  expect_lt(qc$saturation[1], 0.2)
  expect_gt(qc$saturation[2], 0.99)

  # And confirm the consequence, so the diagnostic is tied to the thing it predicts.
  g1 <- rnorm(200); g2 <- rnorm(200)
  res <- multiScaleAC(co, g1, g2, h = c(0.05, 25), perms = 20, verbose = FALSE)
  expect_gt(abs(res@observed[1]), 10 * abs(res@observed[2]))
})

test_that("qc_neighborhoods rejects invalid bandwidths", {
  set.seed(12); co <- matrix(runif(40), ncol = 2)
  expect_error(qc_neighborhoods(co, h = 0), regexp = "positive")
  expect_error(qc_neighborhoods(co, h = -1), regexp = "positive")
  expect_error(qc_neighborhoods(co, h = NA_real_), regexp = "positive|finite")
  expect_error(qc_neighborhoods(1:10, h = 0.1), regexp = "matrix")
})

test_that("show and print agree, and neither dumps the permutation matrix", {
  # There was no show method, so auto-printing at the prompt fell through to the default S4 show,
  # which prints every slot including the full permutation and coordinate matrices.
  set.seed(13); n <- 100
  co <- matrix(runif(2 * n), ncol = 2)
  res <- multiScaleAC(co, rnorm(n), rnorm(n), h = c(0.05, 0.1), perms = 50, verbose = FALSE)

  shown <- capture.output(show(res))
  printed <- capture.output(print(res))
  expect_identical(shown, printed)
  # A handful of summary lines, not hundreds of numbers.
  expect_lt(length(shown), 15L)
  expect_true(any(grepl("multiScaleAC result", shown, fixed = TRUE)))
})
