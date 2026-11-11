# Boundary conditions. Each is checked against the reference implementation where the reference is
# defined, so "what should happen" is not a matter of opinion.

ref_obs <- function(coords, g1, g2, h, perms = 0, seed = 1) {
  soupir_moran_multi(coords, g1, g2, h_vec = h, perms = max(1L, perms), seed = seed)$observed
}

test_that("the minimum viable sample size works", {
  set.seed(1)
  co <- matrix(c(0, 0, 1, 0, 0.5, 1), ncol = 2, byrow = TRUE)   # N = 3
  g1 <- c(1, 2, 3); g2 <- c(3, 1, 2)
  res <- multiScaleAC(co, g1, g2, h = 1, perms = 5, verbose = FALSE)
  expect_length(res@observed, 1L)
  expect_equal(res@observed, ref_obs(co, g1, g2, 1), tolerance = REFERENCE_TOLERANCE)

  # Below the minimum must fail with an explanation, not a subscript error.
  co2 <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  expect_error(multiScaleAC(co2, c(1, 2), c(2, 1), h = 1, perms = 5, verbose = FALSE),
               regexp = "3 observations", ignore.case = TRUE)
})

test_that("a single permutation is handled", {
  set.seed(2); n <- 60
  co <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  res <- multiScaleAC(co, g1, g2, h = c(0.1, 0.2), perms = 1, seed = 5, verbose = FALSE)
  expect_identical(dim(res@permutations), c(2L, 1L))
  expect_true(all(is.finite(res@permutations)))
  expect_equal(res@observed, ref_obs(co, g1, g2, c(0.1, 0.2)),
               tolerance = REFERENCE_TOLERANCE)
})

test_that("duplicated coordinates are warned about but computed", {
  # Coincident cells give each other the maximum kernel weight, which is worth flagging, but it is
  # not an error and the reference handles it.
  set.seed(3); n <- 80
  co <- matrix(runif(2 * n), ncol = 2)
  co[2, ] <- co[1, ]
  co[7, ] <- co[6, ]
  g1 <- rnorm(n); g2 <- rnorm(n)
  expect_warning(res <- multiScaleAC(co, g1, g2, h = 0.1, perms = 5, verbose = FALSE),
                 regexp = "duplicat", ignore.case = TRUE)
  expect_equal(res@observed, ref_obs(co, g1, g2, 0.1), tolerance = REFERENCE_TOLERANCE)
})

test_that("all-identical coordinates do not crash", {
  # Every distance is 0, so every off-diagonal weight is 1 and the kernel is exactly constant.
  set.seed(4); n <- 20
  co <- matrix(rep(c(0.5, 0.5), each = n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  res <- suppressWarnings(multiScaleAC(co, g1, g2, h = 0.1, perms = 5, verbose = FALSE))
  expect_true(is.finite(res@observed) || is.na(res@observed))
  expect_equal(res@observed, suppressWarnings(ref_obs(co, g1, g2, 0.1)),
               tolerance = REFERENCE_TOLERANCE)
})

test_that("an extreme bandwidth range is handled at both ends", {
  set.seed(5); n <- 100
  co <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  # 1e-6 underflows to no neighbours; 1e6 saturates completely.
  res <- suppressWarnings(
    multiScaleAC(co, g1, g2, h = c(1e-6, 0.1, 1e6), perms = 5, verbose = FALSE)
  )
  expect_length(res@observed, 3L)
  expect_true(all(is.finite(res@observed)))
  expect_equal(res@observed, suppressWarnings(ref_obs(co, g1, g2, c(1e-6, 0.1, 1e6))),
               tolerance = REFERENCE_TOLERANCE)
})

test_that("bandwidths are reported in the order supplied", {
  # An unsorted vector must not be silently reordered, or a user's h and their curve stop
  # corresponding.
  set.seed(6); n <- 80
  co <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  hv <- c(0.2, 0.05, 0.1)
  res <- multiScaleAC(co, g1, g2, h = hv, perms = 5, verbose = FALSE)
  expect_identical(res@h, hv)
  expect_equal(res@observed, ref_obs(co, g1, g2, hv), tolerance = REFERENCE_TOLERANCE)
})

test_that("non-positive bandwidths are dropped with a warning", {
  set.seed(7); n <- 60
  co <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  expect_warning(res <- multiScaleAC(co, g1, g2, h = c(0.1, 0, 0.2), perms = 5, verbose = FALSE),
                 regexp = "zero bandwidth", ignore.case = TRUE)
  expect_identical(res@h, c(0.1, 0.2))
  expect_warning(multiScaleAC(co, g1, g2, h = c(-1, 0.1), perms = 5, verbose = FALSE),
                 regexp = "negative", ignore.case = TRUE)
  expect_error(suppressWarnings(multiScaleAC(co, g1, g2, h = c(-1, 0), perms = 5, verbose = FALSE)),
               regexp = "No valid bandwidth", ignore.case = TRUE)
})

test_that("the result object satisfies its own validity rules", {
  set.seed(8); n <- 70L
  co <- matrix(runif(2 * n), ncol = 2)
  res <- multiScaleAC(co, rnorm(n), rnorm(n), h = c(0.1, 0.2), perms = 11, verbose = FALSE)
  expect_true(methods::validObject(res))
  expect_identical(res@n_cells, n)
  expect_identical(nrow(res@coords), n)
  expect_identical(length(res@h), length(res@observed))
  expect_identical(nrow(res@permutations), length(res@h))
  # Significance slots must all still be empty on a fresh result.
  for (s in c("p_values", "z_scores", "perm_mean", "perm_sd",
              "null_lower", "null_upper", "significant", "p_adjusted")) {
    expect_length(methods::slot(res, s), 0L)
  }
})

test_that("a different seed changes the permutations but not the observed curve", {
  set.seed(9); n <- 80
  co <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  a <- multiScaleAC(co, g1, g2, h = c(0.1, 0.2), perms = 30, seed = 1, verbose = FALSE)
  b <- multiScaleAC(co, g1, g2, h = c(0.1, 0.2), perms = 30, seed = 2, verbose = FALSE)
  expect_identical(a@observed, b@observed)
  expect_false(identical(a@permutations, b@permutations))
})

test_that("multiScaleAC does not disturb the caller's random stream unpredictably", {
  # It calls set.seed() internally, which is the reference's behaviour and is what makes results
  # reproducible. Document that consequence with a test rather than leaving it to be discovered:
  # the stream afterwards is a deterministic function of `seed`, not of the caller's state.
  set.seed(123); n <- 50
  co <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)

  set.seed(999)
  invisible(multiScaleAC(co, g1, g2, h = 0.1, perms = 5, seed = 7, verbose = FALSE))
  after_a <- runif(3)

  set.seed(111)   # deliberately different caller state
  invisible(multiScaleAC(co, g1, g2, h = 0.1, perms = 5, seed = 7, verbose = FALSE))
  after_b <- runif(3)

  expect_identical(after_a, after_b)
})
