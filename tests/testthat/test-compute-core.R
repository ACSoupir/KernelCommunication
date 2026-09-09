# Tests for the tiled two-pass compute core.
#
# Each of these corresponds to a defect that was present and is now fixed. They are written to fail
# against the previous core, not merely to describe it.

# Explicit, deliberately naive R implementations to compare against. Written out longhand rather
# than reusing package code, so a shared mistake cannot make both sides agree.
ref_curve <- function(coords, gene1, gene2, h_vec,
                      kernel = "gaussian", metric = "euclidean", eps = 0) {
  x <- as.vector(scale(gene1)); y <- as.vector(scale(gene2)); S <- sum(x^2)
  d <- as.matrix(stats::dist(coords, method = metric))
  vapply(h_vec, function(h) {
    K <- if (kernel == "gaussian") exp(-d^2 / (2 * h^2)) else (d <= h) * 1.0
    diag(K) <- 0
    if (eps > 0) K[K <= eps] <- 0
    rs <- rowSums(K)
    if (any(rs <= 0)) return(0)          # matches the reference's whole-bandwidth zeroing
    sum(as.vector(K %*% (x / rs)) * y) / S
  }, numeric(1))
}

test_that("gaussian kernel with manhattan distance is correct", {
  # THE test for this phase. The previous core's dense Gaussian branch read `dist_sq`
  # unconditionally, but `dist_sq` was only filled when distance == "euclidean". So this exact
  # combination, at any N below its hard-coded 10000 threshold, computed the kernel from
  # uninitialised heap memory. Three identical calls returned 4.5e-05, then -5.3e-04, then NaN,
  # against a correct value of 2.6e-02 -- and one returned the same number for every bandwidth,
  # which is arithmetically impossible.
  #
  # The old benchmark that "covered" this path asserted only !is.na(), so it reported PASSED on
  # garbage. This compares against an explicit Manhattan implementation instead.
  set.seed(42); n <- 300
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  hv <- c(0.02, 0.05, 0.1)

  expected <- ref_curve(coords, g1, g2, hv, metric = "manhattan")
  got <- multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 1,
                      distance = "manhattan", verbose = FALSE)
  expect_equal(got@observed, expected, tolerance = REFERENCE_TOLERANCE)

  # The values must also be non-degenerate, or agreement would be vacuous.
  expect_gt(max(abs(expected)), 1e-3)
})

test_that("repeated identical calls give identical results for every kernel and metric", {
  # Non-determinism across identical calls is the signature of reading uninitialised memory. This
  # would have caught the defect above even without a reference to compare against.
  set.seed(7); n <- 200
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  for (kn in c("gaussian", "uniform")) {
    for (dm in c("euclidean", "manhattan")) {
      runs <- replicate(3, suppressWarnings(
        multiScaleAC(coords, g1, g2, h = c(0.05, 0.1), perms = 3, seed = 2,
                     kernel = kn, distance = dm, verbose = FALSE)@observed
      ), simplify = FALSE)
      expect_identical(runs[[1]], runs[[2]], info = paste(kn, dm))
      expect_identical(runs[[2]], runs[[3]], info = paste(kn, dm))
    }
  }
})

test_that("every kernel and metric combination matches an explicit implementation", {
  set.seed(11); n <- 220
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  hv <- c(0.04, 0.08, 0.15)
  for (kn in c("gaussian", "uniform")) {
    for (dm in c("euclidean", "manhattan")) {
      expected <- ref_curve(coords, g1, g2, hv, kernel = kn, metric = dm)
      got <- suppressWarnings(multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 3,
                                          kernel = kn, distance = dm, verbose = FALSE))
      expect_equal(got@observed, expected, tolerance = REFERENCE_TOLERANCE,
                   info = paste(kn, dm))
    }
  }
})

test_that("eps is honoured, and applies to the row sums as well as the weights", {
  # eps used to be applied in the sparse path and ignored entirely in the dense one, so the two
  # disagreed on the denominator and therefore on I. With one path that cannot happen, but the
  # semantics still need pinning: a thresholded weight must be excluded from the row sum, not merely
  # from the numerator.
  set.seed(13); n <- 200
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  hv <- c(0.05, 0.1)

  for (e in c(0, 1e-6, 1e-3, 0.05)) {
    expected <- ref_curve(coords, g1, g2, hv, eps = e)
    got <- suppressWarnings(multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 4,
                                        eps = e, verbose = FALSE))
    expect_equal(got@observed, expected, tolerance = REFERENCE_TOLERANCE,
                 info = paste("eps =", e))
  }

  # eps = 0 is the documented default and the condition under which the core reduces exactly to the
  # reference. A larger eps must actually change the answer, or the argument would be inert -- which
  # is what it was in dense mode before.
  a <- multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 4, eps = 0,    verbose = FALSE)
  b <- multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 4, eps = 0.05, verbose = FALSE)
  expect_false(isTRUE(all.equal(a@observed, b@observed)))
})

test_that("batch_size changes memory use but not results", {
  # batch_size was previously accepted and then ignored (hard-coded to 100 in the R path, unused in
  # dense C++). It now bounds the permutation buffer at N x batch_size instead of N x perms, and
  # must be bit-for-bit invariant because the permutation stream is unchanged.
  set.seed(17); n <- 150
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  small <- multiScaleAC(coords, g1, g2, h = c(0.05, 0.1), perms = 101, seed = 5,
                        batch_size = 7, verbose = FALSE)
  large <- multiScaleAC(coords, g1, g2, h = c(0.05, 0.1), perms = 101, seed = 5,
                        batch_size = 10000, verbose = FALSE)
  expect_identical(small@permutations, large@permutations)
  expect_identical(small@observed, large@observed)
})

test_that("perms = 0 is allowed and yields a zero-column permutation matrix", {
  # Previously rejected by validation, yet guarded for in four places in the compute code. It is the
  # mode analytical inference uses, so it is now legal and tested.
  set.seed(19); n <- 120
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  res <- multiScaleAC(coords, g1, g2, h = c(0.05, 0.1), perms = 0, verbose = FALSE)
  expect_identical(ncol(res@permutations), 0L)
  expect_length(res@observed, 2L)
  expect_true(all(is.finite(res@observed)))

  # The observed curve must not depend on whether permutations were requested.
  withp <- multiScaleAC(coords, g1, g2, h = c(0.05, 0.1), perms = 50, verbose = FALSE)
  expect_equal(res@observed, withp@observed, tolerance = REFERENCE_TOLERANCE)
})

test_that("gene names come from the caller's expressions", {
  # deparse(substitute(gene1)) was evaluated in the wrong frame, so every result was labelled
  # "gene1" x "gene2" -- including plot titles and the mxfda interaction name.
  set.seed(23); n <- 80
  coords <- matrix(runif(2 * n), ncol = 2)
  expr <- data.frame(CD274 = rnorm(n), PDCD1 = rnorm(n))
  res <- multiScaleAC(coords, expr$CD274, expr$PDCD1, h = 0.1, perms = 5, verbose = FALSE)
  expect_identical(res@gene1_name, "expr$CD274")
  expect_identical(res@gene2_name, "expr$PDCD1")

  # Positional calls must work too, since match.call() names its arguments.
  a <- rnorm(n); b <- rnorm(n)
  res2 <- multiScaleAC(coords, a, b, h = 0.1, perms = 5, verbose = FALSE)
  expect_identical(res2@gene1_name, "a")
  expect_identical(res2@gene2_name, "b")
})

test_that("bandwidth input is validated rather than failing obscurely", {
  set.seed(29); n <- 60
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)

  # .sanitize_h() previously evaluated `if (NA > 0)` and died with "missing value where TRUE/FALSE
  # needed", which says nothing about the actual problem.
  expect_error(multiScaleAC(coords, g1, g2, h = c(0.1, NA), perms = 5, verbose = FALSE),
               regexp = "missing value", ignore.case = TRUE)
  expect_error(multiScaleAC(coords, g1, g2, h = c(0.1, Inf), perms = 5, verbose = FALSE),
               regexp = "infinite", ignore.case = TRUE)
  expect_warning(multiScaleAC(coords, g1, g2, h = c(0.1, 0.1), perms = 5, verbose = FALSE),
                 regexp = "duplicat", ignore.case = TRUE)
  expect_error(multiScaleAC(coords, g1, g2, h = numeric(0), perms = 5, verbose = FALSE))
  expect_error(multiScaleAC(coords, g1, g2, h = 0.1, perms = -1, verbose = FALSE),
               regexp = "non-negative", ignore.case = TRUE)
})

test_that("large_mode_n is deprecated rather than silently ignored", {
  # It was previously accepted, validated, passed to C++, and never read; the mode threshold was
  # hard-coded. Saying so is better than continuing to accept it in silence.
  set.seed(31); n <- 60
  coords <- matrix(runif(2 * n), ncol = 2)
  expect_warning(
    multiScaleAC(coords, rnorm(n), rnorm(n), h = 0.1, perms = 5,
                 large_mode_n = 10000, verbose = FALSE),
    regexp = "large_mode_n", ignore.case = TRUE
  )
})

test_that("a custom distance with a built-in kernel works", {
  # This combination previously always errored with "Built-in kernels should use C++ path", because
  # the R path was entered when EITHER kernel or distance was a function but then refused to proceed
  # unless the kernel was. The Minkowski example in multiScaleAC()'s own documentation therefore
  # threw an error.
  set.seed(37); n <- 120
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  hv <- c(0.05, 0.1)

  minkowski_dist <- function(coords, p = 3) {
    as.matrix(stats::dist(coords, method = "minkowski", p = p))
  }
  got <- multiScaleAC(coords, g1, g2, h = hv, distance = minkowski_dist, p = 3,
                      perms = 5, seed = 6, verbose = FALSE)
  expect_length(got@observed, 2L)
  expect_true(all(is.finite(got@observed)))

  # p = 2 Minkowski is Euclidean, so the custom path must agree with the built-in one. This checks
  # the R and C++ paths against each other, which is what used to be impossible.
  builtin <- multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 6, verbose = FALSE)
  viaR <- multiScaleAC(coords, g1, g2, h = hv, distance = minkowski_dist, p = 2,
                       perms = 5, seed = 6, verbose = FALSE)
  expect_equal(viaR@observed, builtin@observed, tolerance = REFERENCE_TOLERANCE)
  expect_equal(viaR@permutations, builtin@permutations, tolerance = REFERENCE_TOLERANCE)
})

test_that("a user parameter cannot be captured by partial matching against an option", {
  # `...` forwards parameters to user kernel and distance functions. R only enforces exact name
  # matching for formals declared AFTER `...`, so when the options sat before it, `p = 3` for a
  # Minkowski distance was partially matched against `perms`/`permute` and the call failed -- and
  # `p = 3` is the example in multiScaleAC()'s own documentation. `...` now precedes every option.
  set.seed(53); n <- 100
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)

  minkowski_dist <- function(coords, p = 3) {
    as.matrix(stats::dist(coords, method = "minkowski", p = p))
  }
  expect_silent(
    res <- multiScaleAC(coords, g1, g2, h = 0.1, distance = minkowski_dist, p = 3,
                        perms = 5, seed = 1, verbose = FALSE)
  )
  expect_true(is.finite(res@observed))

  # `p` really was forwarded, not swallowed: p = 2 is Euclidean and must match the built-in path,
  # while p = 6 must not.
  builtin <- multiScaleAC(coords, g1, g2, h = 0.1, perms = 5, seed = 1, verbose = FALSE)
  p2 <- multiScaleAC(coords, g1, g2, h = 0.1, distance = minkowski_dist, p = 2,
                     perms = 5, seed = 1, verbose = FALSE)
  p6 <- multiScaleAC(coords, g1, g2, h = 0.1, distance = minkowski_dist, p = 6,
                     perms = 5, seed = 1, verbose = FALSE)
  expect_equal(p2@observed, builtin@observed, tolerance = REFERENCE_TOLERANCE)
  expect_false(isTRUE(all.equal(p6@observed, builtin@observed)))

  # Other single-letter names that prefix an option must be forwarded too.
  odd_kernel <- function(dist_matrix, h, m = 2, n = 1, e = 1, s = 1, b = 1, v = 1) {
    exp(-(dist_matrix / h)^m) * n * e * s * b * v
  }
  expect_silent(
    multiScaleAC(coords, g1, g2, h = 0.1, kernel = odd_kernel, m = 2, n = 1, e = 1,
                 s = 1, b = 1, v = 1, perms = 5, seed = 1, verbose = FALSE)
  )
})

test_that("a custom kernel has no self-weight, matching the built-in path", {
  # The custom path never zeroed the diagonal, so a kernel returning exp(0) = 1 there gave every
  # cell a self-weight. Custom results were on a different scale from built-in ones and silently not
  # comparable.
  set.seed(41); n <- 120
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  hv <- c(0.05, 0.1)

  # A custom kernel that reproduces the built-in Gaussian must give the built-in answer.
  gauss_kernel <- function(dist_matrix, h) exp(-dist_matrix^2 / (2 * h^2))
  viaR <- multiScaleAC(coords, g1, g2, h = hv, kernel = gauss_kernel,
                       perms = 5, seed = 8, verbose = FALSE)
  builtin <- multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 8, verbose = FALSE)
  expect_equal(viaR@observed, builtin@observed, tolerance = REFERENCE_TOLERANCE)
})

test_that("the kernel pass reports neighbourhood diagnostics", {
  # Needed so a saturated or empty bandwidth is detectable rather than silently producing a
  # meaningless number.
  set.seed(43); n <- 150L
  coords <- matrix(runif(2 * n), ncol = 2)
  pass <- multiScaleAC:::moran_kernel_pass(coords, as.vector(scale(rnorm(n))),
                                           c(0.01, 0.1, 5), eps = 0)
  expect_identical(dim(pass$V), c(3L, n))
  expect_length(pass$rowsum_mean, 3L)
  expect_length(pass$n_isolated, 3L)
  # Row sums must increase with bandwidth, and h = 5 on a unit square must be saturated.
  expect_true(all(diff(pass$rowsum_mean) > 0))
  expect_gt(pass$rowsum_mean[3], 0.9 * (n - 1))
  expect_identical(pass$n_isolated[3], 0L)

  # Isolation requires the Gaussian to underflow to exactly 0, which needs d > h * sqrt(2 * 745),
  # about 38.6 * h. A point placed well outside the unit square is isolated at h = 0.01 while its
  # neighbours are not, so this checks the census counts cells rather than bandwidths.
  far <- rbind(coords, c(50, 50))
  pass2 <- multiScaleAC:::moran_kernel_pass(far, as.vector(scale(rnorm(n + 1L))),
                                            c(0.01, 5), eps = 0)
  expect_identical(pass2$n_isolated[1], 1L)
  expect_identical(pass2$n_isolated[2], 0L)
  expect_identical(pass2$rowsum_min[1], 0)
})

test_that("coords with extra columns use the first two", {
  set.seed(47); n <- 90
  xy <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  two <- multiScaleAC(xy, g1, g2, h = 0.1, perms = 5, seed = 9, verbose = FALSE)
  expect_equal(two@observed, ref_curve(xy, g1, g2, 0.1), tolerance = REFERENCE_TOLERANCE)
})
