# Every built-in kernel and distance combination, against an explicit R implementation written out
# longhand so a shared mistake cannot make both sides agree.
#
# All of these run inside the same tiled core, which is the point: a kernel added as a separate code
# path could bypass the eps predicate or the zero-row-sum policy, which is exactly how the previous
# implementation's divergences arose.

ref_any <- function(coords, g1, g2, h_vec, kernel, metric, p = 2, eps = 0) {
  x <- as.vector(scale(g1)); y <- as.vector(scale(g2)); S <- sum(x^2)
  d <- if (metric == "minkowski") {
    as.matrix(stats::dist(coords, method = "minkowski", p = p))
  } else {
    as.matrix(stats::dist(coords, method = metric))
  }
  vapply(h_vec, function(h) {
    K <- switch(kernel,
                gaussian    = exp(-d^2 / (2 * h^2)),
                uniform     = (d <= h) * 1.0,
                exponential = exp(-d / h))
    diag(K) <- 0
    if (eps > 0) K[K <= eps] <- 0
    rs <- rowSums(K)
    if (any(rs <= 0)) return(0)
    sum(as.vector(K %*% (x / rs)) * y) / S
  }, numeric(1))
}

KERNELS <- c("gaussian", "uniform", "exponential")
METRICS <- c("euclidean", "manhattan")

test_that("every kernel and metric combination matches an explicit implementation", {
  set.seed(1); n <- 220
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  hv <- c(0.04, 0.08, 0.15)

  for (kn in KERNELS) for (dm in METRICS) {
    got <- suppressWarnings(
      multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 2,
                   kernel = kn, distance = dm, verbose = FALSE)
    )
    expect_equal(got@observed, ref_any(coords, g1, g2, hv, kn, dm),
                 tolerance = REFERENCE_TOLERANCE, info = paste(kn, dm))
  }
})

test_that("the Minkowski distance is correct and reduces to the special cases", {
  set.seed(2); n <- 200
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  hv <- c(0.05, 0.1)

  for (pp in c(1, 1.5, 2, 3, 6)) {
    got <- multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 3,
                        distance = "minkowski", p = pp, verbose = FALSE)
    expect_equal(got@observed, ref_any(coords, g1, g2, hv, "gaussian", "minkowski", p = pp),
                 tolerance = REFERENCE_TOLERANCE, info = paste("p =", pp))
  }

  # p = 1 is Manhattan and p = 2 is Euclidean. Those are routed to the cheap paths, so this checks the
  # routing gives the same answer rather than a subtly different one.
  euc <- multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 3,
                      distance = "euclidean", verbose = FALSE)
  mk2 <- multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 3,
                      distance = "minkowski", p = 2, verbose = FALSE)
  expect_identical(mk2@observed, euc@observed)

  man <- multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 3,
                      distance = "manhattan", verbose = FALSE)
  mk1 <- multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 3,
                      distance = "minkowski", p = 1, verbose = FALSE)
  expect_identical(mk1@observed, man@observed)
})

test_that("an invalid Minkowski order is rejected", {
  set.seed(3); n <- 80
  coords <- matrix(runif(2 * n), ncol = 2)
  expect_error(multiScaleAC(coords, rnorm(n), rnorm(n), h = 0.1, perms = 5,
                            distance = "minkowski", p = 0.5, verbose = FALSE),
               regexp = "p.*1", ignore.case = TRUE)
})

test_that("the exponential and Gaussian kernels have the shape relationship they should", {
  # A behavioural check, and one worth stating precisely because the intuitive version is wrong.
  #
  # exp(-d/h) < exp(-d^2/2h^2) for d < 2h, and greater beyond it: the exponential decays FASTER near
  # the origin and only has the heavier tail further out, crossing over at exactly d = 2h. So
  # "exponential has a heavier tail" does NOT mean it gives a wider neighbourhood at the same h --
  # the near field dominates the row sum, and the exponential's is smaller. An earlier version of
  # this test asserted the opposite and failed, which is why it is spelled out here.
  h <- 0.05
  d <- c(0.1 * h, 0.5 * h, h, 1.9 * h, 2 * h, 3 * h, 6 * h)
  g <- exp(-d^2 / (2 * h^2))
  e <- exp(-d / h)
  expect_true(all(e[d < 2 * h] < g[d < 2 * h]))
  expect_true(all(e[d > 2 * h] > g[d > 2 * h]))
  expect_equal(e[d == 2 * h], g[d == 2 * h])       # crossover

  # And the consequence, so the kernels cannot be silently transposed in the C++.
  set.seed(4); n <- 200
  coords <- matrix(runif(2 * n), ncol = 2)
  qg <- qc_neighborhoods(coords, h = h, kernel = "gaussian")
  qe <- qc_neighborhoods(coords, h = h, kernel = "exponential")
  expect_lt(qe$mean_neighbors, qg$mean_neighbors)
})

test_that("every built-in combination honours eps identically", {
  # A new kernel must not be able to bypass the eps predicate, which is the reason they all live in
  # one core rather than in separate branches.
  set.seed(5); n <- 180
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  hv <- c(0.05, 0.1)

  for (kn in KERNELS) for (e in c(0, 1e-4, 0.05)) {
    got <- suppressWarnings(
      multiScaleAC(coords, g1, g2, h = hv, perms = 5, seed = 6,
                   kernel = kn, eps = e, verbose = FALSE)
    )
    expect_equal(got@observed,
                 ref_any(coords, g1, g2, hv, kn, "euclidean", eps = e),
                 tolerance = REFERENCE_TOLERANCE, info = paste(kn, "eps =", e))
  }
})

test_that("every built-in combination honours the zero-row-sum policy", {
  # A uniform kernel at a small bandwidth genuinely isolates cells, which is the case the policy
  # exists for. Every kernel must route through it.
  set.seed(6); n <- 100
  coords <- rbind(matrix(runif(2 * (n - 1)), ncol = 2), c(50, 50))
  g1 <- rnorm(n); g2 <- rnorm(n)

  for (kn in KERNELS) {
    expect_warning(
      res <- multiScaleAC(coords, g1, g2, h = 0.05, perms = 5, seed = 7,
                          kernel = kn, verbose = FALSE),
      regexp = "no neighbours|no neighbors", ignore.case = TRUE, info = kn
    )
    expect_true(all(is.finite(res@observed)), info = kn)
    expect_identical(res@observed, 0, info = kn)
  }
})

test_that("exact moments are available for every kernel", {
  # The moment identities require row sums of 1, which every built-in kernel provides. If a kernel
  # were added that could not be normalised, this would catch it.
  set.seed(7); n <- 200
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  # Bandwidths wide enough that no cell is isolated under any of the three kernels, since an
  # isolated cell zeroes the whole bandwidth and there would then be no distribution to compare.
  for (kn in KERNELS) {
    res <- multiScaleAC(coords, g1, g2, h = c(0.15, 0.3), perms = 3000, seed = 8,
                        kernel = kn, verbose = FALSE)
    expect_true(all(res@observed != 0), info = paste(kn, "must not be a zeroed bandwidth"))
    emp <- apply(res@permutations, 1, stats::var)
    expect_equal(res@moment_var, emp, tolerance = 0.15, info = kn)
  }
})

test_that("a custom kernel function still reaches the R fallback and agrees with the core", {
  # The R path exists for arbitrary user functions. A custom function reproducing a built-in kernel
  # must give the built-in answer, which is what makes the two paths comparable at all.
  set.seed(8); n <- 150
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  hv <- c(0.05, 0.1)

  exp_kernel <- function(dist_matrix, h) exp(-dist_matrix / h)
  viaR <- multiScaleAC(coords, g1, g2, h = hv, kernel = exp_kernel,
                       perms = 5, seed = 9, verbose = FALSE)
  builtin <- multiScaleAC(coords, g1, g2, h = hv, kernel = "exponential",
                          perms = 5, seed = 9, verbose = FALSE)
  expect_equal(viaR@observed, builtin@observed, tolerance = REFERENCE_TOLERANCE)
})
