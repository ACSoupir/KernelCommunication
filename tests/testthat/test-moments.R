# Exact moments of the permutation null.
#
# These formulas are the substance of the permutation-reduction work, so they are re-derived against
# brute-force permutations on every run rather than trusted. The reference values are computed by
# explicitly permuting, which is slow but unambiguous.

test_that("the exact moments match brute-force permutation moments", {
  # Enough permutations that the Monte Carlo error on the variance is well under the tolerance,
  # and on a case with real skewness so the third moment is resolvable.
  set.seed(1)
  n <- 150
  coords <- rbind(matrix(rnorm(2 * (n - 15), 0.5, 0.12), ncol = 2),
                  matrix(runif(30), ncol = 2))
  g1 <- rlnorm(n, 0, 1)
  g2 <- c(rep(0, n - 25), rlnorm(25, 1.5, 1))[sample.int(n)]
  hv <- c(0.03, 0.1)

  res <- multiScaleAC(coords, g1, g2, h = hv, perms = 40000, seed = 3, verbose = FALSE)

  for (i in seq_along(hv)) {
    null <- res@permutations[i, ]
    emp_var <- stats::var(null)
    emp_skew <- mean((null - mean(null))^3) / stats::sd(null)^3

    # Monte Carlo error on a variance from M draws is about sqrt(2/M) relative; at M = 40000 that is
    # 0.7%, so 5% is a comfortable but still meaningful bound.
    expect_equal(res@moment_var[i], emp_var, tolerance = 0.05,
                 info = paste("variance, h =", hv[i]))
    # Skewness is noisier; assert it where it is large enough to be resolved.
    if (abs(res@moment_skewness[i]) > 0.3) {
      expect_equal(res@moment_skewness[i], emp_skew, tolerance = 0.15,
                   info = paste("skewness, h =", hv[i]))
    }
  }

  # The exact mean is 0. The empirical mean is an estimate of it, so it must sit inside its own
  # Monte Carlo error rather than being equal.
  mc_se <- sqrt(res@moment_var / ncol(res@permutations))
  expect_true(all(abs(rowMeans(res@permutations)) < 4 * mc_se))
})

test_that("the closed-form cross-bandwidth covariance matches the empirical one", {
  set.seed(2)
  n <- 150
  coords <- matrix(runif(2 * n), ncol = 2)
  res <- multiScaleAC(coords, rnorm(n), rnorm(n), h = c(0.05, 0.1, 0.2),
                      perms = 30000, seed = 4, verbose = FALSE)
  emp <- stats::cov(t(res@permutations))
  expect_equal(dim(res@moment_cov), c(3L, 3L))
  expect_equal(res@moment_cov, emp, tolerance = 0.08, ignore_attr = TRUE)
  # The diagonal must agree with the separately computed variances.
  expect_equal(diag(res@moment_cov), res@moment_var, tolerance = 1e-10)
})

test_that("the exact variance reduces to the known special case", {
  # With both genes z-scored, Var(I) = ||t(W) x||^2 / (n-1)^2. That is the formula the retracted
  # analytical-variance note gave; the mathematics there was correct even though its evidence was
  # not (docs/RETRACTIONS.md section 4). Confirming it keeps that distinction honest.
  set.seed(3)
  n <- 120
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  h <- 0.1

  x <- as.vector(scale(g1))
  d2 <- as.matrix(stats::dist(coords))^2
  K <- exp(-d2 / (2 * h^2)); diag(K) <- 0
  v <- as.vector(K %*% (x / rowSums(K)))
  expected <- sum(v^2) / (n - 1)^2

  res <- multiScaleAC(coords, g1, g2, h = h, perms = 0, verbose = FALSE)
  expect_equal(res@moment_var, expected, tolerance = 1e-10)
})

test_that("moments are computed without any permutations", {
  # The point of the exercise: the moments cost nothing but the kernel sweep.
  set.seed(4)
  n <- 100
  coords <- matrix(runif(2 * n), ncol = 2)
  res <- multiScaleAC(coords, rnorm(n), rnorm(n), h = c(0.05, 0.1), perms = 0, verbose = FALSE)
  expect_identical(ncol(res@permutations), 0L)
  expect_length(res@moment_var, 2L)
  expect_true(all(is.finite(res@moment_var)))
  expect_true(all(res@moment_var > 0))
  expect_true(methods::validObject(res))
})

test_that("moments can be declined", {
  set.seed(5)
  n <- 80
  coords <- matrix(runif(2 * n), ncol = 2)
  res <- multiScaleAC(coords, rnorm(n), rnorm(n), h = 0.1, perms = 10,
                      moments = FALSE, verbose = FALSE)
  expect_length(res@moment_var, 0L)
  expect_error(compute_significance(res, method = "moment"), regexp = "no analytical moments")
})

test_that("leverage falls as bandwidth grows", {
  # The diagnostic must behave the way the theory says, or it is not measuring what it claims. Small
  # bandwidths concentrate the statistic on few cells; large ones spread it.
  set.seed(6)
  n <- 300
  coords <- matrix(runif(2 * n), ncol = 2)
  res <- multiScaleAC(coords, rnorm(n), rnorm(n), h = c(0.02, 0.05, 0.1, 0.25),
                      perms = 0, verbose = FALSE)
  expect_true(all(diff(res@leverage) < 0))
  expect_true(all(res@leverage > 0 & res@leverage <= 1))
})

test_that("moment p-values agree with permutations when permutations can resolve that deep", {
  # The claim that matters. Not "fast" -- correct below the permutation floor.
  set.seed(7)
  n <- 250
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n)
  sm <- exp(-as.matrix(stats::dist(coords))^2 / (2 * 0.1^2))
  g2 <- as.vector(sm %*% g1) / rowSums(sm) + rnorm(n, sd = 1.5)
  hv <- c(0.08, 0.15)

  big <- multiScaleAC(coords, g1, g2, h = hv, perms = 50000, seed = 8, verbose = FALSE)
  p_perm <- compute_significance(big, method = "permutation")@p_values
  p_mom <- compute_significance(big, method = "moment")@p_values

  # Compared on the log scale, because agreement in the tail is what matters. The bound comes from
  # the measured calibration range of the default gamma method rather than being chosen: over
  # unflagged bandwidths it ran from 0.86 to 3.11 times the truth, so log10 differences up to about
  # 0.5 are expected behaviour, not a defect. See inst/validation/results/README.md.
  expect_lt(max(abs(log10(p_mom) - log10(p_perm))), 0.55)

  # Direction check: where they differ, the moment route must not be materially anti-conservative,
  # which is the failure mode that would matter.
  expect_true(all(p_mom > 0.4 * p_perm))
})

test_that("moment inference resolves below the permutation floor", {
  # With few permutations the smallest reportable p-value is 1/(perms+1). The moment route has no such
  # floor, which is the practical benefit.
  set.seed(9)
  n <- 250
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n)
  sm <- exp(-as.matrix(stats::dist(coords))^2 / (2 * 0.1^2))
  g2 <- as.vector(sm %*% g1) / rowSums(sm) + rnorm(n, sd = 0.2)   # strong signal
  res <- multiScaleAC(coords, g1, g2, h = 0.1, perms = 99, seed = 10, verbose = FALSE)

  p_perm <- compute_significance(res, method = "permutation")@p_values
  p_mom <- compute_significance(res, method = "moment")@p_values
  expect_equal(p_perm, 1 / 100)          # exactly at the floor: no permutation was as extreme
  expect_lt(p_mom, 1 / 100)              # the moment route goes deeper
})

test_that("auto mode uses moments where trusted and permutations where not", {
  set.seed(11)
  n <- 400
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rlnorm(n, 0, 1.2)
  g2 <- c(rep(0, n - 40), rlnorm(40, 1.5, 1))[sample.int(n)]
  # A wide bandwidth range so some bandwidths are flagged and some are not.
  res <- multiScaleAC(coords, g1, g2, h = c(0.01, 0.02, 0.1, 0.25), perms = 999,
                      seed = 12, verbose = FALSE)
  sig <- compute_significance(res, method = "auto")

  expect_length(sig@p_method, 4L)
  expect_true(all(sig@p_method %in% c("permutation", "moment:gamma")))
  # The provenance must match the diagnostic, not be arbitrary.
  flagged <- multiScaleAC:::.flag_untrustworthy_from_slots(res)
  expect_identical(sig@p_method == "permutation", flagged)
  expect_true(any(flagged))     # this configuration is meant to produce a mix
  expect_true(any(!flagged))
})

test_that("permute = flagged permutes only the flagged bandwidths", {
  set.seed(13)
  n <- 400
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rlnorm(n, 0, 1.2)
  g2 <- c(rep(0, n - 40), rlnorm(40, 1.5, 1))[sample.int(n)]
  hv <- c(0.01, 0.02, 0.1, 0.25)

  full <- multiScaleAC(coords, g1, g2, h = hv, perms = 499, seed = 14,
                       permute = "all", verbose = FALSE)
  part <- multiScaleAC(coords, g1, g2, h = hv, perms = 499, seed = 14,
                       permute = "flagged", verbose = FALSE)

  flagged <- multiScaleAC:::.flag_untrustworthy_from_slots(full)
  permuted <- apply(part@permutations, 1, function(r) any(is.finite(r)))
  expect_identical(permuted, flagged)
  expect_lt(sum(permuted), length(hv))     # something really was skipped

  # Where it did permute, the values must be identical to the full run: the same seed and the same
  # stream, just fewer rows contracted.
  expect_equal(part@permutations[flagged, , drop = FALSE],
               full@permutations[flagged, , drop = FALSE],
               tolerance = 1e-12)
  # And the observed curve is untouched by any of this.
  expect_identical(part@observed, full@observed)
})

test_that("every tail approximation is monotone and stays in (0, 1)", {
  # An Edgeworth expansion in particular can leave [0, 1] and is not monotone; the wrappers guard
  # both, and that guard is what is being tested.
  set.seed(15)
  n <- 200
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rlnorm(n, 0, 1)
  g2 <- c(rep(0, n - 30), rlnorm(30, 1.5, 1))[sample.int(n)]
  res <- multiScaleAC(coords, g1, g2, h = 0.08, perms = 0, verbose = FALSE)

  sd_i <- sqrt(res@moment_var)
  zs <- seq(-4, 6, by = 0.25)
  for (m in c("gamma", "edgeworth", "normal", "saddlepoint")) {
    p <- vapply(zs, function(z) {
      suppressWarnings(multiScaleAC:::.moment_tail(
        z * sd_i, res@moment_var, res@moment_skewness, res@moment_excess_kurtosis, method = m
      ))
    }, numeric(1))
    expect_true(all(p > 0 & p <= 1), info = m)
    # Non-increasing in z, allowing for the floor clamp.
    expect_true(all(diff(p) <= 1e-12), info = paste(m, "monotonicity"))
  }
})

test_that("two-sided moment p-values account for both tails of an asymmetric null", {
  set.seed(16)
  n <- 200
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rlnorm(n, 0, 1)
  g2 <- c(rep(0, n - 30), rlnorm(30, 1.5, 1))[sample.int(n)]
  res <- multiScaleAC(coords, g1, g2, h = c(0.05, 0.1), perms = 0, verbose = FALSE)
  two <- compute_significance(res, method = "moment", alternative = "two.sided")
  up <- compute_significance(res, method = "moment", alternative = "greater")
  lo <- compute_significance(res, method = "moment", alternative = "less")
  expect_true(all(two@p_values >= pmax(up@p_values, lo@p_values) - 1e-12))
  expect_true(all(two@p_values <= 1))
})

test_that("a permutation-free result still gets a null envelope", {
  set.seed(17)
  n <- 150
  coords <- matrix(runif(2 * n), ncol = 2)
  res <- multiScaleAC(coords, rnorm(n), rnorm(n), h = c(0.05, 0.1), perms = 0, verbose = FALSE)
  sig <- compute_significance(res, method = "moment")
  expect_true(all(is.finite(sig@null_lower)))
  expect_true(all(is.finite(sig@null_upper)))
  expect_true(all(sig@null_lower < sig@null_upper))
  # The envelope must bracket the exact null mean of 0.
  expect_true(all(sig@null_lower < 0 & sig@null_upper > 0))
})

test_that("the moment envelope is close to the permutation envelope", {
  set.seed(18)
  n <- 200
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n); g2 <- rnorm(n)
  hv <- c(0.05, 0.1, 0.2)
  withp <- compute_significance(
    multiScaleAC(coords, g1, g2, h = hv, perms = 20000, seed = 19, verbose = FALSE),
    method = "permutation")
  nop <- compute_significance(
    multiScaleAC(coords, g1, g2, h = hv, perms = 0, verbose = FALSE),
    method = "moment")
  expect_equal(nop@null_upper, withp@null_upper, tolerance = 0.1)
  expect_equal(nop@null_lower, withp@null_lower, tolerance = 0.1)
})
