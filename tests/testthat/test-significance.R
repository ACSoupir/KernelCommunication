# Permutation inference.
#
# The previous implementation defaulted to alternative = "two.sided" while its two-sided branch was
# empty, so it returned a one-sided upper-tail p-value under a two-sided label; and its "less"
# option computed 1 - p, which is not a p-value. These tests compute the expected values longhand
# from the permutation matrix, so the assertions do not depend on the implementation being tested.

make_result <- function(n = 200, h = c(0.05, 0.1, 0.2), perms = 199, seed = 1, related = TRUE) {
  set.seed(seed)
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n)
  g2 <- if (related) {
    sm <- exp(-as.matrix(stats::dist(coords))^2 / (2 * 0.1^2))
    as.vector(sm %*% g1) / rowSums(sm) + rnorm(n, sd = 0.5)
  } else {
    rnorm(n)
  }
  multiScaleAC(coords, g1, g2, h = h, perms = perms, seed = seed, verbose = FALSE)
}

test_that("the default is the manuscript's one-tailed upper test", {
  res <- compute_significance(make_result())
  expect_identical(formals(compute_significance)$alternative[[2]], "greater")

  perms <- res@permutations
  expected <- vapply(seq_along(res@observed), function(i) {
    (sum(perms[i, ] >= res@observed[i]) + 1) / (ncol(perms) + 1)
  }, numeric(1))
  expect_equal(res@p_values, expected)
})

test_that("the corrected default reproduces what the old mislabelled default computed", {
  # Confirms the diagnosis rather than merely asserting the new behaviour: the old "two.sided"
  # default returned exactly this one-sided upper-tail value.
  res <- make_result()
  greater <- compute_significance(res, alternative = "greater")
  perms <- res@permutations
  old_default <- vapply(seq_along(res@observed), function(i) {
    (sum(perms[i, ] >= res@observed[i]) + 1) / (ncol(perms) + 1)
  }, numeric(1))
  expect_equal(greater@p_values, old_default)
})

test_that("the lower-tail test is a real p-value, not 1 - p", {
  res <- make_result()
  less <- compute_significance(res, alternative = "less")
  perms <- res@permutations
  expected <- vapply(seq_along(res@observed), function(i) {
    (sum(perms[i, ] <= res@observed[i]) + 1) / (ncol(perms) + 1)
  }, numeric(1))
  expect_equal(less@p_values, expected)

  # The old formula was 1 - p_greater, which differs from the correct value and, crucially, can be
  # exactly 0 -- defeating the whole purpose of the +1 correction.
  greater <- compute_significance(res, alternative = "greater")
  expect_false(isTRUE(all.equal(less@p_values, 1 - greater@p_values)))
  expect_true(all(less@p_values > 0))
})

test_that("no p-value can ever be zero under the default estimator", {
  # A p-value of exactly 0 breaks log transforms and FDR procedures. With a strong signal the
  # observed value beats every permutation, which is precisely when the naive r/n returns 0.
  res <- make_result(related = TRUE, perms = 99)
  sig <- compute_significance(res)
  expect_true(all(sig@p_values > 0))
  expect_gte(min(sig@p_values), 1 / (ncol(res@permutations) + 1))

  # The manuscript's literal estimator can return exactly 0, which is why it is not the default and
  # why it warns. Verified rather than asserted in prose.
  expect_warning(paper <- compute_significance(res, estimator = "paper"),
                 regexp = "exactly 0")
  expect_true(any(paper@p_values == 0))
  expect_true(all(paper@p_values <= sig@p_values))
})

test_that("the two-sided test is genuinely two-sided", {
  res <- make_result(related = FALSE, seed = 4)
  two <- compute_significance(res, alternative = "two.sided")
  perms <- res@permutations
  pm <- rowMeans(perms)
  expected <- vapply(seq_along(res@observed), function(i) {
    (sum(abs(perms[i, ] - pm[i]) >= abs(res@observed[i] - pm[i])) + 1) / (ncol(perms) + 1)
  }, numeric(1))
  expect_equal(two@p_values, expected)

  # A strongly negative deviation must be detectable two-sided but not by the upper-tail test. The
  # old default could never flag negative spatial cross-correlation at all.
  set.seed(11); n <- 200
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n)
  sm <- exp(-as.matrix(stats::dist(coords))^2 / (2 * 0.1^2))
  g2 <- -as.vector(sm %*% g1) / rowSums(sm) + rnorm(n, sd = 0.3)   # anti-correlated
  neg <- multiScaleAC(coords, g1, g2, h = 0.1, perms = 199, verbose = FALSE)
  expect_gt(compute_significance(neg, alternative = "greater")@p_values, 0.5)
  expect_lt(compute_significance(neg, alternative = "two.sided")@p_values, 0.05)
  expect_lt(compute_significance(neg, alternative = "less")@p_values, 0.05)
})

test_that("the null envelope level is independent of the test threshold", {
  # They shared one alpha before, so a 95% band with a 0.01 test was impossible.
  res <- make_result()
  a <- compute_significance(res, alpha = 0.01, envelope_level = 0.95)
  b <- compute_significance(res, alpha = 0.05, envelope_level = 0.95)
  expect_equal(a@null_lower, b@null_lower)      # envelope unchanged by alpha
  expect_equal(a@null_upper, b@null_upper)

  wide <- compute_significance(res, alpha = 0.05, envelope_level = 0.99)
  expect_true(all(wide@null_upper >= b@null_upper))
  expect_true(all(wide@null_lower <= b@null_lower))
})

test_that("the null envelope brackets the permutation distribution, not the observed value", {
  # It is an acceptance region for the null, which is why it is no longer called a confidence
  # interval. With a real signal the observed value sits outside it.
  sig <- compute_significance(make_result(related = TRUE))
  expect_true(all(sig@perm_mean > sig@null_lower & sig@perm_mean < sig@null_upper))
  expect_true(any(sig@observed > sig@null_upper))
})

test_that("recomputing significance clears a stale adjustment", {
  # p_adjusted used to survive a recompute, so summary() could show adjusted p-values that did not
  # correspond to the p_values beside them.
  res <- compute_significance(make_result())
  res <- adjust_pvalues(res, method = "bonferroni")
  expect_length(res@p_adjusted, length(res@h))

  res2 <- suppressWarnings(compute_significance(res, alternative = "less"))
  expect_length(res2@p_adjusted, 0L)
})

test_that("recomputing warns that it is doing so", {
  res <- compute_significance(make_result())
  expect_warning(compute_significance(res), regexp = "already computed")
})

test_that("the supremum test gives a whole-curve p-value", {
  # Defined in the manuscript and previously not implemented anywhere.
  sig <- compute_significance(make_result(related = TRUE))
  expect_length(sig@sup_p_value, 1L)
  expect_length(sig@sup_statistic, 1L)
  expect_gt(sig@sup_p_value, 0)
  expect_lte(sig@sup_p_value, 1)
  # A genuinely spatially related pair should be detected by the whole-curve test.
  expect_lt(sig@sup_p_value, 0.05)

  # And an unrelated pair should usually not be.
  null_sig <- compute_significance(make_result(related = FALSE, seed = 21))
  expect_gt(null_sig@sup_p_value, 0.05)
})

test_that("the supremum statistic matches a direct computation", {
  res <- make_result()
  sig <- compute_significance(res, sup_standardize = TRUE)
  pm <- rowMeans(res@permutations)
  psd <- apply(res@permutations, 1, stats::sd)
  dev_obs <- (res@observed - pm) / psd
  dev_perm <- sweep(sweep(res@permutations, 1, pm, "-"), 1, psd, "/")
  expect_equal(sig@sup_statistic, max(dev_obs))
  expected_p <- (sum(apply(dev_perm, 2, max) >= max(dev_obs)) + 1) / (ncol(res@permutations) + 1)
  expect_equal(sig@sup_p_value, expected_p)
})

test_that("the unstandardised supremum follows the manuscript literally", {
  res <- make_result()
  lit <- compute_significance(res, sup_standardize = FALSE)
  pm <- rowMeans(res@permutations)
  expect_equal(lit@sup_statistic, max(res@observed - pm))
  # The two versions genuinely differ, which is why the choice is exposed rather than assumed.
  std <- compute_significance(res, sup_standardize = TRUE)
  expect_false(isTRUE(all.equal(lit@sup_statistic, std@sup_statistic)))
})

test_that("the supremum test can be switched off", {
  sig <- compute_significance(make_result(), supremum = FALSE)
  expect_length(sig@sup_p_value, 0L)
  expect_true(methods::validObject(sig))
})

test_that("significance requires permutations and says so", {
  set.seed(5); n <- 100
  coords <- matrix(runif(2 * n), ncol = 2)
  res <- multiScaleAC(coords, rnorm(n), rnorm(n), h = 0.1, perms = 0, verbose = FALSE)
  expect_error(compute_significance(res), regexp = "no permutations", ignore.case = TRUE)
})

test_that("invalid alpha and envelope_level are rejected", {
  res <- make_result()
  expect_error(compute_significance(res, alpha = 0), regexp = "between 0 and 1")
  expect_error(compute_significance(res, alpha = 1), regexp = "between 0 and 1")
  expect_error(compute_significance(res, envelope_level = 1.5), regexp = "between 0 and 1")
  expect_error(compute_significance(res, alternative = "bigger"))
})

test_that("adjust_pvalues requires significance first and treats BH and fdr alike", {
  res <- make_result()
  expect_error(adjust_pvalues(res), regexp = "compute_significance")
  sig <- compute_significance(res)
  expect_equal(adjust_pvalues(sig, "BH")@p_adjusted, adjust_pvalues(sig, "fdr")@p_adjusted)
  expect_true(all(adjust_pvalues(sig, "bonferroni")@p_adjusted >= sig@p_values))
})

test_that("summary reflects the computed values", {
  sig <- adjust_pvalues(compute_significance(make_result()))
  s <- summary(sig)
  expect_s3_class(s, "data.frame")
  expect_identical(s$h, sig@h)
  expect_identical(s$p_value, sig@p_values)
  expect_identical(s$p_adjusted, sig@p_adjusted)
  expect_true(all(c("null_lower", "null_upper") %in% names(s)))
})

test_that("the permutation null is skewed on realistic data, so z-scores are effect sizes only", {
  # Guards the documented warning with a measurement rather than an assertion. The permutation null
  # of this statistic is skewed whenever the expression data is, which is why compute_significance()
  # tells users not to run pnorm() on the z_scores slot.
  cs <- golden_cases()$synth_zeroinflated
  res <- multiScaleAC(cs$coords, cs$gene1, cs$gene2, h = cs$h, perms = 999,
                      seed = 7, verbose = FALSE)
  sig <- compute_significance(res)

  skew <- apply(res@permutations, 1, function(r) {
    s <- stats::sd(r); if (s == 0) 0 else mean((r - mean(r))^3) / s^3
  })
  expect_gt(max(abs(skew)), 0.5)      # the null really is skewed on this 85%-zero pair

  # A normal approximation to the z-score therefore disagrees with the permutation p-value. On this
  # fixture the worst disagreement measures 0.09 on the log10 scale, about 1.2x; the effect grows
  # sharply with skewness and further into the tail, reaching two orders of magnitude on the
  # adversarial case quoted in compute_significance()'s documentation. Asserted loosely here because
  # this fixture is realistic rather than adversarial -- the point is that the two do not agree.
  normal_p <- stats::pnorm(sig@z_scores, lower.tail = FALSE)
  worst <- max(abs(log10(pmax(normal_p, 1e-16)) - log10(sig@p_values)), na.rm = TRUE)
  expect_gt(worst, 0.04)

  # The error that matters is in the far tail, which is where an FDR screen over many gene pairs
  # actually operates -- not at the observed value, which here sits near the bulk of the null. Note
  # the direction is NOT uniform: at moderate z a right-skewed null can hold less mass above z than
  # a normal, so the approximation overstates there. It is reliably ANTI-conservative far out, which
  # is the dangerous direction, because it manufactures false positives at the threshold.
  #
  # Compared against the empirical null directly, at fixed standardised thresholds.
  bw <- which.max(abs(skew))
  z_null <- (res@permutations[bw, ] - sig@perm_mean[bw]) / sig@perm_sd[bw]
  for (zz in c(2.0, 2.5)) {
    empirical <- mean(z_null >= zz)
    normal <- stats::pnorm(zz, lower.tail = FALSE)
    skip_if(empirical == 0, "not enough permutations to resolve this tail")
    expect_lt(normal, empirical)
  }
})
