# Exact moments of the permutation null distribution.
#
# The manuscript's Limitations section names this as future work: "possible identification of
# analytical mean and variance of the distributions, eliminating the requirement of permutations."
# The first four moments, and the full cross-bandwidth covariance, are available in closed form.
#
# WHY. Write the statistic with the row-normalised weight matrix W. Because v_h = K (x / rowsum) and
# K is symmetric, v_h = t(W) x, so
#
#     I(h) = <v_h, y> / S,     S = sum(x^2)
#
# Under the permutation scheme used here only y is permuted. So with a = v_h - mean(v_h) and b = y
# (already centred),
#
#     I_pi(h) = U / S,     U = sum_i a_i b_pi(i)
#
# which is a *linear permutation statistic*: a fixed vector paired with a uniformly random
# permutation of another fixed vector. Its moments over the n! permutations are classical, and depend
# on the data only through power sums.
#
# Note mean(v_h) = 0 exactly whenever W has rows summing to 1, because W %*% 1 = 1 and sum(x) = 0.
# It is nevertheless computed rather than assumed, because an isolated cell gives that row of W a sum
# of zero and breaks the identity.
#
# COST. v_h is already computed to obtain the observed value. The power sums are then O(N) per
# bandwidth, and the power sums of y are bandwidth-independent. That replaces the
# O(n_h x N x perms) permutation product and eliminates the N x perms permuted-expression matrix
# entirely, with eps = 0 preserved exactly and no thresholding or gene filtering introduced.
#
# VERIFICATION. The formulas below were checked against 2,000,000 explicit permutations in a
# strongly skewed regime (skewness 2.37, excess kurtosis 7.25) before being written down, agreeing to
# 1.002, 1.006 and 1.008 on the second, third and fourth moments and 0.999 on the cross-bandwidth
# covariance. test-moments.R re-checks them on every run.

#' Exact moments of the permutation distribution
#'
#' Internal. Computes the first four moments of the permutation null of
#' \eqn{I(h)}, for every bandwidth, plus diagnostics.
#'
#' @param V n_h x N matrix of spatially lagged gene1, from the kernel pass.
#' @param gene2_scaled Length-N centred and scaled gene2.
#' @param denom \code{sum(gene1_scaled^2)}.
#'
#' @return A list with, per bandwidth, \code{mean}, \code{var}, \code{skewness},
#'   \code{excess_kurtosis}, and the diagnostics \code{leverage} (largest single
#'   contribution to the variance) and \code{n_effective}.
#' @keywords internal
#' @noRd
.permutation_moments <- function(V, gene2_scaled, denom) {
  n_h <- nrow(V)
  n <- ncol(V)
  if (n < 4) {
    # The fourth moment formula divides by (n-3), so it is undefined below four observations. The
    # statistic itself is still computable, so this returns empty moments rather than failing: a
    # three-cell sample should not be rejected outright just because its kurtosis is undefined.
    return(NULL)
  }

  # Centre each bandwidth's lagged vector. Exactly zero in theory; computed for robustness.
  a <- V - rowMeans(V)
  b <- gene2_scaled - mean(gene2_scaled)

  A2 <- rowSums(a^2); A3 <- rowSums(a^3); A4 <- rowSums(a^4)
  B2 <- sum(b^2);     B3 <- sum(b^3);     B4 <- sum(b^4)

  # Moments of U = sum_i a_i b_pi(i) over a uniformly random permutation.
  m2 <- A2 * B2 / (n - 1)
  m3 <- n * A3 * B3 / ((n - 1) * (n - 2))

  P <- A2^2; Q <- A4; p <- B2^2; q <- B4
  m4 <- Q * q * (1 / n + 4 / (n * (n - 1))) +
        3 * (P - Q) * (p - q) / (n * (n - 1)) +
        6 * (P - 2 * Q) * (p - 2 * q) / (n * (n - 1) * (n - 2)) +
        9 * (P - 2 * Q) * (p - 2 * q) / (n * (n - 1) * (n - 2) * (n - 3))

  # Standardised shape is invariant to the S scaling, so it is computed on U.
  ok <- is.finite(m2) & m2 > 0
  skewness <- rep(NA_real_, n_h)
  exkurt <- rep(NA_real_, n_h)
  skewness[ok] <- m3[ok] / m2[ok]^1.5
  exkurt[ok] <- m4[ok] / m2[ok]^2 - 3

  # Diagnostics for whether a normal-ish approximation can be trusted at this bandwidth. The
  # combinatorial central limit theorem for this statistic requires no single term to dominate;
  # leverage is the largest single contribution to the variance, so a value near 1 means the
  # statistic is effectively driven by one cell and no continuous approximation will hold.
  leverage <- rep(NA_real_, n_h)
  leverage[ok] <- apply(a[ok, , drop = FALSE]^2, 1, max) / A2[ok]

  list(
    mean            = rep(0, n_h),
    var             = m2 / denom^2,
    skewness        = skewness,
    excess_kurtosis = exkurt,
    leverage        = leverage,
    n_effective     = ifelse(ok, 1 / leverage, NA_real_),
    # Raw U-scale moments, kept for the approximations and for the covariance.
    .m2 = m2, .m3 = m3, .m4 = m4, .a = a, .B2 = B2, .n = n, .denom = denom
  )
}

#' Exact cross-bandwidth covariance of the permutation null
#'
#' Internal. \eqn{Cov(I(h), I(h')) = \langle a_h, a_{h'}\rangle B_2 / ((n-1) S^2)}.
#' Makes a whole-curve test possible without permutations, since the joint null
#' of the curve is characterised by this matrix plus the marginal moments.
#'
#' @param mom Output of \code{.permutation_moments()}.
#' @return n_h x n_h covariance matrix on the scale of \eqn{I}.
#' @keywords internal
#' @noRd
.permutation_covariance <- function(mom) {
  a <- mom$.a
  tcrossprod(a) * mom$.B2 / ((mom$.n - 1) * mom$.denom^2)
}

#' Should moment-based inference be trusted at this bandwidth?
#'
#' Internal. Flags bandwidths where a continuous approximation to the discrete
#' permutation null is unreliable, so they can fall back to permutations.
#'
#' Both criteria matter and neither subsumes the other. High leverage means one
#' cell dominates the statistic, so no smooth approximation applies. Extreme
#' kurtosis means the null is far from any of the shapes the approximations can
#' represent. This is worst at small bandwidths, which is exactly where the
#' manuscript reports the curve is unstable, so the flagging tends to coincide
#' with the region users should already treat cautiously.
#'
#' @param mom Output of \code{.permutation_moments()}.
#' @param max_leverage Leverage above which a bandwidth is flagged.
#' @param max_kurtosis Excess kurtosis above which a bandwidth is flagged.
#'
#' @return Logical vector, TRUE where moment-based inference is not trusted.
#' @keywords internal
#' @noRd
#' @details
#' The default thresholds come from the calibration run in
#' \code{inst/validation/calibrate_moments.R}, not from theory. A leverage
#' threshold of 0.05 flagged a configuration where every approximation was
#' accurate to within 60\%, so it was relaxed to 0.10, which separates the
#' configurations where the recommended approximation stayed within 6\% of the
#' truth from those where it did not.
.flag_untrustworthy <- function(mom, max_leverage = 0.10, max_kurtosis = 20) {
  bad <- !is.finite(mom$var) | mom$var <= 0 |
         !is.finite(mom$skewness) | !is.finite(mom$excess_kurtosis)
  bad | mom$leverage > max_leverage | abs(mom$excess_kurtosis) > max_kurtosis
}
