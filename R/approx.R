# Turning exact moments into tail probabilities.
#
# The moments in R/moments.R are exact. Mapping them to a p-value is not, and is the genuinely open
# part of this work. Four methods are provided and the default is chosen by measurement, not
# assertion: inst/validation/calibrate_moments.R compares each against large-permutation ground
# truth and writes a table that inst/validation/results/ records.
#
# RESULT OF THAT MEASUREMENT (200,000 permutations, four configurations from benign to adversarial,
# scored at target upper-tail probabilities down to 1e-4 because that is where an FDR screen over
# ~1100 gene pairs operates). Ratios are approximate / empirical; below 1 is anti-conservative, the
# direction that manufactures false positives. Over the bandwidths that would actually use moment
# inference, i.e. those the diagnostics do not flag:
#
#   method        median   range           anti-conservative (< 0.5)
#   gamma          1.05    [0.86, 3.11]    0 of 24
#   edgeworth      1.01    [0.78, 3.11]    0 of 24
#   saddlepoint    1.00    [0.28, 1.53]    2 of 24
#   normal         1.01    [0.12, 5.47]    2 of 24
#
# THE DEFAULT IS "gamma", a three-moment shifted gamma. It is the only method with no
# anti-conservative failures and a bounded worst case; where it errs it errs conservatively.
#
# Two findings worth recording, because both contradict what theory would suggest:
#
# 1. The saddlepoint is the best-founded method here and measured WORSE than the gamma. On the real
#    ccRCCExample data at h = 0.08 it returned ratios of 0.59, 0.38 and 0.28 at targets of 1e-2, 1e-3
#    and 1e-4, understating the tail. It is retained because it is excellent on the strongly skewed
#    configurations (0.91 to 1.02 where the gamma ranged 0.83 to 1.20), but it is not the default, and
#    the discrepancy on near-symmetric continuous data is unexplained -- possibly the 128-atom
#    representation of a continuous gene, possibly the conditioning. Treat it as experimental.
# 2. The normal approximation fails on REAL DATA, not merely on contrived examples. On ccRCCExample at
#    h = 0.02 (exact skewness 1.91) it returned 0.01 at a target of 1e-2 and 0.00 at 1e-3 and 1e-4.
#    Anyone converting the z_scores slot with pnorm() on data like this would generate false
#    positives at exactly the threshold they were screening at.
#
# Edgeworth is included for comparison, not recommended: it is not monotone and it collapses in the
# deep tail under skew, returning 32.52 then 0.00 on adjacent targets of one configuration.
#
# Full table: inst/validation/results/calibration.csv, regenerate with
# Rscript inst/validation/calibrate_moments.R

#' Upper-tail probability from moments
#'
#' Internal. All methods return \eqn{P(I \ge q)} under the permutation null.
#'
#' @param q Observed value(s), on the scale of \eqn{I}, one per bandwidth.
#' @param var,skewness,excess_kurtosis Exact moments on the scale of \eqn{I}, as
#'   stored on the result object.
#' @param method One of "gamma", "saddlepoint", "edgeworth", "normal".
#' @param lower Return the lower tail instead.
#' @param n_atoms Atoms used to represent the distribution of gene2 in the
#'   saddlepoint. Ignored by the other methods.
#' @param a n_h x N centred lagged matrix, and \code{b} the centred scaled gene2,
#'   with \code{denom}. Required only by the saddlepoint, which cannot be
#'   evaluated from summary moments alone. The other three methods are
#'   scale-free functions of the standardised deviate and the shape, so they work
#'   from the stored slots without the original data.
#' @param b,denom See \code{a}.
#'
#' @return Numeric vector of tail probabilities, clamped to (0, 1).
#' @keywords internal
#' @noRd
.moment_tail <- function(q, var, skewness, excess_kurtosis,
                         method = c("gamma", "saddlepoint", "edgeworth", "normal"),
                         lower = FALSE, n_atoms = 128L,
                         a = NULL, b = NULL, denom = NULL) {
  method <- match.arg(method)
  n_h <- length(var)

  sd_i <- sqrt(var)
  skew <- skewness
  exk <- excess_kurtosis
  qq <- q

  # Lower tail by reflection, so each method only implements one direction.
  if (lower) {
    qq <- -qq
    skew <- -skew
  }

  if (method == "saddlepoint" && (is.null(a) || is.null(b) || is.null(denom))) {
    warning("the saddlepoint needs the lagged matrix and gene2, which are not stored on the result ",
            "object; falling back to the gamma approximation", call. = FALSE)
    method <- "gamma"
  }

  out <- rep(NA_real_, n_h)
  for (i in seq_len(n_h)) {
    if (!is.finite(sd_i[i]) || sd_i[i] <= 0 || !is.finite(qq[i])) next
    z <- qq[i] / sd_i[i]
    out[i] <- switch(
      method,
      normal     = stats::pnorm(z, lower.tail = FALSE),
      gamma      = .tail_gamma(z, skew[i]),
      edgeworth  = .tail_edgeworth(z, skew[i], exk[i]),
      saddlepoint = .tail_saddlepoint(qq[i] * denom, if (lower) -a[i, ] else a[i, ], b,
                                      n_atoms = n_atoms,
                                      fallback_z = z, fallback_skew = skew[i])
    )
  }

  # Every method is clamped away from 0 and 1. A p-value of exactly 0 is not a usable estimate, and
  # the practical floor on resolution is no better than permutations would give.
  pmin(pmax(out, .Machine$double.eps), 1)
}

#' Three-moment shifted gamma tail
#'
#' Matches mean, variance and skewness with a shifted gamma (Pearson type III).
#' Monotone in \code{z} by construction and confined to [0, 1], unlike an
#' Edgeworth expansion.
#'
#' @param z Standardised deviate.
#' @param skew Exact skewness of the null.
#' @return Upper-tail probability.
#' @keywords internal
#' @noRd
.tail_gamma <- function(z, skew) {
  if (!is.finite(skew) || abs(skew) < 1e-8) {
    return(stats::pnorm(z, lower.tail = FALSE))
  }
  # Standardised shifted gamma: shape k = 4/skew^2, and the standardised variable is
  # (2/skew) * (gamma(k, scale = 1) / sqrt(k) - sqrt(k)) for positive skew.
  k <- 4 / skew^2
  if (!is.finite(k) || k <= 0) return(stats::pnorm(z, lower.tail = FALSE))

  # Map z back to the gamma scale. For positive skew the support is bounded below at -2/skew.
  if (skew > 0) {
    x <- k + z * sqrt(k)
    if (x <= 0) return(1)                      # below the support: all mass is above
    stats::pgamma(x, shape = k, lower.tail = FALSE)
  } else {
    x <- k - z * sqrt(k)
    if (x <= 0) return(.Machine$double.eps)    # above the support
    stats::pgamma(x, shape = k, lower.tail = TRUE)
  }
}

#' Two-term Edgeworth tail
#'
#' Included for comparison and because it is the textbook first move, not because
#' it is recommended. An Edgeworth expansion is not monotone and can leave
#' [0, 1]; both are guarded here, but guarding a non-monotone expansion does not
#' make it accurate, and the measurements in this file's header show it
#' oscillating between 0.6x and 1.6x of the truth.
#'
#' @param z Standardised deviate.
#' @param skew,exk Exact skewness and excess kurtosis.
#' @return Upper-tail probability.
#' @keywords internal
#' @noRd
.tail_edgeworth <- function(z, skew, exk) {
  if (!is.finite(skew)) skew <- 0
  if (!is.finite(exk)) exk <- 0
  p <- stats::pnorm(z, lower.tail = FALSE) +
       stats::dnorm(z) * (skew / 6 * (z^2 - 1) +
                          exk / 24 * (z^3 - 3 * z) +
                          skew^2 / 72 * (z^5 - 10 * z^3 + 15 * z))
  # Fall back rather than return an impossible value, which is a real failure mode of this expansion
  # in the far tail where it is most needed.
  if (!is.finite(p) || p <= 0 || p >= 1) {
    return(.tail_gamma(z, skew))
  }
  p
}

#' Double saddlepoint tail for a linear permutation statistic
#'
#' Follows the conditional construction of Robinson (1980) and Booth and Butler
#' (1990): represent the permutation distribution of \eqn{\sum a_i b_{\pi(i)}} as
#' the distribution of \eqn{\sum a_i Y_i} with \eqn{Y_i} drawn independently from
#' the empirical distribution of \eqn{b}, conditioned on \eqn{\sum Y_i =
#' \sum b_i}. A double saddlepoint in the two corresponding tilting parameters
#' then gives the tail.
#'
#' Every cumulant evaluation goes through a log-sum-exp, because the naive form
#' overflows: a direct implementation of this failed at moderate \eqn{z} on the
#' adversarial case used to size these methods.
#'
#' The distribution of \eqn{b} is represented by at most \code{n_atoms} weighted
#' atoms. Distinct values are used when there are few of them, which is the common
#' case for zero-inflated expression, and quantile bins otherwise. Each cumulant
#' evaluation is then O(N x atoms) rather than O(N^2).
#'
#' @param u Observed value on the U scale.
#' @param a Centred lagged vector for this bandwidth.
#' @param b Centred and scaled gene2.
#' @param n_atoms Maximum atoms representing the distribution of \code{b}.
#' @param fallback_z,fallback_skew Used if the saddlepoint fails to converge.
#'
#' @return Upper-tail probability, or the gamma approximation if the solve fails.
#' @keywords internal
#' @noRd
.tail_saddlepoint <- function(u, a, b, n_atoms = 128L,
                              fallback_z = NA_real_, fallback_skew = NA_real_) {
  if (is.null(b)) return(.tail_gamma(fallback_z, fallback_skew))

  atoms <- .atomize(b, n_atoms)
  bv <- atoms$value
  lw <- log(atoms$weight)

  # log-sum-exp over the atoms for each i, given the tilt c_i.
  # K(c) = sum_i log( sum_j w_j exp(c_i b_j) )
  cgf_parts <- function(cvec) {
    M <- outer(cvec, bv) + rep(lw, each = length(cvec))
    mx <- apply(M, 1, max)
    lse <- mx + log(rowSums(exp(M - mx)))
    list(lse = lse, M = M, mx = mx)
  }
  # First two moments of each tilted atom distribution.
  tilted <- function(cvec) {
    pr <- cgf_parts(cvec)
    Wt <- exp(pr$M - pr$lse)                 # rows sum to 1
    m1 <- as.vector(Wt %*% bv)
    m2 <- as.vector(Wt %*% (bv^2))
    list(lse = pr$lse, m1 = m1, v = pmax(m2 - m1^2, 0))
  }

  n <- length(a)
  # Solve for (lam, tau) with grad K = (u, 0), and for tau0 with dK/dtau(0, tau) = 0.
  solve2 <- function(target_u) {
    th <- c(0, 0)
    for (it in 1:100) {
      tl <- tilted(th[1] * a + th[2])
      g <- c(sum(a * tl$m1) - target_u, sum(tl$m1))
      H <- matrix(c(sum(a^2 * tl$v), sum(a * tl$v),
                    sum(a * tl$v),   sum(tl$v)), 2, 2)
      if (!all(is.finite(g)) || !all(is.finite(H))) return(NULL)
      step <- tryCatch(solve(H, g), error = function(e) NULL)
      if (is.null(step)) return(NULL)
      # Damp the step so a poor starting point cannot diverge.
      sc <- min(1, 4 / max(abs(step)))
      th <- th - sc * step
      if (max(abs(g)) < 1e-11 * max(1, abs(target_u))) {
        return(list(theta = th, lse = tl$lse, H = H))
      }
    }
    NULL
  }
  solve_null <- function() {
    tau <- 0
    for (it in 1:100) {
      tl <- tilted(rep(tau, n))
      g <- sum(tl$m1)
      H <- sum(tl$v)
      if (!is.finite(g) || !is.finite(H) || H <= 0) return(NULL)
      step <- g / H
      tau <- tau - min(1, 4 / max(abs(step), 1e-12)) * step
      if (abs(g) < 1e-12) return(list(tau = tau, lse = tl$lse, H = H))
    }
    NULL
  }

  fit <- solve2(u)
  nul <- solve_null()
  if (is.null(fit) || is.null(nul)) return(.tail_gamma(fallback_z, fallback_skew))

  lam <- fit$theta[1]; tau <- fit$theta[2]
  K_hat <- sum(fit$lse)
  K_null <- sum(nul$lse)

  # w = sign(lam) sqrt(2 [ (K(0,tau0) - tau0 t0) - (K(lam,tau) - lam u - tau t0) ]) with t0 = 0.
  inner <- 2 * (K_null - (K_hat - lam * u))
  if (!is.finite(inner) || inner < 0) return(.tail_gamma(fallback_z, fallback_skew))
  w <- sign(lam) * sqrt(inner)

  detH <- fit$H[1, 1] * fit$H[2, 2] - fit$H[1, 2] * fit$H[2, 1]
  if (!is.finite(detH) || detH <= 0 || nul$H <= 0) {
    return(.tail_gamma(fallback_z, fallback_skew))
  }
  u_hat <- lam * sqrt(detH / nul$H)

  # Near lam = 0 the Lugannani-Rice correction is singular; the distribution is essentially at its
  # centre there, so the normal term alone is the right answer.
  if (abs(w) < 1e-6 || abs(u_hat) < 1e-10) return(0.5)

  p <- stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * (1 / u_hat - 1 / w)
  if (!is.finite(p) || p <= 0 || p >= 1) return(.tail_gamma(fallback_z, fallback_skew))
  p
}

#' Represent a vector by weighted atoms
#'
#' Internal. Uses distinct values when there are few (common for zero-inflated
#' expression, where most cells share the value zero) and quantile bins
#' otherwise, so the saddlepoint cumulant evaluation is O(N x atoms).
#'
#' @param b Numeric vector.
#' @param n_atoms Maximum number of atoms.
#' @return List with \code{value} and \code{weight}, weights summing to 1.
#' @keywords internal
#' @noRd
.atomize <- function(b, n_atoms = 128L) {
  ub <- unique(b)
  if (length(ub) <= n_atoms) {
    tb <- table(b)
    return(list(value = as.numeric(names(tb)), weight = as.numeric(tb) / length(b)))
  }
  # Equal-count bins, represented by their means so the first moment is preserved exactly.
  cuts <- stats::quantile(b, probs = seq(0, 1, length.out = n_atoms + 1L), names = FALSE)
  cuts[1] <- -Inf; cuts[length(cuts)] <- Inf
  grp <- cut(b, breaks = unique(cuts), include.lowest = TRUE, labels = FALSE)
  val <- tapply(b, grp, mean)
  wt <- tapply(b, grp, length) / length(b)
  list(value = as.numeric(val), weight = as.numeric(wt))
}
