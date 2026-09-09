#' Compute Significance from Permutation Results
#'
#' Adds permutation p-values, a null envelope, and a whole-curve test to a
#' \code{\linkS4class{kernelcomm_result}}.
#'
#' @section Tail direction:
#' The default is \code{"greater"}, a one-tailed upper test, matching the
#' companion manuscript: "counting the number of permutations greater than the
#' observed (one tailed) and dividing by the total number of permutations".
#'
#' This is a corrected default rather than a changed one. The previous version
#' defaulted to \code{"two.sided"} but its two-sided branch was empty, so it
#' returned a one-sided upper-tail p-value under a two-sided label. Its
#' \code{"less"} option computed \code{1 - p}, which is not a p-value: it can be
#' exactly 0, and it double-counts ties. Both are fixed, and \code{"two.sided"}
#' now genuinely tests both tails.
#'
#' @section Estimator:
#' \code{"add_one"} (default) computes \eqn{(r + 1) / (n + 1)}, where \eqn{r}
#' counts permutations at least as extreme as the observed value. \code{"paper"}
#' computes \eqn{r / n}, matching the manuscript literally. \code{"add_one"} is
#' the default because \eqn{r / n} can return exactly 0, which is not a valid
#' p-value estimate and breaks any log transform or FDR procedure downstream; the
#' two differ negligibly except in the far tail, which is exactly where a zero
#' would matter.
#'
#' @section Null envelope, not a confidence interval:
#' \code{null_lower} and \code{null_upper} are quantiles of the \emph{permutation
#' distribution}. They describe where \eqn{I(h)} would fall under complete spatial
#' randomness, so the observed curve leaving the band is evidence against the
#' null. They are not a confidence interval for the observed statistic, and were
#' previously named \code{ci_lower}/\code{ci_upper} and plotted as though they
#' were. The envelope level is now separate from the test threshold, which shared
#' one \code{alpha} before, making a 95\% band with a 0.01 test impossible.
#'
#' @section Why z-scores should not be turned into p-values:
#' \code{z_scores} are provided as a comparable effect size, standardising the
#' observed value by the permutation mean and standard deviation. Do not convert
#' them with \code{pnorm()}. The permutation null of this statistic is skewed
#' whenever the expression data is skewed, which for zero-inflated single-cell
#' data is the normal case. On a deliberately adversarial example (300 cells, 85\%
#' zeros, null skewness 2.4) a normal approximation returned an upper-tail
#' probability 100 times too small at the 0.0001 level, which would manufacture
#' false positives at exactly the threshold an FDR screen operates on. Use the
#' permutation p-values, or the moment-based methods, which model the skewness.
#'
#' @section Permutation or moment based:
#' \code{method = "permutation"} (default) uses the permutation distribution.
#' \code{"moment"} uses the exact analytical moments of the permutation null,
#' which \code{\link{multiScaleAC}} computes for free, mapped to a tail
#' probability by \code{moment_dist}. \code{"auto"} uses moments where the
#' diagnostics say a continuous approximation is trustworthy and permutations
#' elsewhere, recording the choice per bandwidth in \code{p_method}.
#'
#' The moments are exact. Mapping them to a tail probability is an approximation,
#' and the default mapping was selected by measurement against 200,000
#' permutations rather than by theory; see
#' \code{inst/validation/calibrate_moments.R}. Permutations remain the default
#' because they need no such approximation.
#'
#' @param result Object of class \code{\linkS4class{kernelcomm_result}}
#' @param method Character. \code{"permutation"}, \code{"moment"}, or
#'   \code{"auto"}.
#' @param moment_dist Character, how moments become a tail probability.
#'   \code{"gamma"} (default) matches three moments with a shifted gamma and was
#'   the only method with no anti-conservative failures in calibration.
#'   \code{"saddlepoint"} is better founded and excellent under strong skew but
#'   measured worse on near-symmetric real data; treat it as experimental.
#'   \code{"edgeworth"} and \code{"normal"} are provided for comparison and are
#'   not recommended: both fail badly in the deep tail under skew.
#' @param alternative Character. \code{"greater"} (default, the manuscript's
#'   one-tailed upper test), \code{"less"}, or \code{"two.sided"}.
#' @param alpha Significance threshold for the \code{significant} slot.
#' @param envelope_level Coverage of the null envelope, e.g. 0.95 for the 2.5th
#'   to 97.5th percentiles of the permutation distribution. Independent of
#'   \code{alpha}.
#' @param estimator Character. \code{"add_one"} for \eqn{(r+1)/(n+1)} or
#'   \code{"paper"} for \eqn{r/n}.
#' @param supremum Logical. Also compute the whole-curve supremum test.
#' @param sup_standardize Logical. When \code{TRUE} (default), the supremum is
#'   taken over deviations standardised by the permutation standard deviation at
#'   each bandwidth, so every bandwidth contributes comparably. When
#'   \code{FALSE}, raw deviations are used, which follows the manuscript
#'   literally but is dominated by whichever bandwidth has the largest scale.
#'
#' @return The same object with these slots populated: \code{p_values},
#'   \code{z_scores}, \code{perm_mean}, \code{perm_sd}, \code{null_lower},
#'   \code{null_upper}, \code{significant}, and when \code{supremum = TRUE},
#'   \code{sup_statistic} and \code{sup_p_value}.
#'
#' @seealso \code{\link{adjust_pvalues}} for multiple-testing correction.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' coords <- matrix(runif(2 * n), ncol = 2)
#' g1 <- rnorm(n)
#' # gene2 spatially related to gene1, so there is something to detect
#' sm <- exp(-as.matrix(dist(coords))^2 / (2 * 0.1^2))
#' g2 <- as.vector(sm %*% g1) / rowSums(sm) + rnorm(n, sd = 0.5)
#'
#' res <- multiScaleAC(coords, g1, g2, h = c(0.05, 0.1, 0.2), perms = 199)
#' res <- compute_significance(res)
#' summary(res)
#'
#' # Whole-curve significance, rather than one test per bandwidth
#' res@sup_p_value
compute_significance <- function(result,
                                 method = c("permutation", "moment", "auto"),
                                 moment_dist = c("gamma", "saddlepoint", "edgeworth", "normal"),
                                 alternative = c("greater", "less", "two.sided"),
                                 alpha = 0.05,
                                 envelope_level = 0.95,
                                 estimator = c("add_one", "paper"),
                                 supremum = TRUE,
                                 sup_standardize = TRUE) {

  if (!inherits(result, "kernelcomm_result")) {
    stop("result must be a kernelcomm_result object")
  }
  method <- match.arg(method)
  moment_dist <- match.arg(moment_dist)
  alternative <- match.arg(alternative)
  estimator <- match.arg(estimator)

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single value strictly between 0 and 1")
  }
  if (!is.numeric(envelope_level) || length(envelope_level) != 1 ||
      envelope_level <= 0 || envelope_level >= 1) {
    stop("envelope_level must be a single value strictly between 0 and 1")
  }

  perms <- result@permutations
  observed <- result@observed
  n_perms <- ncol(perms)
  n_h <- length(observed)
  have_moments <- length(result@moment_var) == n_h

  if (length(result@p_values) > 0) {
    warning("significance was already computed; recomputing with the new parameters",
            call. = FALSE)
  }

  # Which bandwidths can be served by which route.
  # A bandwidth has permutations if its row is not all NA -- permute = "flagged" deliberately leaves
  # rows unpermuted.
  has_perm <- if (n_perms > 0) apply(perms, 1, function(r) any(is.finite(r))) else rep(FALSE, n_h)

  if (method %in% c("moment", "auto") && !have_moments) {
    stop("this result carries no analytical moments. Recompute with multiScaleAC(..., ",
         "moments = TRUE), which is the default.")
  }
  if (method == "permutation" && !any(has_perm)) {
    stop("this result has no permutations, so permutation inference is not possible. ",
         "Recompute with perms > 0, or use method = \"moment\".")
  }

  use_moment <- switch(
    method,
    permutation = rep(FALSE, n_h),
    moment      = rep(TRUE, n_h),
    # Moments where the diagnostics trust a continuous approximation, permutations elsewhere -- and
    # moments anyway where no permutations were drawn, since that is better than nothing.
    auto        = !.flag_untrustworthy_from_slots(result) | !has_perm
  )
  if (any(!use_moment & !has_perm)) {
    if (method == "auto" && have_moments) {
      use_moment[!has_perm] <- TRUE
    } else {
      stop(sum(!use_moment & !has_perm), " bandwidth(s) have neither permutations nor a usable ",
           "moment route. Recompute with perms > 0 or permute = \"all\".")
    }
  }

  # Permutation moments, where available, come from the permutations; otherwise from the exact
  # analytical values, which is strictly better since they are not estimates.
  perm_mean <- rep(0, n_h)
  perm_sd <- rep(NA_real_, n_h)
  if (any(has_perm)) {
    perm_mean[has_perm] <- rowMeans(perms[has_perm, , drop = FALSE])
    perm_sd[has_perm] <- apply(perms[has_perm, , drop = FALSE], 1, stats::sd)
  }
  if (have_moments) {
    # The exact mean is 0 and the exact sd is sqrt(moment_var); use them wherever there are no
    # permutations to estimate from.
    perm_sd[!has_perm] <- sqrt(result@moment_var[!has_perm])
  }

  p_values <- rep(NA_real_, n_h)
  p_method <- rep(NA_character_, n_h)

  # --- permutation route ------------------------------------------------------
  # Counts of permutations at least as extreme as the observed value. `>=` and `<=` include ties,
  # which is the conservative and correct choice for a permutation test.
  perm_idx <- which(!use_moment)
  if (length(perm_idx)) {
    r <- numeric(length(perm_idx))
    for (k in seq_along(perm_idx)) {
      i <- perm_idx[k]
      pr <- perms[i, ]
      pr <- pr[is.finite(pr)]
      r[k] <- switch(
        alternative,
        greater    = sum(pr >= observed[i]),
        less       = sum(pr <= observed[i]),
        # Two-sided by distance from the null centre. Preferred over doubling the smaller one-sided
        # p, which can exceed 1 and behaves poorly under a skewed null -- and this null is skewed
        # whenever the expression data is.
        two.sided  = sum(abs(pr - perm_mean[i]) >= abs(observed[i] - perm_mean[i]))
      )
    }
    p_values[perm_idx] <- switch(estimator,
                                 add_one = (r + 1) / (n_perms + 1),
                                 paper   = r / n_perms)
    p_method[perm_idx] <- "permutation"
  }

  # --- moment route -----------------------------------------------------------
  if (any(use_moment)) {
    mi <- which(use_moment)
    q <- observed[mi]
    v <- result@moment_var[mi]
    sk <- result@moment_skewness[mi]
    kt <- result@moment_excess_kurtosis[mi]

    pm <- switch(
      alternative,
      greater   = .moment_tail(q, v, sk, kt, method = moment_dist),
      less      = .moment_tail(q, v, sk, kt, method = moment_dist, lower = TRUE),
      # Two-sided from both tails of an asymmetric null, capped at 1. The null mean is exactly 0, so
      # "distance from the centre" is |q| and the two tails must be evaluated separately rather than
      # one being doubled.
      two.sided = pmin(1, .moment_tail(abs(q), v, sk, kt, method = moment_dist) +
                          .moment_tail(-abs(q), v, sk, kt, method = moment_dist, lower = TRUE))
    )
    p_values[mi] <- pm
    p_method[mi] <- paste0("moment:", moment_dist)
  }

  if (estimator == "paper" && any(p_values == 0, na.rm = TRUE)) {
    warning(sum(p_values == 0), " of ", n_h, " p-values are exactly 0, because no permutation was ",
            "as extreme as the observed value under estimator = \"paper\" (r/n). A p-value of 0 ",
            "breaks log transforms and multiple-testing procedures. Use estimator = \"add_one\", ",
            "which bounds them below at 1/(perms + 1), or increase perms.", call. = FALSE)
  }

  # Effect size only. See the warning in this function's documentation about pnorm().
  z_scores <- rep(NA_real_, n_h)
  ok <- !is.na(perm_sd) & perm_sd > 0
  z_scores[ok] <- (observed[ok] - perm_mean[ok]) / perm_sd[ok]

  # Null envelope. From the permutations where they exist, otherwise from the moment fit, so a
  # permutation-free result can still be plotted with a band.
  tail <- (1 - envelope_level) / 2
  null_lower <- rep(NA_real_, n_h)
  null_upper <- rep(NA_real_, n_h)
  if (any(has_perm)) {
    null_lower[has_perm] <- apply(perms[has_perm, , drop = FALSE], 1, stats::quantile,
                                 probs = tail, names = FALSE, na.rm = TRUE)
    null_upper[has_perm] <- apply(perms[has_perm, , drop = FALSE], 1, stats::quantile,
                                 probs = 1 - tail, names = FALSE, na.rm = TRUE)
  }
  if (any(!has_perm) && have_moments) {
    env <- .moment_envelope(result@moment_var[!has_perm], result@moment_skewness[!has_perm], tail)
    null_lower[!has_perm] <- env$lower
    null_upper[!has_perm] <- env$upper
  }

  result@p_values   <- p_values
  result@p_method   <- p_method
  result@z_scores   <- z_scores
  result@perm_mean  <- perm_mean
  result@perm_sd    <- perm_sd
  result@null_lower <- null_lower
  result@null_upper <- null_upper
  result@significant <- p_values < alpha

  # A stale adjustment from an earlier adjust_pvalues() call must not survive a recompute, or
  # summary() will show adjusted p-values that do not correspond to the current p_values.
  result@p_adjusted <- numeric(0)

  if (supremum && any(has_perm) && all(has_perm)) {
    sup <- .supremum_test(observed, perms, perm_mean, perm_sd,
                          standardize = sup_standardize,
                          estimator = estimator, alternative = alternative)
    result@sup_statistic <- sup$statistic
    result@sup_p_value <- sup$p_value
  } else {
    if (supremum && !all(has_perm)) {
      # The permutation supremum needs every bandwidth permuted, because it maximises across them.
      # A simulation-based version using the exact cross-bandwidth covariance is possible and is
      # noted in inst/validation, but it rests on a Gaussian copula for the joint null, which the
      # calibration work has not validated. Saying so beats shipping an untested global p-value.
      message("the supremum test needs permutations at every bandwidth; skipping it because ",
              "permute = \"flagged\" or \"none\" left some unpermuted")
    }
    result@sup_statistic <- numeric(0)
    result@sup_p_value <- numeric(0)
  }

  methods::validObject(result)
  result
}


#' Null envelope from moments
#'
#' Internal. Inverts the three-moment gamma fit to give envelope bounds when there
#' is no permutation distribution to take quantiles from.
#'
#' @param var,skew Exact moments on the scale of \eqn{I}.
#' @param tail Tail probability on each side.
#' @return List with \code{lower} and \code{upper}.
#' @keywords internal
#' @noRd
.moment_envelope <- function(var, skew, tail) {
  sd_i <- sqrt(var)
  z_of <- function(p, sk) {
    if (!is.finite(sk) || abs(sk) < 1e-8) return(stats::qnorm(p))
    k <- 4 / sk^2
    if (!is.finite(k) || k <= 0) return(stats::qnorm(p))
    # Standardised shifted gamma quantile.
    if (sk > 0) (stats::qgamma(p, shape = k) - k) / sqrt(k)
    else       -(stats::qgamma(1 - p, shape = k) - k) / sqrt(k)
  }
  lo <- vapply(seq_along(var), function(i) sd_i[i] * z_of(tail, skew[i]), numeric(1))
  hi <- vapply(seq_along(var), function(i) sd_i[i] * z_of(1 - tail, skew[i]), numeric(1))
  list(lower = lo, upper = hi)
}


#' Flag untrustworthy bandwidths from a stored result
#'
#' Internal. Same criteria as \code{.flag_untrustworthy()}, applied to the slots
#' of a result object rather than to the internal moment list.
#'
#' @param result A kernelcomm_result with moments computed.
#' @return Logical vector.
#' @keywords internal
#' @noRd
.flag_untrustworthy_from_slots <- function(result, max_leverage = 0.10, max_kurtosis = 20) {
  bad <- !is.finite(result@moment_var) | result@moment_var <= 0 |
         !is.finite(result@moment_skewness) | !is.finite(result@moment_excess_kurtosis)
  bad | result@leverage > max_leverage | abs(result@moment_excess_kurtosis) > max_kurtosis
}


#' Whole-curve supremum test
#'
#' The manuscript defines a global test alongside the per-bandwidth ones: it
#' "finds the bandwidth with greatest deviation from the mean curve for each
#' permutation and compares it to the largest deviation of the observed". It was
#' not implemented. It costs one O(n_h x perms) reduction of the permutation
#' matrix already in hand.
#'
#' This answers a different question from the per-bandwidth p-values, and a more
#' useful one for screening: is this curve as a whole inconsistent with complete
#' spatial randomness? It needs no multiplicity correction across bandwidths,
#' which matters because bandwidths of the same gene pair are strongly positively
#' dependent and correcting across them is hard to justify.
#'
#' @param observed Observed curve.
#' @param perms n_h x perms permutation matrix.
#' @param perm_mean,perm_sd Per-bandwidth permutation moments.
#' @param standardize Divide deviations by \code{perm_sd} before taking the max.
#' @param estimator "add_one" or "paper".
#' @param alternative Tail direction, applied to the deviation.
#'
#' @return List with \code{statistic} and \code{p_value}.
#' @keywords internal
#' @noRd
.supremum_test <- function(observed, perms, perm_mean, perm_sd,
                           standardize = TRUE, estimator = "add_one",
                           alternative = "greater") {
  n_perms <- ncol(perms)

  # Deviation from the null mean curve. Standardising makes each bandwidth contribute comparably;
  # without it the maximum is taken over quantities on different scales and is dominated by
  # whichever bandwidth happens to have the largest variance.
  scale_by <- if (standardize) {
    s <- perm_sd
    s[!is.finite(s) | s <= 0] <- NA_real_
    s
  } else {
    rep(1, length(observed))
  }

  dev_obs <- (observed - perm_mean) / scale_by
  dev_perm <- (perms - perm_mean) / scale_by

  # A bandwidth with zero permutation variance carries no information; exclude it rather than let
  # NA swallow the whole statistic.
  keep <- is.finite(dev_obs)
  if (!any(keep)) return(list(statistic = NA_real_, p_value = NA_real_))
  dev_obs <- dev_obs[keep]
  dev_perm <- dev_perm[keep, , drop = FALSE]

  reduce <- switch(alternative,
                   greater   = function(z) max(z),
                   less      = function(z) min(z),
                   two.sided = function(z) max(abs(z)))

  stat_obs <- reduce(dev_obs)
  stat_perm <- apply(dev_perm, 2, reduce)

  r <- if (alternative == "less") sum(stat_perm <= stat_obs) else sum(stat_perm >= stat_obs)
  p <- switch(estimator, add_one = (r + 1) / (n_perms + 1), paper = r / n_perms)

  list(statistic = stat_obs, p_value = p)
}


#' Adjust P-values for Multiple Testing
#'
#' Applies a multiple-testing correction to the per-bandwidth p-values of a
#' single result.
#'
#' @section What this does and does not correct:
#' This corrects \emph{across the bandwidths of one gene pair}. Those tests are
#' strongly positively dependent, because the kernels are nested on the same data
#' and the same permutations, so a Benjamini-Hochberg or Bonferroni correction
#' across them is conservative and not straightforward to interpret. If what you
#' want is a global statement about the curve, prefer \code{sup_p_value} from
#' \code{\link{compute_significance}}, which tests the whole curve directly and
#' needs no correction across bandwidths.
#'
#' The manuscript applies its FDR threshold across \emph{gene pairs}, not across
#' bandwidths. To reproduce that, collect one p-value per pair, for example
#' \code{sup_p_value} or the value at a chosen bandwidth, and correct across the
#' collection, as with the output of \code{\link{multiScaleAC_batch}}.
#'
#' @param result Object of class \code{\linkS4class{kernelcomm_result}} with
#'   significance already computed.
#' @param method Character, passed to \code{\link[stats]{p.adjust}}. Note that
#'   \code{"BH"} and \code{"fdr"} are the same procedure under two names.
#'
#' @return The same object with \code{p_adjusted} populated.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 150
#' coords <- matrix(runif(2 * n), ncol = 2)
#' res <- multiScaleAC(coords, rnorm(n), rnorm(n), h = c(0.05, 0.1, 0.2), perms = 99)
#' res <- compute_significance(res)
#' res <- adjust_pvalues(res, method = "BH")
#' summary(res)$p_adjusted
adjust_pvalues <- function(result, method = "BH") {
  if (!inherits(result, "kernelcomm_result")) {
    stop("result must be a kernelcomm_result object")
  }
  if (length(result@p_values) == 0) {
    stop("call compute_significance() before adjust_pvalues()")
  }
  method <- match.arg(method, stats::p.adjust.methods)

  result@p_adjusted <- stats::p.adjust(result@p_values, method = method)
  result
}
