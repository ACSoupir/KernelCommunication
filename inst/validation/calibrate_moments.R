#!/usr/bin/env Rscript
#
# Calibration of the moment-based tail approximations against permutation ground truth.
#
# The moments in R/moments.R are exact; mapping them to a tail probability is not. This script
# decides which mapping becomes the default, by measurement. It does NOT assume an answer.
#
# METHOD. For each configuration, draw a large permutation null. For a set of target upper-tail
# probabilities, read off the empirical quantile at that probability, then ask each approximation what
# tail probability it assigns to that quantile. A perfectly calibrated method returns the target, so
# the reported ratio is 1. Ratios below 1 are ANTI-CONSERVATIVE: the method understates how easily the
# null reaches that value, which manufactures false positives. That is the dangerous direction, and it
# is why the ratio is reported rather than an absolute error.
#
# The targets go down to 1e-4 because an FDR screen over the ~1100 ligand-receptor pairs per sample in
# the companion manuscript operates there. A method calibrated at p = 0.05 and wrong at p = 1e-4 is
# not useful for the actual workload.
#
# Run from the package root:
#   Rscript inst/validation/calibrate_moments.R [n_perms]

suppressMessages(library(multiScaleAC))

args <- commandArgs(TRUE)
N_PERM <- if (length(args)) as.integer(args[1]) else 200000L
TARGETS <- c(0.1, 0.01, 1e-3, 1e-4)
METHODS <- c("saddlepoint", "gamma", "edgeworth", "normal")

# The smallest target must be resolvable: at least ~20 permutations in the tail.
stopifnot(N_PERM * min(TARGETS) >= 20)

message(sprintf("Calibration with %d permutations per configuration\n", N_PERM))

# ---------------------------------------------------------------------------
# Configurations, spanning benign to adversarial. The point is to find where the approximations break,
# so the set deliberately includes cases they should struggle with.
# ---------------------------------------------------------------------------
make_configs <- function() {
  cfg <- list()

  set.seed(101)
  n <- 300
  co <- matrix(runif(2 * n), ncol = 2)
  cfg$gaussian_genes <- list(
    label = "gaussian genes, uniform points",
    coords = co, g1 = rnorm(n), g2 = rnorm(n), h = c(0.02, 0.08, 0.25)
  )

  set.seed(102)
  n <- 300
  co <- rbind(matrix(rnorm(2 * (n - 20), 0.5, 0.12), ncol = 2), matrix(runif(40), ncol = 2))
  cfg$zero_inflated <- list(
    label = "85% zero gene2, clustered points",
    coords = co, g1 = rlnorm(n, 0, 1),
    g2 = c(rep(0, n - 45), rlnorm(45, 1.5, 1))[sample.int(n)],
    h = c(0.02, 0.08, 0.25)
  )

  # Deliberately adversarial: very few non-zero cells, so the statistic is dominated by a handful of
  # terms and the combinatorial CLT has little to work with.
  set.seed(103)
  n <- 250
  co <- rbind(matrix(rnorm(2 * (n - 10), 0.5, 0.1), ncol = 2), matrix(runif(20, -1, 2), ncol = 2))
  cfg$extreme_sparse <- list(
    label = "94% zero gene2, lognormal gene1 (adversarial)",
    coords = co, g1 = rlnorm(n, 0, 1.4),
    g2 = c(rep(0, n - 15), rlnorm(15, 2, 1))[sample.int(n)],
    h = c(0.02, 0.1)
  )

  # Real data, subsampled for speed.
  data(ccRCCExample, envir = environment())
  set.seed(104)
  idx <- sort(sample.int(nrow(ccRCCExample), 500))
  sub <- ccRCCExample[idx, ]
  cfg$ccrcc <- list(
    label = "ccRCCExample subsample, CD274 vs PDCD1",
    coords = as.matrix(sub[, c("x", "y")]), g1 = sub$CD274, g2 = sub$PDCD1,
    h = c(0.02, 0.08, 0.25)
  )

  cfg
}

rows <- list()

for (cname in names(configs <- make_configs())) {
  cf <- configs[[cname]]
  message(sprintf("--- %s ---", cf$label))

  res <- multiScaleAC(cf$coords, cf$g1, cf$g2, h = cf$h,
                      perms = N_PERM, seed = 20260908, verbose = FALSE)

  g2s <- as.vector(scale(cf$g2))
  g1s <- as.vector(scale(cf$g1))
  denom <- sum(g1s^2)
  # Rebuild the internals the approximations need. Using the package's own kernel pass rather than a
  # separate implementation, so this measures the shipped code.
  pass <- multiScaleAC:::moran_kernel_pass(as.matrix(cf$coords[, 1:2]), g1s, cf$h, eps = 0)
  mom <- multiScaleAC:::.permutation_moments(pass$V, g2s, denom)

  for (i in seq_along(cf$h)) {
    null <- res@permutations[i, ]
    for (tg in TARGETS) {
      # Empirical quantile at this upper-tail probability.
      q <- as.numeric(stats::quantile(null, probs = 1 - tg, names = FALSE))
      # Recompute the achieved empirical tail at that q, since ties and discreteness mean it is not
      # exactly `tg`. This is the reference the approximations are scored against.
      emp <- mean(null >= q)

      one <- list(config = cname, h = cf$h[i], target = tg, empirical = emp,
                  skewness = mom$skewness[i], kurtosis = mom$excess_kurtosis[i],
                  leverage = mom$leverage[i],
                  flagged = multiScaleAC:::.flag_untrustworthy(mom)[i])
      for (m in METHODS) {
        p <- multiScaleAC:::.moment_tail(
          q = rep(q, length(cf$h)),
          var = mom$var, skewness = mom$skewness, excess_kurtosis = mom$excess_kurtosis,
          method = m, a = mom$.a, b = g2s, denom = denom
        )[i]
        one[[m]] <- p / emp
      }
      rows[[length(rows) + 1L]] <- one
    }
  }
}

tab <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))

cat("\n\n=====================================================================\n")
cat(sprintf("CALIBRATION: ratio of approximate to empirical upper-tail p (%d perms)\n", N_PERM))
cat("1.00 = perfectly calibrated.  BELOW 1 = anti-conservative = false positives.\n")
cat("=====================================================================\n\n")

fmt <- function(x) ifelse(is.na(x), "     -", sprintf("%6.2f", x))
cat(sprintf("%-16s %-6s %-8s %-7s %-7s | %s\n", "config", "h", "target", "skew", "lev",
            paste(sprintf("%-11s", METHODS), collapse = "")))
for (k in seq_len(nrow(tab))) {
  r <- tab[k, ]
  cat(sprintf("%-16s %-6g %-8.0e %-7.2f %-7.3f | %s%s\n",
              r$config, r$h, r$target, r$skewness, r$leverage,
              paste(sprintf("%-11s", vapply(METHODS, function(m) fmt(r[[m]]), character(1))),
                    collapse = ""),
              if (isTRUE(r$flagged)) "  <- flagged" else ""))
}

cat("\n--- summary over UNFLAGGED configurations only ---\n")
cat("(flagged bandwidths fall back to permutations, so their calibration does not decide the default)\n\n")
un <- tab[!tab$flagged, , drop = FALSE]
cat(sprintf("%-13s %-9s %-9s %-9s %-9s %s\n", "method", "median", "min", "max",
            "worst_low", "n_anticons"))
for (m in METHODS) {
  v <- un[[m]][is.finite(un[[m]])]
  cat(sprintf("%-13s %-9.3f %-9.3f %-9.3f %-9.3f %d/%d\n", m,
              stats::median(v), min(v), max(v), min(v), sum(v < 0.5), length(v)))
}

cat("\n--- deep tail only (target <= 1e-3), unflagged ---\n\n")
deep <- un[un$target <= 1e-3, , drop = FALSE]
for (m in METHODS) {
  v <- deep[[m]][is.finite(deep[[m]])]
  if (!length(v)) next
  cat(sprintf("%-13s median %.3f  range [%.3f, %.3f]\n", m, stats::median(v), min(v), max(v)))
}

outdir <- file.path("inst", "validation", "results")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, file.path(outdir, "calibration.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s\n", file.path(outdir, "calibration.csv")))
