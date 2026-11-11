# Generate the golden reference fixtures.
#
# These fixtures are the regression oracle for the whole package. They are produced by
# soupir_moran_multi() and NOTHING ELSE. That function is pure base R, so this script needs no
# compiled code and cannot be invalidated by a package rename or a rewrite of the compute core.
# soupir_moran_multi() must never be modified; tests/testthat/test-oracle-integrity.R enforces that.
#
# Run from the package root:
#   Rscript tests/testthat/fixtures/make_golden.R
#
# BANDWIDTH UNITS -- read before changing any `h` below.
# ccRCCExample coordinates are in millimetres, spanning 1.02 x 1.02 mm with a median
# nearest-neighbour distance of 0.0073 mm. The manuscript's bandwidth range of 10-250 is in
# micrometres, i.e. h = 0.01 to 0.25 here. Earlier validation in this repository used h = 10..250,
# which drives every off-diagonal weight to >= 0.9999 and makes the kernel a constant matrix -- so
# it appeared to pass while testing nothing. Every case below records
# mean_rowsum so test-fixture-sanity.R can fail if a bandwidth is degenerate.

# Locate the package root by walking up until DESCRIPTION is found, so this script works whether it
# is run from the package root, from tests/, or from the fixtures directory itself.
find_pkg_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = TRUE)
  for (i in 1:6) {
    if (file.exists(file.path(d, "DESCRIPTION"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("could not locate the package root (no DESCRIPTION found above ", start, ")")
}
pkg_root   <- find_pkg_root()
fixture_dir <- file.path(pkg_root, "tests", "testthat", "fixtures")
dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)

# Source the reference implementation directly rather than loading the package, so fixtures never
# depend on a built package or on whatever the rest of R/ currently does.
source(file.path(pkg_root, "R", "soupir_moran_multi.R"))
stopifnot(is.function(soupir_moran_multi))

# Canonical bandwidth grid, in the coordinate units of the data (mm for ccRCCExample; the synthetic
# cases are generated on a comparable unit square so the same grid applies). Reused verbatim by the
# examples and the vignette so the three cannot drift apart.
H_UNIT <- c(0.01, 0.02, 0.04, 0.08, 0.15, 0.25)

PERMS <- 199L   # small enough to store, large enough that a permutation-matrix bug shows up
SEED  <- 333L

# ---------------------------------------------------------------------------
# Diagnostics recorded alongside each case: how many cells actually carry weight
# at each bandwidth. This is what makes a degenerate fixture detectable.
# ---------------------------------------------------------------------------
kernel_diagnostics <- function(coords, h_vec) {
  d2 <- as.matrix(stats::dist(coords))^2
  t(vapply(h_vec, function(h) {
    K <- exp(-d2 / (2 * h^2))
    diag(K) <- 0
    rs <- rowSums(K)
    c(mean_rowsum = mean(rs), min_rowsum = min(rs),
      max_offdiag = max(K), n_isolated = sum(rs == 0))
  }, numeric(4)))
}

make_case <- function(name, coords, gene1, gene2, h_vec, note) {
  message(sprintf("  %-22s n=%-6d n_h=%d", name, nrow(coords), length(h_vec)))
  res  <- soupir_moran_multi(coords, gene1, gene2, h_vec = h_vec, perms = PERMS, seed = SEED)
  diag <- kernel_diagnostics(coords, h_vec)
  list(name = name, note = note,
       coords = coords, gene1 = gene1, gene2 = gene2,
       h = h_vec, perms = PERMS, seed = SEED,
       observed = res$observed, permutations = res$permutations,
       diagnostics = diag)
}

cases <- list()

# --- 1. plain synthetic: uniform points, independent normal genes -----------
set.seed(101)
n <- 250
co <- matrix(runif(2 * n), ncol = 2)
cases$synth_uniform <- make_case(
  "synth_uniform", co, rnorm(n), rnorm(n), H_UNIT,
  "Uniform points on the unit square, independent standard normal genes. Baseline null case."
)

# --- 2. synthetic with genuine spatial cross-correlation -------------------
# gene2 is a spatially smoothed function of gene1, so I(h) is clearly non-zero and a sign or
# normalisation error cannot hide behind a near-zero statistic.
set.seed(202)
n <- 300
co <- matrix(runif(2 * n), ncol = 2)
g1 <- rnorm(n)
sm <- exp(-as.matrix(stats::dist(co))^2 / (2 * 0.08^2))
g2 <- as.vector(sm %*% g1) / rowSums(sm) + rnorm(n, sd = 0.25)
cases$synth_correlated <- make_case(
  "synth_correlated", co, g1, g2, H_UNIT,
  "gene2 is a spatially smoothed version of gene1 (scale 0.08) plus noise, so I(h) is strongly non-zero."
)

# --- 3. skewed, zero-inflated gene2 (realistic single-cell expression) -----
# This is the regime where the permutation null is skewed and a normal approximation fails, so it is
# the case the analytical-moment work is calibrated against.
set.seed(303)
n <- 300
co <- rbind(matrix(rnorm(2 * (n - 20), 0.5, 0.12), ncol = 2),
            matrix(runif(40), ncol = 2))
g1 <- rlnorm(n, 0, 1.0)
g2 <- c(rep(0, n - 45), rlnorm(45, 1.5, 1.0))[sample.int(n)]
cases$synth_zeroinflated <- make_case(
  "synth_zeroinflated", co, g1, g2, H_UNIT,
  "Clustered points, lognormal gene1, 85% zero gene2. Skewed permutation null; drives the moment work."
)

# --- 4. real data: ccRCCExample subsample ----------------------------------
# Subsampled so the fixture stays small and the reference (which builds a dense N x N matrix) stays
# fast. Real CD274 / PDCD1 SCTransform residuals, real millimetre coordinates.
load(file.path(pkg_root, "data", "ccRCCExample.rda"))
set.seed(404)
idx <- sort(sample.int(nrow(ccRCCExample), 1200))
sub <- ccRCCExample[idx, ]
cases$ccrcc_subsample <- make_case(
  "ccrcc_subsample", as.matrix(sub[, c("x", "y")]), sub$CD274, sub$PDCD1, H_UNIT,
  "1200-cell subsample of ccRCCExample; CD274 vs PDCD1; coordinates in millimetres."
)

# --- 5. isolated cell: pins the NaN / zero-rowsum semantics ----------------
# One point placed far away. At the smallest bandwidths exp(-d^2/2h^2) underflows to exactly 0 for
# that row, so its row sum is 0, gene1/rowsum is Inf, v_h contains NaN, and the reference's
# `observed[is.nan(observed)] <- 0` zeroes THE WHOLE BANDWIDTH -- not just that cell. Any
# reimplementation must reproduce that, so it is captured explicitly here.
set.seed(505)
n <- 120
co <- rbind(matrix(runif(2 * (n - 1)), ncol = 2), c(500, 500))
cases$isolated_cell <- make_case(
  "isolated_cell", co, rnorm(n), rnorm(n), c(0.02, 0.05, 0.1, 0.25),
  paste("One point at (500,500), far outside the unit square. Its row sum underflows to 0 at every",
        "bandwidth here, so the reference returns 0 for the entire bandwidth. Pins nan_policy='zero'.")
)

# --- 6. zero-variance gene: scale() yields NaN throughout ------------------
set.seed(606)
n <- 100
co <- matrix(runif(2 * n), ncol = 2)
cases$zero_variance <- make_case(
  "zero_variance", co, rep(1.5, n), rnorm(n), c(0.05, 0.1, 0.25),
  "gene1 is constant, so scale() returns NaN and the reference returns 0 for every bandwidth."
)

golden <- list(
  cases       = cases,
  H_UNIT      = H_UNIT,
  perms       = PERMS,
  seed        = SEED,
  generated_by = "soupir_moran_multi() -- reference implementation, base R only",
  r_version   = R.version.string,
  platform    = R.version$platform
)

out <- file.path(fixture_dir, "golden.rds")
saveRDS(golden, out, version = 2)
message(sprintf("\nwrote %s (%.0f KB, %d cases)", out, file.size(out) / 1024, length(cases)))

# Report the diagnostics so a degenerate bandwidth is visible at generation time, not just in tests.
message("\nEffective neighbours (mean row sum) per case and bandwidth:")
for (cs in cases) {
  message(sprintf("  %-20s N=%-5d %s", cs$name, nrow(cs$coords),
                  paste(sprintf("h=%g:%.1f", cs$h, cs$diagnostics[, "mean_rowsum"]),
                        collapse = "  ")))
}
