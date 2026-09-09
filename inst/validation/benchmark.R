#!/usr/bin/env Rscript
#
# Time and peak memory, against the reference implementation.
#
# This replaces three earlier root-level scripts which could not run at all: they accessed the result
# with `$`, and multiScaleAC() returns an S4 object. It also replaces validation_results.txt, whose
# headline "71.65x speedup" compared four separate soupir_moran_bv() calls, each rebuilding the
# kernel, against one multiScaleAC() call that builds the distance data once and reuses it across all
# bandwidths -- an amortisation gain reported as a language gain.
#
# Peak memory cannot be measured meaningfully from inside R, because the allocator does not return
# freed pages to the OS promptly. Each configuration therefore runs in a fresh process under
# /usr/bin/time, and this script orchestrates that.
#
# Run from the package root:
#   Rscript inst/validation/benchmark.R

WORKER <- Sys.getenv("MSAC_BENCH_WORKER", "")

# ---------------------------------------------------------------------------
# Worker: one configuration in a clean process.
# ---------------------------------------------------------------------------
if (nzchar(WORKER)) {
  args <- commandArgs(TRUE)
  mode <- args[1]; N <- as.integer(args[2])
  suppressMessages(library(multiScaleAC))
  data(ccRCCExample, package = "multiScaleAC", envir = environment())

  set.seed(1)
  if (N <= nrow(ccRCCExample)) {
    d <- ccRCCExample[sort(sample.int(nrow(ccRCCExample), N)), ]
  } else {
    # Replicate with a small jitter to reach sizes the bundled sample cannot, so the memory scaling
    # claim is tested rather than extrapolated.
    idx <- sample.int(nrow(ccRCCExample), N, replace = TRUE)
    d <- ccRCCExample[idx, ]
    d$x <- d$x + stats::rnorm(N, 0, 1e-4)
    d$y <- d$y + stats::rnorm(N, 0, 1e-4)
  }
  co <- as.matrix(d[, c("x", "y")]); g1 <- d$CD274; g2 <- d$PDCD1
  hv <- c(0.01, 0.02, 0.04, 0.08, 0.15, 0.25)
  PERMS <- 999L

  t0 <- proc.time()[3]
  obs <- switch(
    mode,
    reference = soupir_moran_multi(co, g1, g2, h_vec = hv, perms = PERMS, seed = 333)$observed,
    permutation = multiScaleAC(co, g1, g2, h = hv, perms = PERMS, seed = 333,
                               verbose = FALSE)@observed,
    moment = multiScaleAC(co, g1, g2, h = hv, permute = "none", seed = 333,
                          verbose = FALSE)@observed,
    stop("unknown mode: ", mode)
  )
  el <- proc.time()[3] - t0
  cat(sprintf("RESULT\t%s\t%d\t%.3f\t%.8e\t%d\n",
              mode, N, el, obs[1], multiScaleAC:::moran_max_threads()))
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# Orchestrator.
# ---------------------------------------------------------------------------
self <- "inst/validation/benchmark.R"
if (!file.exists(self)) stop("run this from the package root")

timer <- if (file.exists("/usr/bin/time")) "/usr/bin/time -l" else NA_character_
rss_pattern <- "maximum resident set size"

run_one <- function(mode, N) {
  cmd <- sprintf("MSAC_BENCH_WORKER=1 %s Rscript --vanilla %s %s %d 2>&1",
                 if (is.na(timer)) "" else timer, self, mode, N)
  out <- suppressWarnings(system(cmd, intern = TRUE))
  res <- grep("^RESULT\t", out, value = TRUE)
  if (!length(res)) {
    return(list(mode = mode, N = N, seconds = NA_real_, rss_mb = NA_real_,
                observed = NA_real_, threads = NA_integer_,
                note = paste(utils::tail(out, 2), collapse = " ")))
  }
  f <- strsplit(res[1], "\t")[[1]]
  rssl <- grep(rss_pattern, out, value = TRUE)
  rss <- if (length(rssl)) {
    as.numeric(sub("^\\s*([0-9]+).*$", "\\1", rssl[1])) / 1024^2
  } else NA_real_
  list(mode = f[2], N = as.integer(f[3]), seconds = as.numeric(f[4]),
       rss_mb = rss, observed = as.numeric(f[5]), threads = as.integer(f[6]), note = "")
}

CONFIGS <- list(
  list("reference",   2000L), list("permutation", 2000L), list("moment", 2000L),
  list("reference",   7122L), list("permutation", 7122L), list("moment", 7122L),
  # The reference is deliberately not run at 20000: it would need three N x N doubles, about 9.6 GB,
  # and would thrash. That it cannot run is itself the result.
  list("permutation", 20000L), list("moment", 20000L)
)

rows <- lapply(CONFIGS, function(cf) {
  message(sprintf("running %-12s N = %d", cf[[1]], cf[[2]]))
  run_one(cf[[1]], cf[[2]])
})
tab <- do.call(rbind, lapply(rows, as.data.frame))

cat("\n=====================================================================\n")
cat("BENCHMARK: six bandwidths, 999 permutations where applicable\n")
cat(sprintf("OpenMP threads available: %s\n",
            paste(unique(stats::na.omit(tab$threads)), collapse = ", ")))
cat("=====================================================================\n\n")
cat(sprintf("%-13s %-8s %-10s %-11s %s\n", "mode", "N", "seconds", "peak RSS MB", "I(h1)"))
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  cat(sprintf("%-13s %-8d %-10.2f %-11.0f %.6e%s\n", r$mode, r$N, r$seconds, r$rss_mb,
              r$observed, if (nzchar(r$note)) paste("  <-", r$note) else ""))
}

cat("\n--- agreement with the reference, where the reference could run ---\n")
for (nn in unique(tab$N)) {
  ref <- tab$observed[tab$N == nn & tab$mode == "reference"]
  if (!length(ref) || is.na(ref)) next
  for (md in c("permutation", "moment")) {
    v <- tab$observed[tab$N == nn & tab$mode == md]
    if (length(v) && !is.na(v)) {
      cat(sprintf("  N = %-6d %-12s relative difference %.3e\n", nn, md, abs(v - ref) / abs(ref)))
    }
  }
}

cat("\n--- memory scaling ---\n")
for (md in c("reference", "permutation")) {
  s <- tab[tab$mode == md & !is.na(tab$rss_mb), ]
  if (nrow(s) < 2) next
  cat(sprintf("  %-12s N %d -> %d : %.0f MB -> %.0f MB (%.1fx for %.1fx the cells)\n",
              md, min(s$N), max(s$N),
              s$rss_mb[which.min(s$N)], s$rss_mb[which.max(s$N)],
              s$rss_mb[which.max(s$N)] / s$rss_mb[which.min(s$N)],
              max(s$N) / min(s$N)))
}
cat("\n  An O(N^2) implementation must grow roughly with N^2; the tiled core is flat in N,\n")
cat("  because its working set is O(N * length(h)) and the rest is the R process itself.\n")

outdir <- file.path("inst", "validation", "results")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, file.path(outdir, "benchmark.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s\n", file.path(outdir, "benchmark.csv")))
