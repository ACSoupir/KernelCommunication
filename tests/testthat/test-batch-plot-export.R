# Batch processing, plotting, and the mxfda bridge.
#
# Every plot function is called on a real result, because several of them had never been executed:
# plot_batch_comparison() failed on its own default argument as soon as any cell was significant.

mk_result <- function(seed = 1, n = 120, h = c(0.05, 0.1, 0.2), perms = 99, related = FALSE) {
  set.seed(seed)
  coords <- matrix(runif(2 * n), ncol = 2)
  g1 <- rnorm(n)
  g2 <- if (related) {
    sm <- exp(-as.matrix(stats::dist(coords))^2 / (2 * 0.1^2))
    as.vector(sm %*% g1) / rowSums(sm) + rnorm(n, sd = 0.4)
  } else rnorm(n)
  compute_significance(
    multiScaleAC(coords, g1, g2, h = h, perms = perms, seed = seed, verbose = FALSE)
  )
}

mk_batch_data <- function(k = 3, n = 100) {
  lapply(seq_len(k), function(i) {
    set.seed(100 + i)
    data.frame(x = runif(n), y = runif(n), CD274 = rnorm(n), PDCD1 = rnorm(n))
  })
}

# ---------------------------------------------------------------------------
# Batch
# ---------------------------------------------------------------------------

test_that("batch processing runs and labels samples", {
  d <- mk_batch_data(3)
  names(d) <- c("s1", "s2", "s3")
  res <- multiScaleAC_batch(d, h_vec = c(0.05, 0.1), perms = 49, verbose = FALSE)
  expect_length(res, 3L)
  expect_identical(names(res), c("s1", "s2", "s3"))
  expect_identical(vapply(res, function(r) r@sample_id, character(1)),
                   c(s1 = "s1", s2 = "s2", s3 = "s3"))
  # Gene names come from the column arguments, which is more useful than the caller expression here.
  expect_identical(res[[1]]@gene1_name, "CD274")
  # Distinct seeds per sample, so samples are not identically permuted.
  expect_false(identical(res[[1]]@permutations, res[[2]]@permutations))
})

test_that("missing columns are reported before any computation starts", {
  # Columns used to be subset before being checked, so this produced "subscript out of bounds" and
  # the informative messages were unreachable. It also failed only on reaching that sample.
  d <- mk_batch_data(2)
  d[[2]]$CD274 <- NULL
  expect_error(multiScaleAC_batch(d, h_vec = 0.1, perms = 9, verbose = FALSE),
               regexp = "missing column", ignore.case = TRUE)
  # The message must name the offending column.
  expect_error(multiScaleAC_batch(d, h_vec = 0.1, perms = 9, verbose = FALSE),
               regexp = "CD274")
})

test_that("duplicate sample ids are rejected", {
  # Checkpoint files are keyed by sample id, so duplicates would overwrite each other.
  d <- mk_batch_data(2)
  expect_error(
    multiScaleAC_batch(d, sample_ids = c("a", "a"), h_vec = 0.1, perms = 9, verbose = FALSE),
    regexp = "unique"
  )
})

test_that("checkpointing resumes without recomputing and cleans up on success", {
  d <- mk_batch_data(3)
  names(d) <- c("a", "b", "c")
  ck <- tempfile("msac_ckpt_a_")
  on.exit(unlink(ck, recursive = TRUE), add = TRUE)

  # First pass over a subset, to leave checkpoint files behind.
  part <- multiScaleAC_batch(d[1:2], h_vec = c(0.05, 0.1), perms = 49,
                             checkpoint_dir = ck, verbose = FALSE)
  expect_length(list.files(ck), 0L)   # complete batch cleans up after itself

  # Simulate an interrupted batch: write one sample's checkpoint by hand, then resume the full set.
  dir.create(ck, showWarnings = FALSE, recursive = TRUE)
  saveRDS(part[["a"]], file.path(ck, "msac_a.rds"))
  full <- multiScaleAC_batch(d, h_vec = c(0.05, 0.1), perms = 49,
                             checkpoint_dir = ck, verbose = FALSE)
  expect_length(full, 3L)
  # The resumed sample must be the checkpointed object, not a recomputation.
  expect_identical(full[["a"]]@observed, part[["a"]]@observed)
  expect_false(any(vapply(full, is.null, logical(1))))
  expect_length(list.files(ck), 0L)
})

test_that("a fully checkpointed batch does not reprocess anything in reverse", {
  # `for (i in start_idx:n_samples)` counted DOWN when start_idx was n_samples + 1, silently
  # reprocessing the whole batch backwards. The loop is now over an explicit index set.
  d <- mk_batch_data(2)
  names(d) <- c("p", "q")
  ck <- tempfile("msac_ckpt_b_")
  dir.create(ck, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(ck, recursive = TRUE), add = TRUE)

  first <- multiScaleAC_batch(d, h_vec = 0.1, perms = 19, verbose = FALSE)
  saveRDS(first[["p"]], file.path(ck, "msac_p.rds"))
  saveRDS(first[["q"]], file.path(ck, "msac_q.rds"))

  again <- multiScaleAC_batch(d, h_vec = 0.1, perms = 19, checkpoint_dir = ck, verbose = FALSE)
  expect_length(again, 2L)
  expect_identical(names(again), c("p", "q"))
  expect_identical(again[["p"]]@observed, first[["p"]]@observed)
  expect_identical(again[["q"]]@observed, first[["q"]]@observed)
})

test_that("parallel and sequential batches agree", {
  skip_on_os("windows")
  d <- mk_batch_data(3)
  seqr <- multiScaleAC_batch(d, h_vec = c(0.05, 0.1), perms = 49, n_cores = 1, verbose = FALSE)
  parr <- multiScaleAC_batch(d, h_vec = c(0.05, 0.1), perms = 49, n_cores = 2, verbose = FALSE)
  for (i in seq_along(seqr)) {
    expect_equal(parr[[i]]@observed, seqr[[i]]@observed, tolerance = 1e-14)
    expect_equal(parr[[i]]@permutations, seqr[[i]]@permutations, tolerance = 1e-14)
  }
})

test_that("batch forwards extra arguments to multiScaleAC", {
  d <- mk_batch_data(2)
  # A uniform kernel at h = 0.05 on 100 points leaves some cells with no neighbours at all, so the
  # zero-row-sum warning is expected here and is asserted rather than suppressed: it confirms that
  # nan_policy reaches through the batch driver to each sample.
  expect_warning(
    res <- multiScaleAC_batch(d, h_vec = c(0.05, 0.1), perms = 0, permute = "none",
                              kernel = "uniform", verbose = FALSE),
    regexp = "no neighbours|no neighbors", ignore.case = TRUE
  )
  expect_identical(ncol(res[[1]]@permutations), 0L)
  # A uniform kernel really was used: it differs from the Gaussian default.
  gauss <- multiScaleAC_batch(d, h_vec = c(0.05, 0.1), perms = 0, permute = "none",
                              verbose = FALSE)
  expect_false(isTRUE(all.equal(res[[1]]@observed, gauss[[1]]@observed)))
})

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

test_that("plot_scale_curve works for every type", {
  res <- mk_result(related = TRUE)
  for (ty in c("observed", "difference", "zscore")) {
    p <- plot_scale_curve(res, type = ty)
    expect_s3_class(p, "ggplot")
    expect_silent(invisible(ggplot2::ggplot_build(p)))
  }
})

test_that("show_permutations is implemented, not merely documented", {
  # This argument previously had no implementation anywhere in the function body.
  res <- mk_result()
  plain <- plot_scale_curve(res, show_permutations = FALSE)
  withp <- plot_scale_curve(res, show_permutations = TRUE, show_permutations_n = 20)
  expect_gt(length(withp$layers), length(plain$layers))
  expect_silent(invisible(ggplot2::ggplot_build(withp)))
  # The cap is respected rather than drawing every permutation.
  perm_layer <- withp$layers[[length(plain$layers) + 1L]]
  expect_lte(length(unique(perm_layer$data$perm)), 20L)
})

test_that("plot_scale_curve draws no normal-theory reference lines", {
  # +-1.96 lines invited reading a normal threshold off a skewed permutation null. The envelope now
  # carries that information, taken from the permutations themselves.
  res <- mk_result()
  b <- ggplot2::ggplot_build(plot_scale_curve(res, type = "zscore"))
  yints <- unlist(lapply(b$plot$layers, function(l) l$data$yintercept))
  expect_false(any(abs(yints - 1.96) < 1e-9, na.rm = TRUE))
})

test_that("plot_permutation_distribution works", {
  res <- mk_result()
  p <- plot_permutation_distribution(res, h_index = 2)
  expect_s3_class(p, "ggplot")
  expect_silent(invisible(ggplot2::ggplot_build(p)))
})

test_that("plot_batch_comparison works on its default value_var when cells are significant", {
  # THE regression test for this phase. The significance-asterisk layer passed a data frame with no
  # `value` column while inheriting `fill = value` from the main plot, so this errored with
  # "object 'value' not found" as soon as anything was significant -- on the default argument.
  results <- lapply(1:4, function(i) mk_result(seed = i, related = TRUE))
  names(results) <- paste0("s", 1:4)
  expect_true(any(vapply(results, function(r) any(r@significant), logical(1))),
              info = "the fixture must contain significant cells for this test to bite")

  p <- plot_batch_comparison(results)
  expect_s3_class(p, "ggplot")
  expect_silent(invisible(ggplot2::ggplot_build(p)))

  for (vv in c("observed", "p_value")) {
    expect_silent(invisible(ggplot2::ggplot_build(plot_batch_comparison(results, value_var = vv))))
  }
})

test_that("cluster_rows actually clusters", {
  # It was documented as clustering by similarity and instead sorted by max absolute z-score, and only
  # when value_var happened to be "z_score".
  results <- lapply(1:5, function(i) mk_result(seed = i, related = (i %% 2 == 0)))
  names(results) <- paste0("s", 1:5)

  clustered <- plot_batch_comparison(results, cluster_rows = TRUE)
  plain <- plot_batch_comparison(results, cluster_rows = FALSE)
  expect_identical(levels(plain$data$sample_id), paste0("s", 1:5))
  # Clustering must produce a permutation of the samples, and generally a different one.
  expect_setequal(levels(clustered$data$sample_id), paste0("s", 1:5))
  expect_silent(invisible(ggplot2::ggplot_build(clustered)))

  # It must work for a non-z_score value_var too, which the old branch ignored entirely.
  expect_silent(invisible(ggplot2::ggplot_build(
    plot_batch_comparison(results, value_var = "observed", cluster_rows = TRUE))))
})

test_that("plot_significance_summary honours alpha and both sort orders", {
  results <- lapply(1:4, function(i) mk_result(seed = i, related = TRUE))
  names(results) <- c("d", "b", "c", "a")

  # alpha used to appear only in the title while the bars came from a slot computed with a possibly
  # different alpha, so the label could contradict the plot.
  strict <- plot_significance_summary(results, alpha = 1e-6)
  loose <- plot_significance_summary(results, alpha = 0.5)
  expect_lte(sum(strict$data$n_significant), sum(loose$data$n_significant))

  # sort_by = "sample_id" was documented with no branch implementing it.
  by_id <- plot_significance_summary(results, sort_by = "sample_id")
  expect_identical(levels(by_id$data$sample_id), c("a", "b", "c", "d"))
  by_count <- plot_significance_summary(results, sort_by = "count")
  expect_true(all(diff(by_count$data$n_significant) <= 0))
  expect_silent(invisible(ggplot2::ggplot_build(by_id)))
})

test_that("plot dispatch through the S3 method works", {
  res <- mk_result()
  expect_s3_class(plot(res), "ggplot")
  expect_s3_class(plot(res, type = "difference"), "ggplot")
})

# ---------------------------------------------------------------------------
# mxfda export
# ---------------------------------------------------------------------------

test_that("export_mxfda emits the documented columns", {
  res <- mk_result()
  mx <- export_mxfda(res, sample_id = "S1")
  expect_s3_class(mx, "data.frame")
  expect_identical(nrow(mx), length(res@h))
  expect_true(all(multiScaleAC:::.mxfda_columns %in% names(mx)))
  expect_identical(mx$r, res@h)
  expect_identical(mx$sumfun, res@observed)
  expect_equal(mx$fundiff, res@observed - res@perm_mean)
  expect_true(all(mx$sample_id == "S1"))
})

test_that("export_mxfda_batch binds results with and without adjusted p-values", {
  # rbind() used to fail here: p_adjusted was added as a column only when it existed, so a list where
  # some results had been adjusted and others had not produced frames with different column counts.
  # That is the likely state after lapply(results, compute_significance) without adjust_pvalues().
  a <- adjust_pvalues(mk_result(seed = 1))
  b <- mk_result(seed = 2)
  expect_length(a@p_adjusted, length(a@h))
  expect_length(b@p_adjusted, 0L)

  mx <- export_mxfda_batch(list(A = a, B = b))
  expect_s3_class(mx, "data.frame")
  expect_identical(nrow(mx), length(a@h) + length(b@h))
  expect_true("p_adjusted" %in% names(mx))
  # Labels come from the list names, since both results carry the default sample_id slot.
  expect_setequal(unique(mx$sample_id), c("A", "B"))
  expect_true(all(is.na(mx$p_adjusted[mx$sample_id == "B"])))
  expect_false(any(is.na(mx$p_adjusted[mx$sample_id == "A"])))
})

test_that("export_mxfda requires significance first", {
  set.seed(1); n <- 80
  raw <- multiScaleAC(matrix(runif(2 * n), ncol = 2), rnorm(n), rnorm(n),
                      h = 0.1, perms = 19, verbose = FALSE)
  expect_error(export_mxfda(raw), regexp = "compute_significance")
})

test_that("the export is accepted by mxfda", {
  # The claim that this output is "mxfda-compatible" was never executed. mxfda could not be installed
  # in the environment this was developed in, so this test skips there rather than asserting an
  # unverified compatibility. Run the suite with mxfda installed to actually check it.
  skip_if_not_installed("mxfda")
  res <- mk_result()
  mx <- export_mxfda(res)
  expect_true(all(c("r", "sumfun", "csr", "fundiff") %in% names(mx)))
  # Exercise the real entry point rather than only the column names.
  expect_error(
    mxfda::add_summary_function(mxfda::make_mxfda(
      metadata = data.frame(sample_id = unique(mx$sample_id)),
      spatial = NULL, subject_key = "sample_id", sample_key = "sample_id"
    ), mx, metric = "bivariate"),
    NA
  )
})
