# Shared assembly of a multi-bandwidth result from the kernel pass.
#
# One place owns the zero-row-sum policy and the permutation contraction, so the three-way
# divergence that used to exist between the dense C++ path, the sparse C++ path and the R
# custom-kernel path cannot recur.

#' Contract a kernel pass against gene2 and its permutations
#'
#' Internal. Takes \code{V} (n_h x N) from \code{moran_kernel_pass()} and produces the observed
#' curve and the permutation matrix.
#'
#' Permutations are drawn and contracted in batches. The reference implementation, and the previous
#' version of this package, materialised the full N x perms permuted-expression matrix up front:
#' 57 MB at N = 7122 with 999 permutations, and 800 MB at N = 100000. Because the contraction is a
#' matrix product against a fixed \code{V}, it can be done a batch at a time, which bounds that
#' memory at N x batch_size while producing bit-identical results, since the permutations are drawn
#' from the same seeded stream in the same order.
#'
#' @param V n_h x N matrix of spatially lagged gene1, one row per bandwidth.
#' @param gene2_scaled Length-N centred and scaled gene2.
#' @param denom \code{sum(gene1_scaled^2)}.
#' @param perms Number of permutations; 0 is permitted and yields a zero-column matrix.
#' @param seed Random seed. The permutation stream is drawn exactly as the reference draws it.
#' @param batch_size Permutations contracted per batch.
#'
#' @return A list with \code{observed} (length n_h) and \code{permutations} (n_h x perms).
#' @keywords internal
#' @noRd
.contract_permutations <- function(V, gene2_scaled, denom, perms, seed,
                                   batch_size = 256L) {
  n_h <- nrow(V)
  N   <- ncol(V)

  observed <- as.vector(V %*% gene2_scaled) / denom

  if (perms < 1) {
    return(list(observed = observed,
                permutations = matrix(numeric(0), nrow = n_h, ncol = 0)))
  }

  batch_size <- max(1L, min(as.integer(batch_size), as.integer(perms)))
  out <- matrix(NA_real_, nrow = n_h, ncol = perms)

  # Seed once, then draw the whole stream in order. Batching changes only how many draws are held
  # in memory at a time, never the sequence, so results match an unbatched run exactly.
  set.seed(seed)
  done <- 0L
  while (done < perms) {
    take <- min(batch_size, perms - done)
    Yb <- matrix(NA_real_, nrow = N, ncol = take)
    for (b in seq_len(take)) Yb[, b] <- sample(gene2_scaled, N, replace = FALSE)
    out[, (done + 1L):(done + take)] <- (V %*% Yb) / denom
    done <- done + take
  }

  list(observed = observed, permutations = out)
}

#' Apply the zero-row-sum policy
#'
#' Internal. Decides what a bandwidth reports when at least one cell has no neighbours carrying
#' weight, or when an input gene has no variance.
#'
#' The reference implementation divides gene1 by a zero row sum, so \code{v_h} picks up
#' \code{NaN}; that propagates through the inner product, and \code{observed[is.nan(observed)] <- 0}
#' then zeroes the \emph{entire bandwidth} rather than only the offending cell. The same happens for
#' a zero-variance gene, where \code{scale()} returns \code{NaN} throughout.
#'
#' That is surprising, and arguably not what a user wants, but it is what the reference does, and the
#' reference is the oracle this package is validated against. Silently improving on it would break
#' every fixture while leaving the suite green. So \code{"zero"} reproduces it exactly and warns;
#' the alternatives are offered explicitly.
#'
#' @param res List with \code{observed} and \code{permutations}.
#' @param n_isolated Integer vector, cells with zero row sum per bandwidth.
#' @param h Bandwidth vector, for the warning message.
#' @param degenerate_input TRUE when a gene had no variance, making every bandwidth degenerate.
#' @param policy One of "zero", "drop", "na".
#'
#' @return \code{res}, with the policy applied.
#' @keywords internal
#' @noRd
.apply_nan_policy <- function(res, n_isolated, h, degenerate_input, policy = "zero") {
  policy <- match.arg(policy, c("zero", "drop", "na"))

  affected <- if (degenerate_input) rep(TRUE, length(h)) else n_isolated > 0L
  # Any non-finite value that reached this point is degenerate regardless of the census.
  affected <- affected | !is.finite(res$observed)

  if (!any(affected)) return(res)

  if (degenerate_input) {
    warning("A supplied gene has zero variance, so Moran's I is undefined at every bandwidth. ",
            "Reporting ", switch(policy, zero = "0", drop = "0", na = "NA"),
            " (nan_policy = \"", policy, "\").", call. = FALSE)
  } else {
    warning(sum(affected), " of ", length(h), " bandwidths have at least one cell with no ",
            "neighbours carrying weight (bandwidths: ",
            paste(format(h[affected], digits = 3), collapse = ", "), "; up to ",
            max(n_isolated), " isolated cells). ",
            "The reference implementation reports 0 for the whole bandwidth in this case, which is ",
            "what nan_policy = \"zero\" reproduces. Consider larger bandwidths, or ",
            "nan_policy = \"na\" to make it explicit.", call. = FALSE)
  }

  fill <- switch(policy, zero = 0, drop = 0, na = NA_real_)

  if (policy == "drop") {
    # "drop" already took effect inside the kernel pass, where an isolated cell's contribution was
    # zeroed instead of divided. Whatever remains non-finite here is genuinely undefined.
    bad <- !is.finite(res$observed)
    res$observed[bad] <- 0
    if (ncol(res$permutations)) res$permutations[bad, ] <- 0
    return(res)
  }

  res$observed[affected] <- fill
  if (ncol(res$permutations)) res$permutations[affected, ] <- fill
  res
}

#' Detect an input that makes the statistic undefined everywhere
#'
#' Internal. \code{scale()} of a zero-variance vector is \code{NaN} throughout, and the reference
#' then returns 0 for every bandwidth.
#'
#' @param ... Scaled gene vectors.
#' @return TRUE when any is not entirely finite.
#' @keywords internal
#' @noRd
.has_degenerate_gene <- function(...) {
  any(vapply(list(...), function(g) !all(is.finite(g)), logical(1)))
}
