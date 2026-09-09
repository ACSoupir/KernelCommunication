#' Multi-scale Spatial Autocorrelation Analysis
#'
#' Computes bivariate Moran's I across a range of kernel bandwidths, producing a
#' functional curve \eqn{I(h)} rather than a single value at one chosen scale.
#' Inference is by permutation of the second gene.
#'
#' @param coords numeric matrix of spatial coordinates (N x 2)
#' @param gene1 numeric vector of expression values for the first variable (length N)
#' @param gene2 numeric vector of expression values for the second variable (length N)
#' @param h numeric vector of kernel bandwidths, in the same units as
#'   \code{coords} (must be > 0). See the bandwidth section below; choosing these
#'   on the wrong scale is the most common way to get a meaningless result.
#' @param perms integer number of permutations (default: 1000). \code{0} is
#'   permitted and returns the observed curve with no permutation matrix.
#' @param seed integer random seed for reproducibility (default: 333)
#' @param eps numeric weight threshold (default: 0, meaning no thresholding).
#'   Weights are retained when \code{w > eps} and are excluded from the row sum
#'   otherwise. The default of 0 is deliberate: computing the statistic without
#'   discarding low weights is the point of the method, and is what distinguishes
#'   it from approaches that filter low kernel density estimates first. Raising
#'   \code{eps} is supported but changes the estimand.
#' @param nan_policy character, how to report a bandwidth at which some cell has
#'   no neighbours carrying weight. \code{"zero"} (default) reproduces the
#'   reference implementation, which reports 0 for the whole affected bandwidth;
#'   \code{"drop"} excludes only the isolated cells; \code{"na"} reports
#'   \code{NA} so the situation cannot pass unnoticed. All three warn.
#' @param moments logical, compute the exact analytical moments of the
#'   permutation null (default TRUE). These cost O(N) per bandwidth on top of work
#'   already done, so they are essentially free, and they are \emph{exact} rather
#'   than estimated: the statistic is a linear permutation statistic, whose
#'   moments have closed forms. They enable \code{method = "moment"} in
#'   \code{\link{compute_significance}} and populate \code{leverage}, which
#'   reports whether a continuous approximation to the null is trustworthy at each
#'   bandwidth.
#' @param permute character, how much permuting to do. \code{"all"} (default)
#'   permutes at every bandwidth. \code{"flagged"} permutes only where the
#'   moment diagnostics say a continuous approximation cannot be trusted, using
#'   the exact moments elsewhere; this is what turns permutations from the bulk
#'   cost into a targeted safety net. \code{"none"} skips them entirely and forces
#'   \code{perms = 0}.
#' @param batch_size integer permutations contracted per batch (default: 256).
#'   Bounds peak memory at N x \code{batch_size} instead of N x \code{perms};
#'   results are unaffected.
#' @param n_threads integer number of threads for the kernel sweep (default: 0,
#'   meaning all available). Ignored when the package was built without OpenMP,
#'   which includes builds using stock Apple clang.
#' @param verbose logical show progress (default: NULL = auto for large datasets)
#' @param large_mode_n deprecated and ignored. It previously claimed to select a
#'   sparse code path, but the threshold was hard-coded in C++ and this argument
#'   was never read. There is now a single compute path whose working memory is
#'   O(N x length(h)) regardless of N, so no threshold is needed.
#' @param kernel character or function. Built-in options are \code{"gaussian"}
#'   (the default, and the one the manuscript uses), \code{"uniform"} and
#'   \code{"exponential"}. If a function, a custom kernel with signature
#'   \code{kernel(dist_matrix, h, ...)} returning an N x N weight matrix.
#' @param distance character or function. Built-in options are
#'   \code{"euclidean"} (default), \code{"manhattan"} and \code{"minkowski"};
#'   the last takes its order from \code{p}. If a function, a custom distance
#'   with signature \code{distance(coords, ...)} returning an N x N distance
#'   matrix.
#' @param ... additional parameters forwarded to the distance and kernel
#'   functions, and \code{p} for the built-in Minkowski distance. Note that
#'   \code{...} precedes every option in the signature, so options must be named
#'   exactly; that is what stops a parameter such as \code{p} from being
#'   partially matched against one of them.
#'
#' @return An S4 object of class \linkS4class{kernelcomm_result} with slots:
#' \describe{
#'   \item{h}{numeric vector of bandwidth values used}
#'   \item{observed}{numeric vector of Moran's I for each bandwidth}
#'   \item{permutations}{numeric matrix (n_h x perms) of permuted I values}
#'   \item{coords}{matrix of original coordinates}
#'   \item{gene1_name}{character. Name of gene 1}
#'   \item{gene2_name}{character. Name of gene 2}
#'   \item{sample_id}{character. Sample identifier}
#'   \item{seed}{numeric. Random seed used}
#'   \item{n_cells}{integer. Number of cells}
#'   \item{call}{language. Function call that created this object}
#'   \item{p_values, z_scores, perm_mean, perm_sd, null_lower, null_upper, significant, p_adjusted}{NULL (computed by \code{\link{compute_significance}})}
#' }
#'
#' @details
#' \strong{Built-in Kernels:}
#' \itemize{
#'   \item \code{"gaussian"}: \eqn{w_{ij} = \exp(-d_{ij}^2 / (2h^2))}. The
#'     default, and the kernel the manuscript uses. Unbounded support, so every
#'     pair contributes.
#'   \item \code{"uniform"}: \eqn{w_{ij} = 1} if \eqn{d_{ij} \le h}, else 0.
#'     A hard radius, so cells can genuinely have no neighbours; see
#'     \code{nan_policy}.
#'   \item \code{"exponential"}: \eqn{w_{ij} = \exp(-d_{ij} / h)}. Heavier
#'     tailed than the Gaussian at the same bandwidth, so a given \code{h} reaches
#'     further.
#' }
#'
#' The diagonal is always zero: a cell is not its own neighbour, which for single
#' cells is not a meaningful notion. Row sums are normalised to 1, so \eqn{I(h)}
#' stays bounded and every cell contributes equally regardless of local density.
#'
#' \strong{Built-in Distance Methods:}
#' \itemize{
#'   \item \code{"euclidean"}: \eqn{\sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}}
#'   \item \code{"manhattan"}: \eqn{|x_i - x_j| + |y_i - y_j|}
#'   \item \code{"minkowski"}: \eqn{(|x_i - x_j|^p + |y_i - y_j|^p)^{1/p}},
#'     with \code{p} passed through \code{...}. \code{p = 1} and \code{p = 2}
#'     are routed to the Manhattan and Euclidean paths, which compute the same
#'     thing without calling \code{pow()}.
#' }
#'
#' All combinations run inside the same tiled core, so none of them can bypass the
#' \code{eps} predicate or the zero-row-sum policy. Only a user-supplied
#' \emph{function} takes the R fallback path.
#'
#' \strong{Custom Kernel Functions:}
#' Custom kernel functions must have signature: \code{kernel(dist_matrix, h, ...)} where:
#' \itemize{
#'   \item \code{dist_matrix}: N x N matrix of pairwise distances
#'   \item \code{h}: scalar bandwidth for current computation
#'   \item \code{...}: additional parameters passed from multiScaleAC
#' }
#' The function must return an N x N weight matrix (can be dense matrix or sparse Matrix::dsCMatrix).
#'
#' \strong{Custom Distance Functions:}
#' Custom distance functions must have signature: \code{distance(coords, ...)} where:
#' \itemize{
#'   \item \code{coords}: N x 2 matrix of coordinates
#'   \item \code{...}: additional parameters (e.g., \code{p} for Minkowski distance)
#' }
#' The function must return an N x N distance matrix.
#'
#' \strong{Choosing bandwidths:}
#' Bandwidths are in the units of \code{coords}. This is the most common source
#' of meaningless results. If \code{h} is far too large, every pairwise weight
#' approaches 1, the kernel becomes a constant matrix, and \eqn{I(h)} collapses
#' toward zero for every gene pair; if far too small, no cell has a neighbour
#' carrying weight and the statistic is undefined. The bundled
#' \code{\link{ccRCCExample}} data is in millimetres and spans about 1 mm, so
#' bandwidths of roughly 0.01 to 0.25 are appropriate; the same tissue measured
#' in micrometres would want 10 to 250. Use
#' \code{\link{suggest_bandwidths}} rather than guessing, and
#' \code{\link{qc_neighborhoods}} to confirm the effective neighbourhood size
#' before interpreting the curve.
#'
#' \strong{Performance:}
#' The built-in kernels and distances use a two-pass sweep that computes each
#' pair distance once and reuses it across all bandwidths, with working memory of
#' O(N x length(h)) rather than O(N^2). Custom kernel or distance functions are
#' evaluated in R and do build a full N x N matrix, so they are limited to
#' smaller N.
#'
#' @examples
#' # Bandwidths are in coordinate units. These coordinates span the unit square,
#' # so bandwidths of a few hundredths cover a handful of neighbours each.
#' set.seed(1)
#' n <- 200
#' coords <- matrix(runif(2 * n), ncol = 2)
#' g1 <- rnorm(n)
#' g2 <- rnorm(n)
#'
#' result <- multiScaleAC(coords, g1, g2, h = c(0.05, 0.1, 0.2), perms = 99)
#' result
#'
#' # Check that a bandwidth gives a sensible neighbourhood size before trusting
#' # the curve: a mean close to N means the kernel is saturated and I(h) will
#' # collapse toward zero regardless of the biology.
#' qc_neighborhoods(coords, h = 0.1)
#'
#' \dontrun{
#' # Custom kernel (exponential decay). Evaluated in R, so N x N memory.
#' exp_kernel <- function(dist_matrix, h) exp(-dist_matrix / h)
#' result_custom <- multiScaleAC(coords, g1, g2, h = c(0.05, 0.1),
#'                               kernel = exp_kernel)
#'
#' # Custom distance (Minkowski, p = 3) with the default Gaussian kernel.
#' minkowski_dist <- function(coords, p = 3) {
#'   as.matrix(dist(coords, method = "minkowski", p = p))
#' }
#' result_mink <- multiScaleAC(coords, g1, g2, h = c(0.05, 0.1),
#'                             distance = minkowski_dist, p = 3)
#' }
#'
#' @seealso \code{\link{compute_significance}} to add p-values,
#' \code{\link{suggest_bandwidths}} to choose \code{h},
#' \code{\link{multiScaleAC_batch}} for many samples,
#' \code{\link{soupir_moran_bv}} for single-bandwidth analysis, and
#' \code{\link{soupir_moran_multi}}, the reference implementation, which remains
#' supported and is not deprecated.
#'
#' @export
# NOTE ON THE SIGNATURE: `...` sits immediately after `h`, before every option.
#
# That placement is deliberate and load-bearing. `...` forwards parameters to user-supplied kernel and
# distance functions, and R only requires exact name matching for formals that appear AFTER `...`.
# With the options before it, a user passing `p = 3` for a Minkowski distance -- which is the example
# in this function's own documentation -- had that argument partially matched against `perms` and
# `permute` instead of being forwarded, and the call failed. Moving `...` forward makes every option
# exact-match-only, so no user parameter name can ever be captured by one.
#
# The consequence is that options must be named, which they should be anyway.
multiScaleAC <- function(coords,
                         gene1,
                         gene2,
                         h,
                         ...,
                         perms = 1000,
                         seed = 333,
                         eps = 0,
                         nan_policy = c("zero", "drop", "na"),
                         moments = TRUE,
                         permute = c("all", "flagged", "none"),
                         batch_size = 256,
                         n_threads = 0,
                         verbose = NULL,
                         large_mode_n = NULL,
                         kernel = c("gaussian", "uniform", "exponential"),
                         distance = c("euclidean", "manhattan", "minkowski")) {

  # Capture additional parameters
  extra_params <- list(...)

  nan_policy <- match.arg(nan_policy)
  permute <- match.arg(permute)
  if (permute == "none") perms <- 0L

  # Validate and normalize kernel parameter
  kernel_is_custom <- is.function(kernel)
  if (!kernel_is_custom) {
    kernel <- match.arg(kernel)
  }

  # Validate and normalize distance parameter
  distance_is_custom <- is.function(distance)
  if (!distance_is_custom) {
    distance <- match.arg(distance)
  }

  if (!is.null(large_mode_n)) {
    warning("`large_mode_n` is deprecated and ignored. It previously described a sparse code path ",
            "selected by a threshold hard-coded in C++, and the argument itself was never read. ",
            "There is now one compute path with working memory O(N * length(h)).", call. = FALSE)
  }

  # Validate inputs
  .validate_inputs(coords, gene1, gene2, h, perms,
                   eps, batch_size, seed,
                   kernel, distance, kernel_is_custom, distance_is_custom)

  N <- nrow(coords)

  # Sanitize h (filter non-positive values with warnings)
  h_clean <- .sanitize_h(h)

  if (is.null(verbose)) verbose <- N > 20000L

  # Center and scale genes. A zero-variance gene yields NaN throughout, which the reference maps to
  # an all-zero curve; that is detected here rather than allowed to propagate.
  g1_s <- as.vector(scale(gene1))
  g2_s <- as.vector(scale(gene2))
  degenerate <- .has_degenerate_gene(g1_s, g2_s)
  denom <- sum(g1_s^2)

  if (degenerate) {
    # Nothing to compute: reproduce the reference's all-zero (or NA) curve directly.
    res <- list(observed = rep(NA_real_, length(h_clean)),
                permutations = matrix(NA_real_, length(h_clean), max(0L, as.integer(perms))))
    res <- .apply_nan_policy(res, rep(0L, length(h_clean)), h_clean,
                             degenerate_input = TRUE, policy = nan_policy)
    result <- c(list(h = h_clean), res)
    return(.create_result_object(result, coords, gene1, gene2, seed, call = match.call()))
  }

  if (!kernel_is_custom && !distance_is_custom) {
    if (verbose) {
      message(sprintf("Kernel sweep: N = %d, %d bandwidths, %d thread(s)",
                      N, length(h_clean), if (n_threads > 0) n_threads else moran_max_threads()))
    }
    pass <- moran_kernel_pass(
      coords        = as.matrix(coords[, 1:2, drop = FALSE]),
      gene1         = g1_s,
      h_vec         = as.numeric(h_clean),
      eps           = eps,
      kernel        = kernel,
      dist_method   = distance,
      drop_isolated = TRUE,
      n_threads     = if (n_threads > 0) as.integer(n_threads) else moran_max_threads(),
      p             = if (!is.null(extra_params$p)) as.numeric(extra_params$p) else 2
    )
    V <- pass$V
    n_isolated <- pass$n_isolated
  } else {
    # Custom kernel/distance path (R-based). Returns V in the same n_h x N orientation.
    custom <- .compute_custom_kernel_moran(
      coords = coords, gene1 = g1_s, h_vec = h_clean, eps = eps,
      distance = distance, kernel = kernel,
      extra_params = extra_params, verbose = verbose
    )
    V <- custom$V
    n_isolated <- custom$n_isolated
  }

  # Exact moments of the permutation null. O(N) per bandwidth on top of V, which is already in hand,
  # so this is effectively free and is computed unless explicitly declined.
  mom <- if (isTRUE(moments)) .permutation_moments(V, g2_s, denom) else NULL

  # `permute` decides how much permuting to do. "flagged" is the point of the moment work: draw
  # permutations only where the diagnostics say a continuous approximation to the null cannot be
  # trusted, and rely on the exact moments elsewhere. V is still alive here, so a fallback never
  # rebuilds the kernel.
  keep_h <- switch(
    permute,
    all     = rep(TRUE, length(h_clean)),
    none    = rep(FALSE, length(h_clean)),
    flagged = if (is.null(mom)) rep(TRUE, length(h_clean)) else .flag_untrustworthy(mom)
  )
  if (permute == "flagged" && verbose) {
    message(sprintf("Permuting %d of %d bandwidths flagged as unsuitable for moment-based inference",
                    sum(keep_h), length(keep_h)))
  }

  n_perm_eff <- if (any(keep_h)) perms else 0L
  if (n_perm_eff > 0 && !all(keep_h)) {
    # Contract only the flagged rows, then place them back so the matrix stays n_h x perms with
    # NA marking bandwidths that were deliberately not permuted.
    sub <- .contract_permutations(V[keep_h, , drop = FALSE], g2_s, denom,
                                  n_perm_eff, seed, batch_size)
    res <- list(observed = as.vector(V %*% g2_s) / denom,
                permutations = matrix(NA_real_, length(h_clean), n_perm_eff))
    res$permutations[keep_h, ] <- sub$permutations
  } else {
    res <- .contract_permutations(V, g2_s, denom, n_perm_eff, seed, batch_size)
  }

  res <- .apply_nan_policy(res, n_isolated, h_clean,
                           degenerate_input = FALSE, policy = nan_policy)

  result <- list(h = h_clean, observed = res$observed, permutations = res$permutations)

  out <- .create_result_object(result, coords, gene1, gene2, seed, call = match.call())

  if (!is.null(mom)) {
    out@moment_var <- mom$var
    out@moment_skewness <- mom$skewness
    out@moment_excess_kurtosis <- mom$excess_kurtosis
    out@leverage <- mom$leverage
    out@moment_cov <- .permutation_covariance(mom)
  }
  methods::validObject(out)
  out
}


# Internal helper functions -------------------------------------------------

#' Validate input parameters
#'
#' @noRd
.validate_inputs <- function(coords, gene1, gene2, h, perms,
                             eps, batch_size, seed,
                             kernel, distance, kernel_is_custom, distance_is_custom) {

  # Data checks are delegated to the exported reporter so there is exactly one implementation. This
  # wrapper turns its report into errors and warnings. Previously the two were separate functions
  # with different rules, and the exported one was never called by anything.
  #
  # Zero variance is an error there but only a degenerate case here: the reference implementation
  # returns an all-zero curve rather than failing, and multiScaleAC() reproduces that under
  # nan_policy, so it is downgraded to a warning issued later by .apply_nan_policy().
  report <- validate_inputs(coords, gene1, gene2)
  hard <- grep("zero variance", report$errors, value = TRUE, invert = TRUE)
  if (length(hard) > 0) {
    stop(paste0("invalid input:\n", paste0("  - ", hard, collapse = "\n")), call. = FALSE)
  }
  for (w in report$warnings) warning(w, call. = FALSE)

  # Validate h. NA and Inf are rejected here rather than in .sanitize_h(), where `if (NA > 0)`
  # raised an uninformative error about a missing value.
  if (!is.numeric(h) || length(h) < 1) {
    stop("h must contain at least one bandwidth value")
  }
  if (anyNA(h)) {
    stop("h contains missing values; every bandwidth must be a finite positive number")
  }
  if (any(is.infinite(h))) {
    stop("h contains infinite values; every bandwidth must be a finite positive number")
  }
  if (anyDuplicated(h)) {
    warning("h contains duplicated bandwidths; they will each be computed and reported separately",
            call. = FALSE)
  }

  # Validate perms. Zero is allowed: it yields the observed curve with no permutation matrix, which
  # is the mode the analytical-moment inference uses.
  if (!is.numeric(perms) || length(perms) != 1 || is.na(perms) ||
      perms < 0 || perms != as.integer(perms)) {
    stop("perms must be a single non-negative integer")
  }

  # Validate seed
  if (!is.numeric(seed)) {
    stop("seed must be an integer")
  }

  # Validate eps
  if (!is.numeric(eps) || eps < 0) {
    stop("eps must be a non-negative numeric value")
  }

  # Validate batch_size
  if (!is.numeric(batch_size) || batch_size < 1 || batch_size != as.integer(batch_size)) {
    stop("batch_size must be a positive integer")
  }

  # Validate kernel
  if (!kernel_is_custom) {
    if (!kernel %in% c("gaussian", "uniform", "exponential")) {
      stop("kernel must be 'gaussian', 'uniform', 'exponential', or a function")
    }
  } else {
    # kernel is already verified to be a function
  }

  # Validate distance
  if (!distance_is_custom) {
    if (!distance %in% c("euclidean", "manhattan", "minkowski")) {
      stop("distance must be 'euclidean', 'manhattan', 'minkowski', or a function")
    }
  } else {
    # distance is already verified to be a function
  }
}


#' Sanitize bandwidth vector (filter non-positive values)
#'
#' NA and Inf are rejected earlier, by .validate_inputs(). This function previously evaluated
#' `sum(h < 0)` and then `if (NA > 0)`, which failed with "missing value where TRUE/FALSE needed"
#' rather than saying what was wrong.
#'
#' @noRd
.sanitize_h <- function(h) {
  h_clean <- h

  # Remove negative values
  neg_count <- sum(h < 0)
  if (neg_count > 0) {
    warning(paste0("Removed ", neg_count, " negative bandwidth values"))
    h_clean <- h_clean[h_clean >= 0]
  }

  # Remove zero values
  zero_count <- sum(h_clean == 0)
  if (zero_count > 0) {
    warning(paste0("Removed ", zero_count, " zero bandwidth values (undefined kernel)"))
    h_clean <- h_clean[h_clean > 0]
  }

  # Check at least one valid
  if (length(h_clean) < 1) {
    stop("No valid bandwidth values remaining")
  }

  return(h_clean)
}


#' Compute spatially lagged gene1 with a custom kernel or distance (R-based)
#'
#' Returns \code{V} in the same n_h x N orientation as \code{moran_kernel_pass()}, so the caller
#' contracts it against gene2 and its permutations through exactly the same code. Previously this
#' function duplicated the contraction, the row normalisation and the zero-row-sum handling, and
#' each copy behaved differently from the C++ paths.
#'
#' Three fixes over the previous version:
#'
#' \itemize{
#'   \item It no longer stops with "Built-in kernels should use C++ path" when only the distance is
#'     custom. The path is entered whenever \emph{either} the kernel or the distance is a function,
#'     so a custom distance with the default Gaussian kernel used to be unreachable -- which made
#'     the Minkowski example in this function's own documentation throw an error.
#'   \item The diagonal is now zeroed. It was not, so a custom kernel returning \code{exp(0) = 1} on
#'     the diagonal gave every cell a self-weight, putting custom results on a different scale from
#'     built-in ones and quietly making them incomparable.
#'   \item Row sums of zero are reported rather than silently replaced by 1, so the caller's
#'     nan_policy applies here too.
#' }
#'
#' This path builds a full N x N matrix and is therefore limited to smaller N, which is documented
#' rather than worked around: an arbitrary R function cannot be evaluated tile-by-tile.
#'
#' @noRd
.compute_custom_kernel_moran <- function(coords, gene1, h_vec, eps,
                                         distance, kernel, extra_params,
                                         verbose) {
  p_order <- if (!is.null(extra_params$p)) as.numeric(extra_params$p) else 2

  N <- nrow(coords)
  n_h <- length(h_vec)

  # Distance matrix
  if (is.function(distance)) {
    dist_formals <- names(formals(distance))
    dist_params <- extra_params[names(extra_params) %in% dist_formals]
    dist_matrix <- do.call(distance, c(list(coords), dist_params))
    dist_matrix <- as.matrix(dist_matrix)
  } else {
    dist_matrix <- if (identical(distance, "minkowski")) {
      as.matrix(stats::dist(coords, method = "minkowski", p = p_order))
    } else {
      as.matrix(stats::dist(coords, method = distance))
    }
  }

  if (!is.matrix(dist_matrix) || nrow(dist_matrix) != N || ncol(dist_matrix) != N) {
    stop("the distance function must return a square N x N matrix")
  }
  if (anyNA(dist_matrix)) {
    stop("the distance function returned missing values")
  }

  V <- matrix(0, nrow = n_h, ncol = N)
  n_isolated <- integer(n_h)

  for (h_idx in seq_len(n_h)) {
    h <- h_vec[h_idx]

    if (is.function(kernel)) {
      kernel_formals <- names(formals(kernel))
      kernel_params <- extra_params[names(extra_params) %in% kernel_formals]
      W <- do.call(kernel, c(list(dist_matrix, h), kernel_params))
      W <- as.matrix(W)
      if (!identical(dim(W), c(N, N))) {
        stop("the kernel function must return a square N x N weight matrix")
      }
    } else if (kernel == "gaussian") {
      W <- exp(-dist_matrix^2 / (2 * h^2))
    } else if (kernel == "uniform") {
      W <- (dist_matrix <= h) * 1.0
    } else if (kernel == "exponential") {
      W <- exp(-dist_matrix / h)
    } else {
      stop("unsupported kernel: ", kernel)
    }

    # No self-weight, matching the built-in path and the reference implementation.
    diag(W) <- 0

    # Same eps predicate as the compute core: keep w > eps, and exclude the rest from the row sum.
    if (eps > 0) W[W <= eps] <- 0

    row_sums <- rowSums(W)
    isolated <- row_sums <= 0
    n_isolated[h_idx] <- sum(isolated)

    u <- gene1 / row_sums
    u[isolated] <- 0        # contribute nothing; the caller's nan_policy decides what to report

    V[h_idx, ] <- as.vector(W %*% u)

    if (verbose && n_h > 1) {
      message(sprintf("  bandwidth %d/%d (h = %s)", h_idx, n_h, format(h, digits = 4)))
    }
  }

  list(V = V, n_isolated = n_isolated)
}


# The analytical variance formula that used to sit here commented out --
# Var(I) = ||t(W) x||^2 / (N-1)^2 -- is correct, and now lives in R/moments.R as the special case of
# the exact second moment when both genes are z-scored. See docs/RETRACTIONS.md section 4 for why the
# test report that accompanied it was retracted even though the mathematics was sound.


#' Create kernelcomm_result S4 object
#'
#' Gene names are taken from the recorded call, not from \code{substitute()}.
#' \code{deparse(substitute(gene1))} evaluated in this frame always returned the literal string
#' \code{"gene1"}, because the promise here is the symbol from \code{multiScaleAC()}'s body rather
#' than the caller's expression, and the \code{tryCatch} fallback could therefore never fire. Every
#' plot title and every mxfda interaction name was consequently "gene1 x gene2".
#'
#' @param result List with h, observed, permutations from computation
#' @param coords Matrix of coordinates
#' @param gene1 Gene 1 expression vector
#' @param gene2 Gene 2 expression vector
#' @param seed Random seed used
#' @param call Function call that created this object, from \code{match.call()}
#'
#' @return kernelcomm_result S4 object
#' @noRd
.create_result_object <- function(result, coords, gene1, gene2, seed, call) {
  # match.call() names its arguments, so call$gene1 is the caller's expression regardless of whether
  # they passed positionally. Long expressions are truncated so they stay usable as plot labels.
  arg_label <- function(expr, fallback) {
    if (is.null(expr)) return(fallback)
    txt <- paste(deparse(expr), collapse = "")
    txt <- gsub("\\s+", " ", txt)
    if (nchar(txt) > 40L) txt <- paste0(substr(txt, 1L, 37L), "...")
    if (!nzchar(txt)) fallback else txt
  }

  create_kernelcomm_result(
    h = result$h,
    observed = result$observed,
    permutations = result$permutations,
    coords = as.matrix(coords),
    gene1_name = arg_label(call$gene1, "gene1"),
    gene2_name = arg_label(call$gene2, "gene2"),
    sample_id = "sample1",  # Default, can be overridden by multiScaleAC_batch()
    seed = seed,
    n_cells = nrow(coords),
    call = call
  )
}


#' Batch Processing for Multiple Samples
#'
#' Runs \code{\link{multiScaleAC}} over many samples, sequentially or in
#' parallel, with resumable checkpointing.
#'
#' @section Checkpointing:
#' Each sample is written to its own file in \code{checkpoint_dir}, named after its
#' sample id. Resuming means skipping the samples whose files already exist, so an
#' interrupted batch continues where it stopped and nothing is recomputed.
#'
#' Per-sample files replace an earlier design that kept one shared list and one
#' shared file. That design lost every result under parallel processing, because
#' \code{mclapply} forks and the worker's assignment into the parent's list was
#' discarded, and all workers raced to write the same path.
#'
#' Files are removed only when every sample has succeeded, so a batch with
#' failures stays resumable.
#'
#' @section Multiple testing across samples and pairs:
#' The manuscript applies its FDR threshold across gene pairs. To reproduce that,
#' collect one p-value per pair from these results -- \code{sup_p_value} is the
#' natural choice, since it summarises the whole curve -- and correct across the
#' collection. Do not use \code{\link{adjust_pvalues}} for this: it corrects
#' across the bandwidths within one result, which are strongly dependent.
#'
#' @param data_list List of data frames, each containing coordinates and gene
#'   expression. Columns are checked for every sample before any computation
#'   starts, so a missing column fails immediately rather than hours in.
#' @param sample_ids Character vector of sample identifiers (optional,
#'   auto-generated from names if NULL). Must be unique, since checkpoint files
#'   are keyed by them.
#' @param h_vec Numeric vector of bandwidths, in the units of the coordinate
#'   columns. See \code{\link{multiScaleAC}} on choosing these.
#' @param ... Additional arguments passed to \code{\link{multiScaleAC}}, for
#'   example \code{kernel}, \code{permute} or \code{nan_policy}.
#' @param perms Integer number of permutations per sample.
#' @param seed Integer base seed. Sample \code{i} uses \code{seed + i}, so
#'   samples are not identically permuted while the batch stays reproducible.
#' @param n_cores Integer cores for parallel processing (default 1 = sequential).
#'   Parallelism uses forking and is therefore unavailable on Windows, where this
#'   falls back to sequential with a warning.
#' @param checkpoint_dir Character directory for per-sample checkpoint files
#'   (default NULL = no checkpointing).
#' @param coords_cols Character vector of two coordinate column names.
#' @param gene1_col,gene2_col Character column names for the two genes. These also
#'   become the gene names on each result.
#' @param verbose Logical, report progress.
#'
#' @return Named list of \linkS4class{kernelcomm_result} objects, one per sample.
#'   A sample that failed is \code{NULL} and is reported by a warning.
#' @export
#'
#' @examples
#' # Split the bundled sample into pseudo-samples by field of view.
#' data(ccRCCExample)
#' fovs <- split(ccRCCExample, ccRCCExample$fov)
#' fovs <- fovs[vapply(fovs, nrow, integer(1)) > 100][1:3]
#'
#' # Bandwidths in millimetres, matching the coordinate units of this dataset.
#' results <- multiScaleAC_batch(
#'   fovs, h_vec = c(0.02, 0.05, 0.1), perms = 99, verbose = FALSE
#' )
#' vapply(results, function(r) r@observed[1], numeric(1))
#'
#' \dontrun{
#' # In parallel, with a resumable checkpoint directory.
#' results <- multiScaleAC_batch(
#'   fovs, h_vec = c(0.02, 0.05, 0.1), perms = 999,
#'   n_cores = 4, checkpoint_dir = tempfile("msac_ckpt_")
#' )
#' }
#'
#' @seealso \code{\link{export_mxfda_batch}} for exporting batch results
multiScaleAC_batch <- function(data_list,
                                sample_ids = NULL,
                                h_vec,
                                ...,
                                perms = 1000,
                                seed = 333,
                                n_cores = 1,
                                checkpoint_dir = NULL,
                                coords_cols = c("x", "y"),
                                gene1_col = "CD274",
                                gene2_col = "PDCD1",
                                verbose = TRUE) {

  if (!is.list(data_list) || length(data_list) == 0) {
    stop("data_list must be a non-empty list")
  }
  n_samples <- length(data_list)

  if (is.null(sample_ids)) {
    sample_ids <- names(data_list)
    if (is.null(sample_ids) || all(!nzchar(sample_ids))) {
      sample_ids <- paste0("sample", seq_len(n_samples))
    }
  } else if (length(sample_ids) != n_samples) {
    stop("sample_ids length must match data_list length")
  }
  if (anyDuplicated(sample_ids)) {
    stop("sample_ids must be unique; checkpoint files are keyed by them")
  }

  say <- function(...) if (verbose) message(...)

  # Validate every sample's columns up front, before any computation. Previously the columns were
  # subset before being checked, so a missing column produced "subscript out of bounds" and the three
  # carefully written error messages below were unreachable -- and it failed only once that sample was
  # reached, potentially hours into a batch.
  required <- c(coords_cols, gene1_col, gene2_col)
  for (i in seq_len(n_samples)) {
    have <- colnames(data_list[[i]])
    missing <- setdiff(required, have)
    if (length(missing)) {
      stop("sample ", sample_ids[i], " is missing column(s): ", paste(missing, collapse = ", "),
           ". Available: ", paste(have, collapse = ", "))
    }
  }

  # Checkpointing.
  #
  # The previous design kept one list and did `checkpoint[[i]] <<- result` inside the worker, then
  # wrote the whole list to a single path. Under parallel::mclapply the workers are forked, so those
  # superassignments were discarded, and every child raced to write the same file, leaving whichever
  # finished last. It also resumed with `start_idx:n_samples`, which counts DOWN when start_idx is
  # n_samples + 1, silently reprocessing the entire batch in reverse.
  #
  # Now each sample owns one file, named by its sample_id. Nothing is shared, so forks are safe, and
  # resuming means simply skipping the samples whose files already exist.
  use_ckpt <- !is.null(checkpoint_dir)
  if (use_ckpt) {
    if (!is.character(checkpoint_dir) || length(checkpoint_dir) != 1) {
      stop("checkpoint_dir must be a single directory path")
    }
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }
  ckpt_path <- function(i) {
    file.path(checkpoint_dir, paste0("msac_", .sanitize_filename(sample_ids[i]), ".rds"))
  }

  done <- rep(FALSE, n_samples)
  if (use_ckpt) {
    done <- vapply(seq_len(n_samples), function(i) file.exists(ckpt_path(i)), logical(1))
    if (any(done)) {
      say("Resuming: ", sum(done), " of ", n_samples, " samples already complete")
    }
  }

  process_sample <- function(i) {
    if (done[i]) {
      return(tryCatch(readRDS(ckpt_path(i)), error = function(e) {
        warning("checkpoint for ", sample_ids[i], " is unreadable; recomputing", call. = FALSE)
        NULL
      }))
    }
    say("Processing ", i, "/", n_samples, ": ", sample_ids[i])
    data <- data_list[[i]]

    result <- multiScaleAC(
      coords = as.matrix(data[, coords_cols, drop = FALSE]),
      gene1 = data[[gene1_col]],
      gene2 = data[[gene2_col]],
      h = h_vec,
      # A distinct seed per sample, so samples are not identically permuted, while the batch as a
      # whole stays reproducible from `seed`.
      perms = perms, seed = seed + i, verbose = FALSE,
      ...
    )
    result@sample_id <- sample_ids[i]
    result@gene1_name <- gene1_col
    result@gene2_name <- gene2_col

    if (use_ckpt) saveRDS(result, ckpt_path(i))
    result
  }

  todo <- which(!done)
  parallel_ok <- n_cores > 1 && .Platform$OS.type == "unix"
  if (n_cores > 1 && !parallel_ok) {
    warning("parallel processing uses forking, which is unavailable on this platform; ",
            "running sequentially", call. = FALSE)
  }

  results <- vector("list", n_samples)

  # Recover any completed samples first, whether or not anything remains to do.
  for (i in which(done)) results[[i]] <- process_sample(i)

  if (length(todo)) {
    if (parallel_ok) {
      say("Processing ", length(todo), " samples on ", n_cores, " cores")
      got <- parallel::mclapply(todo, process_sample,
                                mc.cores = min(n_cores, length(todo)),
                                mc.preschedule = TRUE)
      # mclapply reports a failed fork as a try-error rather than raising, so failures must be
      # inspected instead of assumed away.
      for (j in seq_along(todo)) {
        if (inherits(got[[j]], "try-error") || is.null(got[[j]])) {
          warning("sample ", sample_ids[todo[j]], " failed: ",
                  as.character(got[[j]]), call. = FALSE)
        } else {
          results[[todo[j]]] <- got[[j]]
        }
      }
    } else {
      # seq_along, never start_idx:n_samples, which counts backwards when nothing is left to do.
      for (i in todo) results[[i]] <- process_sample(i)
    }
  } else {
    say("Nothing to do: every sample was already complete")
  }

  failed <- vapply(results, is.null, logical(1))
  if (any(failed)) {
    warning(sum(failed), " sample(s) produced no result: ",
            paste(sample_ids[failed], collapse = ", "),
            ". Checkpoint files are kept so the batch can be resumed.", call. = FALSE)
  } else if (use_ckpt) {
    unlink(vapply(seq_len(n_samples), ckpt_path, character(1)))
    say("Checkpoint files removed (batch complete)")
  }

  names(results) <- sample_ids
  say("Batch complete: ", sum(!failed), "/", n_samples, " samples")
  results
}


#' Make a string safe to use in a file name
#'
#' @noRd
.sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  substr(x, 1, 100)
}
