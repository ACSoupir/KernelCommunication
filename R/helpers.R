#' Validate Input Data
#'
#' Checks for common problems before computation and returns a report, rather
#' than raising an error. Use this to inspect data interactively;
#' \code{\link{multiScaleAC}} runs the same checks itself and stops on errors.
#'
#' Previously this function was exported and documented but never called by
#' anything, while \code{multiScaleAC()} used a separate internal validator with
#' different rules. There is now one implementation with two surfaces: this one
#' reports, and the internal wrapper enforces.
#'
#' @param coords Coordinates matrix (N x 2)
#' @param gene1 Gene 1 expression vector
#' @param gene2 Gene 2 expression vector
#'
#' @return List with components:
#'   \item{errors}{Character vector of problems that prevent computation}
#'   \item{warnings}{Character vector of things worth knowing about}
#'   \item{valid}{Logical. TRUE if there are no errors}
#'
#' @seealso \code{\link{qc_neighborhoods}} to check whether a bandwidth is
#'   sensible, and \code{\link{suggest_bandwidths}} to choose one.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' coords <- matrix(runif(200), ncol = 2)
#' report <- validate_inputs(coords, rnorm(100), rnorm(100))
#' report$valid
#'
#' # A constant gene is reported rather than silently producing zeros
#' validate_inputs(coords, rep(1, 100), rnorm(100))$errors
validate_inputs <- function(coords, gene1, gene2) {
  issues <- list(errors = character(), warnings = character())

  add_err <- function(...) issues$errors <<- c(issues$errors, paste0(...))
  add_wrn <- function(...) issues$warnings <<- c(issues$warnings, paste0(...))

  # Structural checks first, and return before anything that would depend on them. The previous
  # version called nrow(coords) before this check had short-circuited, so passing a plain vector
  # produced "argument is of length zero" instead of a useful message.
  if (!is.matrix(coords) || ncol(coords) < 2) {
    add_err("coords must be a matrix with at least 2 columns")
  } else if (!is.numeric(coords)) {
    add_err("coords must be numeric")
  }
  if (!is.numeric(gene1)) add_err("gene1 must be a numeric vector")
  if (!is.numeric(gene2)) add_err("gene2 must be a numeric vector")

  if (length(issues$errors) > 0) {
    issues$valid <- FALSE
    return(issues)
  }

  N <- nrow(coords)
  if (length(gene1) != N) add_err("gene1 length (", length(gene1), ") must match nrow(coords) (", N, ")")
  if (length(gene2) != N) add_err("gene2 length (", length(gene2), ") must match nrow(coords) (", N, ")")
  if (N < 3) add_err("at least 3 observations are required; got ", N)

  if (length(issues$errors) > 0) {
    issues$valid <- FALSE
    return(issues)
  }

  # Missing and non-finite values next. These were previously checked last, after mean(gene1 == 0)
  # and sd(gene1) had already run on them, so NA input died with "missing value where TRUE/FALSE
  # needed" rather than producing the message written for it.
  if (anyNA(coords))                add_err("coords contains missing values")
  if (any(!is.finite(coords)))      add_err("coords contains infinite values")
  if (anyNA(gene1) || anyNA(gene2)) add_err("gene expression contains missing values")
  if (any(!is.finite(gene1[!is.na(gene1)])) || any(!is.finite(gene2[!is.na(gene2)]))) {
    add_err("gene expression contains infinite values")
  }

  if (length(issues$errors) > 0) {
    issues$valid <- FALSE
    return(issues)
  }

  # Now it is safe to compute summaries.
  if (stats::sd(gene1) == 0) add_err("gene1 has zero variance, so Moran's I is undefined")
  if (stats::sd(gene2) == 0) add_err("gene2 has zero variance, so Moran's I is undefined")

  if (anyDuplicated(coords)) {
    add_wrn("found ", sum(duplicated(coords)), " duplicated coordinates; consider adding a small ",
            "jitter, since coincident cells receive the maximum kernel weight from each other")
  }

  zf1 <- mean(gene1 == 0)
  zf2 <- mean(gene2 == 0)
  if (zf1 > 0.9) {
    add_wrn("gene1 is ", round(zf1 * 100), "% zeros; the permutation null will be strongly skewed, ",
            "so treat a normal approximation to it with suspicion")
  }
  if (zf2 > 0.9) {
    add_wrn("gene2 is ", round(zf2 * 100), "% zeros; the permutation null will be strongly skewed, ",
            "so treat a normal approximation to it with suspicion")
  }

  issues$valid <- length(issues$errors) == 0
  issues
}


#' Suggest Bandwidth Values
#'
#' Proposes bandwidths spanning a range of spatial scales, in the units of
#' \code{coords}.
#'
#' \code{"local"} bandwidths are derived from nearest-neighbour distances, so
#' each cell sees a handful of neighbours. \code{"global"} bandwidths are derived
#' from the spatial extent of the sample, so each cell sees a substantial
#' fraction of it. \code{"medium"} spans between the two, which is usually what
#' you want for a functional curve.
#'
#' The previous version returned quantiles of the \emph{same} nearest-neighbour
#' distribution for all three settings, so \code{"global"} was not global: every
#' option produced near-nearest-neighbour bandwidths. It also built a full
#' N x N distance matrix, which is about 400 MB for the bundled 7122-cell example
#' and 80 GB at N = 100000. Nearest-neighbour distances are now computed in one
#' pass with O(N) memory.
#'
#' @param coords Coordinate matrix (N x 2)
#' @param n_bandwidths Number of bandwidths to return
#' @param range Character. One of "local", "medium", "global".
#' @param distance Character. "euclidean" or "manhattan", matching what you will
#'   pass to \code{\link{multiScaleAC}}.
#'
#' @return Numeric vector of suggested bandwidths, in coordinate units,
#'   increasing. Carries a \code{"scale_reference"} attribute recording the
#'   median nearest-neighbour distance and the extent used to derive them.
#'
#' @seealso \code{\link{qc_neighborhoods}} to confirm the resulting
#'   neighbourhood sizes.
#'
#' @export
#'
#' @examples
#' data(ccRCCExample)
#' coords <- as.matrix(ccRCCExample[, c("x", "y")])
#'
#' # Coordinates are in millimetres, so these come back in millimetres.
#' h <- suggest_bandwidths(coords, n_bandwidths = 6)
#' round(h, 4)
#' attr(h, "scale_reference")
suggest_bandwidths <- function(coords, n_bandwidths = 5,
                               range = c("medium", "local", "global"),
                               distance = c("euclidean", "manhattan")) {

  if (!is.matrix(coords) || ncol(coords) < 2) {
    stop("coords must be a matrix with at least 2 columns")
  }
  if (!is.numeric(n_bandwidths) || length(n_bandwidths) != 1 || n_bandwidths < 1) {
    stop("n_bandwidths must be a single positive integer")
  }
  # `range` shadows base::range inside this function, so the argument is resolved immediately and
  # the name is not used again.
  span_choice <- match.arg(range)
  distance <- match.arg(distance)
  n_bandwidths <- as.integer(n_bandwidths)

  xy <- as.matrix(coords[, 1:2, drop = FALSE])
  nn <- nn_distance(xy, dist_method = distance, n_threads = moran_max_threads())
  nn_med <- stats::median(nn[is.finite(nn) & nn > 0])

  # Spatial extent, as the typical distance across the sample rather than the maximum, which a
  # single outlying cell would dominate.
  extent <- if (distance == "euclidean") {
    sqrt(diff(base::range(xy[, 1]))^2 + diff(base::range(xy[, 2]))^2)
  } else {
    diff(base::range(xy[, 1])) + diff(base::range(xy[, 2]))
  }

  # A Gaussian bandwidth of h gives an effective neighbourhood radius of roughly h, so anchoring
  # "local" near the nearest-neighbour distance and "global" near a fraction of the extent gives
  # neighbourhood sizes that actually differ by orders of magnitude.
  lo <- switch(span_choice, local = nn_med,      medium = nn_med,          global = extent / 20)
  hi <- switch(span_choice, local = 4 * nn_med,  medium = extent / 8,      global = extent / 3)

  if (!is.finite(lo) || !is.finite(hi) || lo <= 0 || hi <= lo) {
    stop("could not derive a bandwidth range from these coordinates; are they degenerate?")
  }

  # Geometric spacing, because the statistic varies with the neighbourhood area rather than with h
  # linearly, and because a log-spaced grid is what a functional curve wants.
  h <- if (n_bandwidths == 1L) hi else exp(seq(log(lo), log(hi), length.out = n_bandwidths))

  attr(h, "scale_reference") <- c(median_nn_distance = nn_med, extent = extent)
  h
}


#' Compute Neighborhood Quality Metrics
#'
#' Reports the effective neighbourhood size at one or more bandwidths, so you can
#' tell whether a bandwidth is usable before interpreting a curve computed with
#' it.
#'
#' \strong{What "effective" means.} For the Gaussian kernel every pair has a
#' non-zero weight, so counting neighbours within a hard radius describes a
#' different model from the one being fitted. This function reports the kernel row
#' sum, which is the weight actually carried, and equals a plain neighbour count
#' only for the uniform kernel. The previous version applied a hard
#' \code{distance <= h} cutoff regardless of kernel, so under the default
#' Gaussian setting it could never report an isolated cell.
#'
#' Interpretation: a mean row sum close to \code{N} means the kernel is saturated,
#' every cell is effectively everyone's neighbour, and Moran's I will collapse
#' toward zero whatever the biology. A mean near zero means almost no cell has a
#' neighbour carrying weight and the statistic is unstable or undefined.
#'
#' @param coords Coordinates matrix (N x 2)
#' @param h Bandwidth, or a vector of bandwidths. A vector previously recycled
#'   silently against the distance matrix and returned meaningless numbers.
#' @param kernel Character. "gaussian" or "uniform", matching what you will pass
#'   to \code{\link{multiScaleAC}}.
#' @param distance Character. "euclidean" or "manhattan".
#' @param eps Weight threshold, matching \code{\link{multiScaleAC}}.
#'
#' @return A data frame with one row per bandwidth:
#'   \item{h}{The bandwidth}
#'   \item{mean_neighbors}{Mean kernel row sum, the effective neighbourhood size}
#'   \item{sd_neighbors}{Standard deviation of the row sums}
#'   \item{min_neighbors}{Smallest row sum}
#'   \item{max_neighbors}{Largest row sum}
#'   \item{isolated_cells}{Cells whose row sum is exactly zero}
#'   \item{percent_isolated}{The same as a percentage}
#'   \item{saturation}{Mean row sum divided by N. Values near 1 mean the kernel
#'     is effectively constant and the bandwidth is too large to be informative.}
#'
#' @export
#'
#' @examples
#' data(ccRCCExample)
#' coords <- as.matrix(ccRCCExample[1:1500, c("x", "y")])
#'
#' # Coordinates are in millimetres. Sensible bandwidths, and one far too large:
#' qc_neighborhoods(coords, h = c(0.01, 0.05, 0.25, 25))
#' # The last row has saturation near 1 -- that bandwidth cannot detect anything.
qc_neighborhoods <- function(coords, h,
                             kernel = c("gaussian", "uniform"),
                             distance = c("euclidean", "manhattan"),
                             eps = 0) {
  kernel <- match.arg(kernel)
  distance <- match.arg(distance)

  if (!is.matrix(coords) || ncol(coords) < 2) {
    stop("coords must be a matrix with at least 2 columns")
  }
  if (!is.numeric(h) || length(h) < 1 || anyNA(h) || any(!is.finite(h)) || any(h <= 0)) {
    stop("h must be one or more finite positive bandwidths")
  }

  xy <- as.matrix(coords[, 1:2, drop = FALSE])
  N <- nrow(xy)

  # gene1 is irrelevant to row sums, but the kernel pass takes one; rowsums_only skips the second
  # pass entirely so nothing is computed that is then discarded.
  pass <- moran_kernel_pass(
    coords = xy, gene1 = rep(0, N), h_vec = as.numeric(h),
    eps = eps, kernel = kernel, dist_method = distance,
    n_threads = moran_max_threads(), rowsums_only = TRUE
  )

  rs <- pass$rowsums                      # n_h x N
  data.frame(
    h                = as.numeric(h),
    mean_neighbors   = rowMeans(rs),
    sd_neighbors     = apply(rs, 1, stats::sd),
    min_neighbors    = apply(rs, 1, min),
    max_neighbors    = apply(rs, 1, max),
    isolated_cells   = as.integer(pass$n_isolated),
    percent_isolated = round(100 * pass$n_isolated / N, 2),
    saturation       = round(rowMeans(rs) / N, 4),
    row.names        = NULL
  )
}
