// Multi-bandwidth kernel weighting for bivariate Moran's I.
//
// This replaces src/moran_multi_cpp.cpp. Three things drove the rewrite.
//
// 1. MEMORY. The previous version allocated dist_sq(N,N) and dist(N,N) unconditionally, before it
//    had even decided which mode to use, then added K(N,N) and an N x N temporary in dense mode:
//    four N x N doubles, about 1.6 GB at N = 7122 and 3.2 GB just under its N = 10000 threshold.
//    Its "large dataset" path was worse still, pushing roughly N^2 triplets at the default eps = 0.
//    Nothing here allocates an N x N buffer. Working memory is O(N * n_h), about 1.1 MB at
//    N = 7122 with six bandwidths, which is what makes N in the hundreds of thousands feasible.
//
// 2. A SILENT CORRECTNESS BUG. The old dense Gaussian branch read dist_sq unconditionally, but
//    dist_sq was only populated when distance == "euclidean". So gaussian + manhattan below the
//    N = 10000 threshold computed the kernel from uninitialised heap memory: three identical calls
//    returned 4.5e-05, then -5.3e-04, then NaN, against a correct value of 2.6e-02, and one of them
//    returned the same number for every bandwidth. Here there is no separate per-metric buffer to
//    fall out of step: each pair's squared distance is computed inline from the requested metric,
//    so the defect is not expressible.
//
// 3. DIVERGENT eps AND NaN HANDLING. eps was applied in the sparse path and ignored entirely in the
//    dense one, so the two disagreed on the row sums and therefore on I. Zero row sums were
//    unguarded in dense, guarded in sparse, and patched a third way in the R fallback. There is now
//    one path, so the predicate and the policy cannot diverge from themselves.
//
// STRUCTURE. The statistic is I(h) = <v_h, y> / sum(x^2) with v_h = K(h) (x / rowsum(h)). The row
// sums must be complete before the division, so this cannot be fused into a single sweep over
// pairs; it is genuinely two passes. Pass one accumulates row sums for every bandwidth at once,
// pass two accumulates v_h for every bandwidth at once. Both compute each pair distance once and
// reuse it across all bandwidths, which is where the saving over the reference comes from: the
// reference rebuilds a full kernel matrix per bandwidth.
//
// This function deliberately knows nothing about gene2 or permutations. It returns V, and the
// caller contracts it with the observed gene and with permuted genes. That keeps the permutation
// step in BLAS, lets it run in batches instead of materialising an N x perms matrix, and hands the
// analytical-moment code the same V without recomputing the kernel.
//
// LAYOUT. V and the row sums are stored n_h x N, not N x n_h. Column i then holds all bandwidths
// for cell i contiguously, which is the order the inner loop walks. It also means V %*% y comes out
// as n_h x 1 and V %*% Y_perm as n_h x perms, already the orientation the result object wants, with
// no transpose anywhere.

#include <RcppEigen.h>
#include <cmath>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

// Metric and kernel identifiers. Integers rather than strings so the O(N^2) inner loop never
// compares strings, and an enum rather than booleans so adding a third option cannot silently
// reinterpret an existing one.
enum Metric  { METRIC_EUCLIDEAN = 0, METRIC_MANHATTAN = 1, METRIC_MINKOWSKI = 2 };
enum KernelF { KERNEL_GAUSSIAN = 0, KERNEL_UNIFORM = 1, KERNEL_EXPONENTIAL = 2 };

// Squared distance in the requested metric.
//
// Squared, uniformly, so that a Gaussian kernel reads exp(-d^2 / 2h^2) and a uniform kernel reads
// d <= h under every metric with no per-metric special casing. This is the single definition; the
// previous code kept separate per-metric buffers, and the Gaussian branch read the one that had not
// been filled.
inline double squared_distance(double dx, double dy, int metric, double p, double inv_p) {
  switch (metric) {
    case METRIC_EUCLIDEAN: return dx * dx + dy * dy;
    case METRIC_MANHATTAN: {
      const double d1 = std::fabs(dx) + std::fabs(dy);
      return d1 * d1;
    }
    default: {
      // Minkowski of order p. p = 1 and p = 2 are handled by the cases above, so this path is only
      // taken for other orders and the pow() cost is confined to it.
      const double dp = std::pow(std::pow(std::fabs(dx), p) + std::pow(std::fabs(dy), p), inv_p);
      return dp * dp;
    }
  }
}

inline double kernel_weight(double dsq, int kernel, double inv_two_h_sq, double h_sq,
                            double inv_h) {
  switch (kernel) {
    case KERNEL_GAUSSIAN:    return std::exp(-dsq * inv_two_h_sq);
    case KERNEL_UNIFORM:     return (dsq <= h_sq) ? 1.0 : 0.0;
    // Exponential decay in the distance itself, so it needs the square root the other two avoid.
    //
    // Note the shape relative to the Gaussian at the same h: exp(-d/h) < exp(-d^2/2h^2) for d < 2h
    // and greater beyond it. So it decays FASTER near the origin and has a heavier tail, and the
    // crossover is at d = 2h. In practice the near-field dominates the row sum, so at a given h the
    // exponential gives a SMALLER effective neighbourhood, not a larger one.
    default:                 return std::exp(-std::sqrt(dsq) * inv_h);
  }
}

} // namespace

//' Kernel row sums and spatially lagged gene1, for every bandwidth
//'
//' Internal. Computes, for each bandwidth h, the row sums of the kernel weight matrix and the
//' vector v_h = K(h) (gene1 / rowsum(h)) which the caller contracts with gene2.
//'
//' @param coords N x 2 numeric matrix of coordinates.
//' @param gene1 Length-N numeric vector, already centred and scaled.
//' @param h_vec Numeric vector of bandwidths.
//' @param eps Weights are kept when \code{w > eps}, and excluded from the row sum otherwise.
//'   \code{eps = 0} keeps every non-underflowed weight, which is the documented default and the
//'   condition under which this reduces exactly to the reference implementation.
//' @param kernel "gaussian", "uniform" or "exponential".
//' @param dist_method "euclidean", "manhattan" or "minkowski".
//' @param p Order of the Minkowski distance; ignored by the other metrics.
//' @param drop_isolated When TRUE, cells with a zero row sum contribute nothing rather than
//'   producing a non-finite value. The caller uses this to implement nan_policy.
//' @param n_threads Requested OpenMP threads; ignored when built without OpenMP support.
//' @param rowsums_only When TRUE, skip the second pass and return the full row sum matrix instead
//'   of \code{V}. Used by \code{qc_neighborhoods()}, which needs the neighbourhood sizes but not
//'   the lagged gene.
//'
//' @return A list with \code{V} (n_h x N) or \code{rowsums} (n_h x N), plus
//'   \code{rowsum_mean}, \code{rowsum_min}, \code{n_isolated} and \code{h}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List moran_kernel_pass(const Eigen::Map<Eigen::MatrixXd> coords,
                             const Eigen::Map<Eigen::VectorXd> gene1,
                             const Eigen::Map<Eigen::VectorXd> h_vec,
                             double eps = 0.0,
                             std::string kernel = "gaussian",
                             std::string dist_method = "euclidean",
                             bool drop_isolated = true,
                             int n_threads = 1,
                             bool rowsums_only = false,
                             double p = 2.0) {

  const int N   = static_cast<int>(coords.rows());
  const int n_h = static_cast<int>(h_vec.size());

  if (coords.cols() < 2) Rcpp::stop("`coords` must have at least two columns.");
  if (gene1.size() != N) Rcpp::stop("`gene1` length must equal nrow(coords).");
  if (n_h < 1)           Rcpp::stop("`h_vec` must contain at least one bandwidth.");

  int kern;
  if      (kernel == "gaussian")    kern = KERNEL_GAUSSIAN;
  else if (kernel == "uniform")     kern = KERNEL_UNIFORM;
  else if (kernel == "exponential") kern = KERNEL_EXPONENTIAL;
  else Rcpp::stop("`kernel` must be \"gaussian\", \"uniform\" or \"exponential\".");

  int metric;
  if      (dist_method == "euclidean") metric = METRIC_EUCLIDEAN;
  else if (dist_method == "manhattan") metric = METRIC_MANHATTAN;
  else if (dist_method == "minkowski") {
    if (!(p >= 1.0) || !std::isfinite(p)) {
      Rcpp::stop("`p` must be a finite number >= 1 for the Minkowski distance.");
    }
    // p = 1 and p = 2 are the Manhattan and Euclidean cases; route them to the cheap paths rather
    // than paying for pow() to compute the same thing.
    if (p == 1.0)      metric = METRIC_MANHATTAN;
    else if (p == 2.0) metric = METRIC_EUCLIDEAN;
    else               metric = METRIC_MINKOWSKI;
  }
  else Rcpp::stop("`dist_method` must be \"euclidean\", \"manhattan\" or \"minkowski\".");

  const double inv_p = 1.0 / p;

  // Per-bandwidth constants, hoisted out of the O(N^2) loops.
  std::vector<double> inv_two_h_sq(n_h), h_sq(n_h), inv_h(n_h);
  for (int k = 0; k < n_h; ++k) {
    const double h = h_vec(k);
    if (!(h > 0) || !std::isfinite(h)) {
      Rcpp::stop("every bandwidth in `h_vec` must be finite and strictly positive.");
    }
    inv_two_h_sq[k] = 1.0 / (2.0 * h * h);
    h_sq[k]         = h * h;
    inv_h[k]        = 1.0 / h;
  }

#ifdef _OPENMP
  if (n_threads > 1) omp_set_num_threads(n_threads);
#endif

  const double* cx = coords.data();          // column 0
  const double* cy = coords.data() + N;      // column 1 (column-major storage)

  // ---------------------------------------------------------------------
  // Pass 1: row sums of the kernel, all bandwidths at once.
  //
  // Each iteration owns row i exclusively, so there is no reduction buffer and no race. The
  // diagonal is skipped rather than zeroed afterwards: a self-weight is not meaningful for
  // single cells, and the reference sets diag(K) <- 0 for the same reason.
  //
  // Summation runs over j in ascending order within each row, matching the order R's rowSums()
  // uses on the reference's kernel matrix, so the floating-point result agrees as closely as the
  // two distance formulations allow.
  // ---------------------------------------------------------------------
  Eigen::MatrixXd rowsums(n_h, N);
  rowsums.setZero();

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < N; ++i) {
    const double xi = cx[i], yi = cy[i];
    double* acc = rowsums.data() + static_cast<std::size_t>(i) * n_h;
    for (int j = 0; j < N; ++j) {
      if (j == i) continue;
      const double dsq = squared_distance(xi - cx[j], yi - cy[j], metric, p, inv_p);
      for (int k = 0; k < n_h; ++k) {
        const double w = kernel_weight(dsq, kern, inv_two_h_sq[k], h_sq[k], inv_h[k]);
        if (w > eps) acc[k] += w;
      }
    }
  }

  // Diagnostics, and the isolated-cell census that drives nan_policy upstream.
  Eigen::VectorXd rowsum_mean(n_h), rowsum_min(n_h);
  Eigen::VectorXi n_isolated(n_h);
  for (int k = 0; k < n_h; ++k) {
    double s = 0.0, mn = R_PosInf;
    int iso = 0;
    for (int i = 0; i < N; ++i) {
      const double r = rowsums(k, i);
      s += r;
      if (r < mn) mn = r;
      if (r <= 0.0) ++iso;
    }
    rowsum_mean(k) = s / N;
    rowsum_min(k)  = mn;
    n_isolated(k)  = iso;
  }

  if (rowsums_only) {
    return Rcpp::List::create(
      Rcpp::Named("rowsums")     = rowsums,
      Rcpp::Named("rowsum_mean") = rowsum_mean,
      Rcpp::Named("rowsum_min")  = rowsum_min,
      Rcpp::Named("n_isolated")  = n_isolated,
      Rcpp::Named("h")           = h_vec
    );
  }

  // ---------------------------------------------------------------------
  // gene1 divided by its row sum, per bandwidth.
  //
  // An isolated cell has row sum 0. Dividing anyway is what the reference does, and it lets a
  // single such cell turn the whole bandwidth into NaN, which the reference then maps to 0. That
  // behaviour is reproduced upstream in R where it can be explained and warned about, rather than
  // by propagating a non-finite value through the rest of this routine. Here the contribution is
  // simply zeroed so nothing downstream sees Inf or NaN.
  // ---------------------------------------------------------------------
  Eigen::MatrixXd u(n_h, N);
  for (int i = 0; i < N; ++i) {
    const double g = gene1(i);
    for (int k = 0; k < n_h; ++k) {
      const double r = rowsums(k, i);
      if (r > 0.0) {
        u(k, i) = g / r;
      } else {
        u(k, i) = drop_isolated ? 0.0 : R_NaN;
      }
    }
  }

  // ---------------------------------------------------------------------
  // Pass 2: v_h = K(h) %*% u(h), all bandwidths at once.
  //
  // Distances are recomputed rather than cached. Caching them would mean holding an N x N matrix,
  // which is precisely the allocation this rewrite exists to remove; the recomputation costs one
  // extra sweep of arithmetic and buys a reduction in peak memory of three orders of magnitude.
  // ---------------------------------------------------------------------
  Eigen::MatrixXd V(n_h, N);
  V.setZero();

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < N; ++i) {
    const double xi = cx[i], yi = cy[i];
    double* acc = V.data() + static_cast<std::size_t>(i) * n_h;
    for (int j = 0; j < N; ++j) {
      if (j == i) continue;
      const double dsq = squared_distance(xi - cx[j], yi - cy[j], metric, p, inv_p);
      const double* uj = u.data() + static_cast<std::size_t>(j) * n_h;
      for (int k = 0; k < n_h; ++k) {
        const double w = kernel_weight(dsq, kern, inv_two_h_sq[k], h_sq[k], inv_h[k]);
        if (w > eps) acc[k] += w * uj[k];
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("V")           = V,
    Rcpp::Named("rowsum_mean") = rowsum_mean,
    Rcpp::Named("rowsum_min")  = rowsum_min,
    Rcpp::Named("n_isolated")  = n_isolated,
    Rcpp::Named("h")           = h_vec
  );
}

//' Nearest-neighbour distance for every point
//'
//' Internal. One pass per point, so O(N) memory rather than the O(N^2) a full
//' \code{as.matrix(dist(coords))} would need. \code{suggest_bandwidths()} previously built that
//' full matrix, which made it unusable at the sample sizes this package targets: about 400 MB for
//' the 7122-cell bundled example and 80 GB at N = 100000.
//'
//' @param coords N x 2 numeric matrix.
//' @param dist_method "euclidean" or "manhattan".
//' @param n_threads Requested OpenMP threads; ignored without OpenMP.
//' @return Length-N numeric vector of distances to the closest other point.
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd nn_distance(const Eigen::Map<Eigen::MatrixXd> coords,
                            std::string dist_method = "euclidean",
                            int n_threads = 1) {
  const int N = static_cast<int>(coords.rows());
  if (coords.cols() < 2) Rcpp::stop("`coords` must have at least two columns.");
  if (N < 2) Rcpp::stop("at least two points are needed for a nearest-neighbour distance.");

  int metric;
  if      (dist_method == "euclidean") metric = METRIC_EUCLIDEAN;
  else if (dist_method == "manhattan") metric = METRIC_MANHATTAN;
  else Rcpp::stop("`dist_method` must be \"euclidean\" or \"manhattan\".");

#ifdef _OPENMP
  if (n_threads > 1) omp_set_num_threads(n_threads);
#endif

  const double* cx = coords.data();
  const double* cy = coords.data() + N;

  Eigen::VectorXd out(N);

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < N; ++i) {
    const double xi = cx[i], yi = cy[i];
    double best = R_PosInf;
    for (int j = 0; j < N; ++j) {
      if (j == i) continue;
      const double dsq = squared_distance(xi - cx[j], yi - cy[j], metric, 2.0, 0.5);
      if (dsq < best) best = dsq;
    }
    // squared_distance returns the square of the metric distance under both metrics, so one sqrt
    // recovers the distance itself.
    out(i) = std::sqrt(best);
  }

  return out;
}

//' Report whether the package was compiled with OpenMP support
//'
//' Internal. Used by tests and by the benchmark harness so a timing figure is never reported
//' without saying whether it was threaded.
//'
//' @return An integer: the maximum number of threads available, or 1 without OpenMP.
//' @keywords internal
// [[Rcpp::export]]
int moran_max_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}
