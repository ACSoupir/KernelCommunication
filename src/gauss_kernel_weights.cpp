#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Eigen::SparseMatrix<double> gaussian_kernel_weights(const Eigen::MatrixXd& coords,
                                                    const double h,
                                                    const double eps = 1e-15)
{
  const int n = coords.rows();
  const double inv2h2 = 1.0 / (2.0 * h * h);

  std::vector<Triplet<double>> triplets;
  triplets.reserve(static_cast<size_t>(n) * 16);  // heuristic

  for (int i = 0; i < n; ++i) {
    double rowsum = 0.0;
    std::vector<std::pair<int, double>> row_weights;

    for (int j = 0; j < n; ++j) {
      if (i == j) continue;

      double dx = coords(i, 0) - coords(j, 0);
      double dy = coords(i, 1) - coords(j, 1);
      double d2 = dx * dx + dy * dy;

      double w = std::exp(-d2 * inv2h2);  // Gaussian kernel over squared distance
      if (w > eps) {
        row_weights.emplace_back(j, w);
        rowsum += w;
      }
    }

    for (const auto& [j, w] : row_weights) {
      triplets.emplace_back(i, j, w / rowsum);
    }
  }

  SparseMatrix<double> W(n, n);
  W.setFromTriplets(triplets.begin(), triplets.end());
  return W;
}
