#include <RcppEigen.h>
#include <algorithm>
#include <random>
#include <limits>

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Rcpp::List soupir_moran_bv_multi(const Eigen::MatrixXd& coords,
                          const Rcpp::NumericVector& gene1,
                          const Rcpp::NumericVector& gene2,
                          const Rcpp::NumericVector& h_vec,
                          const int perms,
                          const std::string transformation = "gaussian",
                          const double eps = 0)
{
  const int n = coords.rows();
  const int H = h_vec.size();

  if (gene1.size() != n || gene2.size() != n) {
    stop("gene1 and gene2 must have the same length as the number of rows in coords");
  }

  // scale gene1
  double mean1 = Rcpp::mean(gene1);
  double sd1 = Rcpp::sd(gene1);
  VectorXd g1s(n);
  for (int i = 0; i < n; ++i)
    g1s[i] = (gene1[i] - mean1) / sd1;

  // scale gene2
  double mean2 = Rcpp::mean(gene2);
  double sd2 = Rcpp::sd(gene2);
  VectorXd g2s(n);
  for (int i = 0; i < n; ++i)
    g2s[i] = (gene2[i] - mean2) / sd2;

  double denom = g1s.squaredNorm();

  Rcpp::NumericVector observed(H);
  Rcpp::NumericMatrix permutations(perms, H);

  std::mt19937 rng(static_cast<unsigned int>(R::unif_rand() * std::numeric_limits<unsigned int>::max()));
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);

  for (int h_idx = 0; h_idx < H; ++h_idx) {
    double h = h_vec[h_idx];
    double inv2h2 = 1.0 / (2.0 * h * h);
    double invh = 1.0 / h;

    std::vector<Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(n) * 16);

    for (int i = 0; i < n; ++i) {
      double rowsum = 0.0;
      std::vector<std::pair<int, double>> row_weights;

      for (int j = 0; j < n; ++j) {
        if (i == j) continue;
        double dx = coords(i, 0) - coords(j, 0);
        double dy = coords(i, 1) - coords(j, 1);
        double d2 = dx * dx + dy * dy;
        double d = std::sqrt(d2);

        double w;
        if (transformation == "gaussian") {
          w = std::exp(-d2 * inv2h2);
        } else if (transformation == "exponential") {
          w = std::exp(-d * invh);
        } else if (transformation == "inverse") {
          w = 1.0 / (1.0 + d);
        } else if (transformation == "uniform") {
          w = 1.0;
        } else {
          stop("Unsupported transformation type: " + transformation);
        }

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

    VectorXd Wy = W * g2s;
    double obs = g1s.cwiseProduct(Wy).sum() / denom;
    observed[h_idx] = obs;

    for (int p = 0; p < perms; ++p) {
      std::shuffle(idx.begin(), idx.end(), rng);
      VectorXd g2_perm(n);
      for (int i = 0; i < n; ++i)
        g2_perm[i] = g2s[idx[i]];

      VectorXd Wy_perm = W * g2_perm;
      permutations(p, h_idx) = g1s.cwiseProduct(Wy_perm).sum() / denom;
    }
  }

  return Rcpp::List::create(
    _["observed"] = observed,
    _["permutations"] = permutations
  );
}
