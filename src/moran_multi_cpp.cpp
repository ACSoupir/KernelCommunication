#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace Eigen;

// Use RowMajor for efficient row operations
using SparseMatrixRM = Eigen::SparseMatrix<double, Eigen::RowMajor>;

// [[Rcpp::export]]
Rcpp::List moran_multi_cpp(
    const Eigen::MatrixXd& coords,      // N x 2
    const Eigen::VectorXd& gene1,       // N
    const Eigen::VectorXd& gene2,       // N
    const Eigen::VectorXd& h_vec,       // n_h bandwidths
    const Eigen::MatrixXd& Y_perm,      // N x perms (pre-computed in R)
    int perms,                          // total permutations
    int seed,                           // random seed (not used, kept for API)
    double eps,                         // epsilon threshold
    int batch_size,                     // permutations per batch
    bool verbose,                       // progress reporting
    int large_mode_n,                   // N threshold for large mode
    std::string kernel,                 // "gaussian" or "uniform"
    std::string dist_method             // "euclidean" or "manhattan"
) {

  const int N = coords.rows();
  const int n_h = h_vec.size();

  // Pre-compute invariant quantities
  const double denom = gene1.squaredNorm();

  // Determine if we should use dense mode (faster for small N)
  // Dense mode is faster for N < 10000 due to sparse matrix overhead
  bool use_dense = (N < 10000);

  // Compute distance matrix (once)
  Eigen::MatrixXd dist_sq(N, N);
  Eigen::MatrixXd dist(N, N);

  if (dist_method == "euclidean") {
    // Squared Euclidean distances
    for (int i = 0; i < N; i++) {
      for (int j = i; j < N; j++) {
        double dx = coords(i, 0) - coords(j, 0);
        double dy = coords(i, 1) - coords(j, 1);
        double d2 = dx * dx + dy * dy;
        dist_sq(i, j) = d2;
        dist_sq(j, i) = d2;
      }
    }
  } else if (dist_method == "manhattan") {
    // Manhattan distances
    for (int i = 0; i < N; i++) {
      for (int j = i; j < N; j++) {
        double dx = std::abs(coords(i, 0) - coords(j, 0));
        double dy = std::abs(coords(i, 1) - coords(j, 1));
        double d = dx + dy;
        dist(i, j) = d;
        dist(j, i) = d;
      }
    }
  }

  // Output vectors
  Eigen::VectorXd observed(n_h);
  Eigen::MatrixXd all_perm_vals(n_h, perms);

  // Progress reporting
  int progress_interval = (n_h > 1) ? std::max(1, n_h / 20) : 0;

  // Loop over bandwidths
  for (int h_idx = 0; h_idx < n_h; h_idx++) {
    double h = h_vec(h_idx);

    // Declare variables for both modes
    Eigen::MatrixXd W_norm;  // Used in dense mode
    SparseMatrixRM W;        // Used in sparse mode
    Eigen::VectorXd rowsums(N);
    Eigen::VectorXd colsums(N);

    if (use_dense) {
      // DENSE MODE: Match soupir_moran_multi algorithm
      Eigen::MatrixXd K(N, N);
      rowsums.setZero();
      colsums.setZero();

      if (kernel == "gaussian") {
        const double inv2h2 = 1.0 / (2.0 * h * h);

        // Vectorized kernel computation using Eigen array operations
        Eigen::MatrixXd exp_arg = -dist_sq * inv2h2;
        K = exp_arg.array().exp();
        K.diagonal().setZero();
        
        // Compute rowsums
        rowsums = K.rowwise().sum();
        
      } else if (kernel == "uniform") {
        const double h_sq = h * h;

        // Vectorized kernel computation
        if (dist_method == "euclidean") {
          K = (dist_sq.array() <= h_sq).cast<double>();
        } else {
          K = (dist.array() <= h).cast<double>();
        }
        K.diagonal().setZero();
        
        // Compute rowsums
        rowsums = K.rowwise().sum();
      }

      // Compute weighted gene1
      Eigen::VectorXd gene1_weighted = gene1.array() / rowsums.array();
      
      // v_h = K * gene1_weighted
      Eigen::VectorXd v_h = K * gene1_weighted;

      // Compute observed Moran's I: I = v_h^T * gene2 / denom
      observed(h_idx) = (v_h.dot(gene2)) / denom;

      // Compute permutations: I_perm = v_h^T * Y_perm / denom
      if (perms > 0) {
        Eigen::RowVectorXd perm_vals = (v_h.transpose() * Y_perm) / denom;
        all_perm_vals.row(h_idx) = perm_vals;
      }

    } else {
      // SPARSE MODE: Original sparse implementation
      // Build sparse kernel weight matrix
      std::vector<Triplet<double>> triplets;
      triplets.reserve(static_cast<size_t>(N) * 16);

      rowsums.setZero();
      colsums.setZero();

      int isolated_count = 0;

      if (kernel == "gaussian") {
        const double inv2h2 = 1.0 / (2.0 * h * h);

        // Exploit symmetry: compute upper triangle and mirror
        for (int i = 0; i < N; i++) {
          for (int j = i + 1; j < N; j++) {
            double d2 = (dist_method == "euclidean") ? dist_sq(i, j) : dist(i, j) * dist(i, j);
            double w = std::exp(-d2 * inv2h2);

            if (w > eps) {
              rowsums(i) += w;
              rowsums(j) += w;  // Symmetric
              triplets.push_back(Triplet<double>(i, j, w));
              triplets.push_back(Triplet<double>(j, i, w));
            }
          }
        }
      } else if (kernel == "uniform") {
        const double h_sq = h * h;

        // Exploit symmetry: compute upper triangle and mirror
        for (int i = 0; i < N; i++) {
          for (int j = i + 1; j < N; j++) {
            double d2 = (dist_method == "euclidean") ? dist_sq(i, j) : dist(i, j) * dist(i, j);

            if (d2 <= h_sq) {
              rowsums(i) += 1.0;
              rowsums(j) += 1.0;  // Symmetric
              triplets.push_back(Triplet<double>(i, j, 1.0));
              triplets.push_back(Triplet<double>(j, i, 1.0));
            }
          }
        }
      }

      // Build sparse matrix (row-major for efficient row operations)
      SparseMatrixRM W(N, N);
      W.setFromTriplets(triplets.begin(), triplets.end());

      // Row-normalize (RowMajor: InnerIterator iterates over columns in each row)
      for (int i = 0; i < N; i++) {
        if (rowsums(i) > 0) {
          for (SparseMatrixRM::InnerIterator it(W, i); it; ++it) {
            it.valueRef() /= rowsums(i);
          }
        } else {
          isolated_count++;
        }
      }

      // Warn about isolated points
      if (isolated_count > 0 && verbose) {
        Rcpp::Rcout << "Warning: " << isolated_count
                    << " cells have no neighbors within bandwidth\n";
      }

      // Compute observed Moran's I
      // U = W * gene2
      Eigen::VectorXd U = W * gene2;
      observed(h_idx) = (gene1.dot(U)) / denom;

      // Compute permutations (batched) - use pre-computed Y_perm
      if (perms > 0) {
        // Batched multiplication
        int n_batches = (perms + batch_size - 1) / batch_size;

        for (int b = 0; b < n_batches; b++) {
          int batch_start = b * batch_size;
          int batch_end = std::min(batch_start + batch_size, perms);
          int current_batch_size = batch_end - batch_start;

          Eigen::MatrixXd Y_batch = Y_perm.block(0, batch_start, N, current_batch_size);
          Eigen::MatrixXd U_batch = W * Y_batch;

          for (int p = 0; p < current_batch_size; p++) {
            all_perm_vals(h_idx, batch_start + p) = gene1.dot(U_batch.col(p)) / denom;
          }
        }
      }
    }  // end of sparse mode

    // Progress reporting
    if (verbose && n_h > 1) {
      if ((h_idx + 1) % progress_interval == 0 || h_idx == n_h - 1) {
        int pct = (h_idx + 1) * 100 / n_h;
        Rcpp::Rcout << "\r  Computing bandwidths: " << pct
                    << "% (" << (h_idx + 1) << "/" << n_h << ")"
                    << std::flush;
        if (h_idx == n_h - 1) {
          Rcpp::Rcout << "\n" << std::flush;
        }
      }
    }
  }

  // Build output list
  Rcpp::List output;
  output["h"] = h_vec;
  output["observed"] = observed;
  output["permutations"] = all_perm_vals;

  return output;
}
