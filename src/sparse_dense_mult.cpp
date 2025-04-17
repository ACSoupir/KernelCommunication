#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd sparse_dense_mult(const Eigen::SparseMatrix<double> &Ws,
                                  const Eigen::MatrixXd &P) {
  // Multiply sparse matrix Ws with dense matrix P
  return Ws * P;
}
