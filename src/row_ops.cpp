#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
SEXP muscal_apply_row_system_cpp(List A_list, Eigen::Map<Eigen::MatrixXd> S) {
  const int n = S.rows();
  const int k = S.cols();
  Eigen::MatrixXd out(n, k);

  for (int i = 0; i < n; ++i) {
    Eigen::Map<Eigen::MatrixXd> Ai = as<Eigen::Map<Eigen::MatrixXd> >(A_list[i]);
    out.row(i) = (Ai * S.row(i).transpose()).transpose();
  }

  return wrap(out);
}

// [[Rcpp::export]]
SEXP muscal_apply_row_operator_cpp(List A_list,
                                   Eigen::Map<Eigen::MatrixXd> S,
                                   SEXP graph_laplacian,
                                   const double graph_lambda,
                                   const double ridge) {
  const int n = S.rows();
  const int k = S.cols();
  Eigen::MatrixXd out(n, k);

  for (int i = 0; i < n; ++i) {
    Eigen::Map<Eigen::MatrixXd> Ai = as<Eigen::Map<Eigen::MatrixXd> >(A_list[i]);
    out.row(i) = (Ai * S.row(i).transpose()).transpose();
  }

  if (graph_laplacian != R_NilValue && graph_lambda > 0.0) {
    if (Rf_isS4(graph_laplacian)) {
      Eigen::MappedSparseMatrix<double> L = as<Eigen::MappedSparseMatrix<double> >(graph_laplacian);
      out.noalias() += graph_lambda * (L * S);
    } else {
      Eigen::Map<Eigen::MatrixXd> L = as<Eigen::Map<Eigen::MatrixXd> >(graph_laplacian);
      out.noalias() += graph_lambda * (L * S);
    }
  }

  if (ridge > 0.0) {
    out.noalias() += ridge * S;
  }

  return wrap(out);
}

// [[Rcpp::export]]
List muscal_rowsum_counts_cpp(Eigen::Map<Eigen::MatrixXd> X,
                              IntegerVector idx,
                              const int N) {
  const int n = X.rows();
  const int k = X.cols();
  if (idx.size() != n) {
    stop("idx length must match the number of rows in X.");
  }

  Eigen::MatrixXd sums = Eigen::MatrixXd::Zero(N, k);
  IntegerVector counts(N);

  for (int i = 0; i < n; ++i) {
    const int row = idx[i] - 1;
    if (row < 0 || row >= N) {
      stop("idx values must be in 1..N.");
    }
    sums.row(row) += X.row(i);
    counts[row] += 1;
  }

  return List::create(
    _["sums"] = sums,
    _["counts"] = counts
  );
}
