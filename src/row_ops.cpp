#include <RcppEigen.h>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

static Eigen::MatrixXd apply_row_system_checked(const List& A_list,
                                               const Eigen::Map<Eigen::MatrixXd>& S) {
  const int n = S.rows(), k = S.cols();
  if (A_list.size() != n) stop("A_list length must match the rows of S.");
  if (!S.allFinite()) stop("S must contain finite values.");
  Eigen::MatrixXd out(n, k);
  for (int i = 0; i < n; ++i) {
    SEXP item = A_list[i];
    if (!Rf_isMatrix(item) || (TYPEOF(item) != REALSXP && TYPEOF(item) != INTSXP))
      stop("Each A_list entry must be a finite numeric k-by-k matrix.");
    NumericMatrix A = as<NumericMatrix>(item);
    if (A.nrow() != k || A.ncol() != k)
      stop("Each A_list entry must have dimensions ncol(S) by ncol(S).");
    Eigen::Map<Eigen::MatrixXd> Ai(A.begin(), k, k);
    if (!Ai.allFinite()) stop("A_list matrices must contain finite values.");
    out.row(i) = (Ai * S.row(i).transpose()).transpose();
  }
  return out;
}

// [[Rcpp::export]]
SEXP muscal_apply_row_system_cpp(List A_list, Eigen::Map<Eigen::MatrixXd> S) {
  return wrap(apply_row_system_checked(A_list, S));
}

// [[Rcpp::export]]
SEXP muscal_apply_row_operator_cpp(List A_list,
                                  Eigen::Map<Eigen::MatrixXd> S,
                                  SEXP graph_laplacian,
                                  const double graph_lambda,
                                  const double ridge) {
  if (!std::isfinite(graph_lambda) || graph_lambda < 0 ||
      !std::isfinite(ridge) || ridge < 0) stop("graph_lambda and ridge must be finite and nonnegative.");
  Eigen::MatrixXd out = apply_row_system_checked(A_list, S);
  if (graph_laplacian != R_NilValue) {
    if (Rf_isS4(graph_laplacian)) {
      if (!Rf_inherits(graph_laplacian, "dgCMatrix")) stop("Sparse graph must be a dgCMatrix; coerce before calling the kernel.");
      Eigen::SparseMatrix<double> L = as<Eigen::SparseMatrix<double> >(graph_laplacian);
      if (L.rows() != S.rows() || L.cols() != S.rows()) stop("Graph dimensions must match nrow(S).");
      for (int col = 0; col < L.outerSize(); ++col)
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, col); it; ++it)
          if (!std::isfinite(it.value())) stop("Graph must contain finite values.");
      if (graph_lambda > 0) out.noalias() += graph_lambda * (L * S);
    } else {
      if (!Rf_isMatrix(graph_laplacian) ||
          (TYPEOF(graph_laplacian) != REALSXP && TYPEOF(graph_laplacian) != INTSXP))
        stop("Graph must be a finite numeric matrix.");
      NumericMatrix graph = as<NumericMatrix>(graph_laplacian);
      if (graph.nrow() != S.rows() || graph.ncol() != S.rows()) stop("Graph dimensions must match nrow(S).");
      Eigen::Map<Eigen::MatrixXd> L(graph.begin(), graph.nrow(), graph.ncol());
      if (!L.allFinite()) stop("Graph must contain finite values.");
      if (graph_lambda > 0) out.noalias() += graph_lambda * (L * S);
    }
  }
  if (ridge > 0) out.noalias() += ridge * S;
  return wrap(out);
}

// [[Rcpp::export]]
List muscal_rowsum_counts_cpp(Eigen::Map<Eigen::MatrixXd> X, SEXP idx, SEXP N) {
  if ((TYPEOF(N) != REALSXP && TYPEOF(N) != INTSXP) || Rf_xlength(N) != 1)
    stop("N must be a finite nonnegative integer scalar.");
  const double groups = as<double>(N);
  if (!std::isfinite(groups) || groups < 0 || groups != std::floor(groups) ||
      groups > std::numeric_limits<int>::max()) stop("N must be a finite nonnegative integer scalar.");
  if (!X.allFinite()) stop("X must contain finite values.");
  if ((TYPEOF(idx) != REALSXP && TYPEOF(idx) != INTSXP) || Rf_xlength(idx) != X.rows())
    stop("idx must be numeric with length nrow(X).");
  NumericVector indices = as<NumericVector>(idx);
  for (R_xlen_t i = 0; i < indices.size(); ++i) {
    const double row = indices[i];
    if (!std::isfinite(row) || row != std::floor(row) || row < 1 || row > groups)
      stop("idx values must be finite integers in 1..N.");
  }
  const int n = X.rows(), k = X.cols(), count = static_cast<int>(groups);
  Eigen::MatrixXd sums = Eigen::MatrixXd::Zero(count, k);
  IntegerVector counts(count);
  for (int i = 0; i < n; ++i) {
    const int row = static_cast<int>(indices[i]) - 1;
    sums.row(row) += X.row(i);
    counts[row] += 1;
  }
  return List::create(_["sums"] = sums, _["counts"] = counts);
}
