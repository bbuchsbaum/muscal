// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' @title Fast ridge regression solver
//' @description Solves the ridge regression problem using RcppArmadillo.
//' @param Z The predictor matrix.
//' @param X The response matrix.
//' @param lambda The regularization parameter.
//' @return The coefficient matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::mat ridge_solve(const arma::mat& Z, const arma::mat& X, double lambda) {
    arma::mat ZtZ = Z.t() * Z;
    ZtZ.diag() += lambda;
    arma::mat ZtX = Z.t() * X;
    return arma::solve(ZtZ, ZtX, arma::solve_opts::fast);
} 