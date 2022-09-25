// RcppArmadillo utility function to calculate observed Fisher 
// information matrix of multinomial fit, with 
// probs=fitted probabilities (with 1st category/column dropped)
// Z = model matrix
// row_totals = row totals
// We do this using Kronecker products, as in
// https://ieeexplore.ieee.org/abstract/document/1424458
// B. Krishnapuram; L. Carin; M.A.T. Figueiredo; A.J. Hartemink
// Sparse multinomial logistic regression: fast algorithms and
// generalization bounds
// IEEE Transactions on Pattern Analysis and Machine
// Intelligence ( Volume: 27, Issue: 6, June 2005)

#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat calc_infmatrix_RcppArma(arma::mat probs, arma::mat Z, arma::vec row_totals) {
  int n = Z.n_rows;
  int p = Z.n_cols;
  int k = probs.n_cols;
  int ncoefs = k * p;
  arma::mat info = arma::zeros<arma::mat>(ncoefs, ncoefs);
  arma::mat diag_probs;
  arma::mat tcrossprod_probs;
  arma::mat tcrossprod_Z;
  arma::mat kronecker_prod;
  for (int i = 0; i < n; i++) {
    diag_probs = arma::diagmat(probs.row(i));
    tcrossprod_probs = arma::trans(probs.row(i)) * probs.row(i);
    tcrossprod_Z = (arma::trans(Z.row(i)) * Z.row(i)) * row_totals(i);
    kronecker_prod = arma::kron(diag_probs - tcrossprod_probs, tcrossprod_Z);
    info += kronecker_prod;
  }
  return info;
}