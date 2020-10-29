#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// ROUTINES
// (01) routine_mds : multidimensional scaling


// (01) routine_mds : multidimensional scaling ---------------------------------
// [[Rcpp::export]]
arma::mat routine_mds(arma::mat& dmat, int ndim){
  int N = dmat.n_rows;
  arma::mat D2 = arma::pow(dmat, 2.0);
  arma::mat J  = arma::eye<arma::mat>(N,N) - (arma::ones<arma::mat>(N,N)/(static_cast<double>(N)));
  arma::mat B  = -0.5*J*D2*J;
  
  arma::vec eigval;
  arma::mat eigvec;
  
  arma::eig_sym(eigval, eigvec, B);
  arma::mat Y = eigvec.tail_cols(ndim)*arma::diagmat(arma::sqrt(eigval.tail(ndim)));
  return(Y);
}