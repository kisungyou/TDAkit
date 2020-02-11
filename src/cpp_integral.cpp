#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double simple_integral(arma::mat &lmat, arma::vec tseq){
  // parameters and setup
  int TT = lmat.n_rows;
  int KK = lmat.n_cols;
  
  arma::vec outval(KK,fill::zeros); // recording norms
  arma::vec tmpvec(TT,fill::zeros);
  for (int k=0;k<KK;k++){ // iterate over columns
    tmpvec = lmat.col(k); 
    double record = 0.0;
    for (int t=0;t<(TT-1);t++){
      record += (tseq(t+1)-tseq(t))*(tmpvec(t+1)+tmpvec(t))/2.0;
    }
    outval(k) = record;
  }
  return(arma::accu(outval));
}
// [[Rcpp::export]]
double simple_integral_1d(arma::vec lvec, arma::vec tseq){
  // parameters and setup
  int TT = tseq.n_elem;
  double output = 0.0;
  for (int t=0;t<(TT-1);t++){
    output += (tseq(t+1)-tseq(t))*(lvec(t+1)+lvec(t))/2.0;
  }
  return(output);
}