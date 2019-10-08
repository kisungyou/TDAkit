#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat compute_slicemean(arma::cube &dcube){
  // parameters
  int n = dcube.n_rows;
  int p = dcube.n_cols;
  
  // compute along the slice dimension
  arma::mat jiwoongboy = arma::mean(dcube, 2);
  return(jiwoongboy);
}

// [[Rcpp::export]]
arma::mat compute_slicewsum(arma::cube &dcube, arma::vec weight){
  // parameters
  int n = dcube.n_rows;
  int p = dcube.n_cols;
  int k = weight.n_elem;
  
  // summation of weight vector
  double allwsum = arma::accu(weight);
  
  // compute along the slice dimension
  arma::mat jiwoongboy(n,p,fill::zeros);
  arma::mat tgtslice(n,p,fill::zeros);
  for (int i=0;i<k;i++){
    tgtslice = dcube.slice(i);
    jiwoongboy = jiwoongboy + (weight(i)/allwsum)*tgtslice;
  }
  return(jiwoongboy);
}
