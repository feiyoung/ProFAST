// Revised log:
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
//#include<boost/math/tools/minima.hpp>

#define INT_MIN1 (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;
// using boost::math::tools::brent_find_minima;

// [[Rcpp::export]]
arma::mat pdistance_cpp(const arma::mat& Ar, const arma::mat& Br, const float& eta=1e-10){
  
  vec An = sum(Ar % Ar, 1);
  vec Bn = sum(Br % Br, 1);
  mat C = -2* Ar * Br.t();
  C.each_col() += An;
  C.each_row() += Bn.t();
  return sqrt(C + eta);
}

// [[Rcpp::export]]
arma::mat gene_embed_cpp(const arma::mat& X,  const arma::mat& ce_cell) {
  // X is a gene by cell expression matrix
  mat weight = X;
  colvec sumW = sum(weight, 1);
  weight.each_col() /= sumW;
  
  return weight * ce_cell;
}

