#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;
using namespace std;



// [[Rcpp::export]]
arma::sp_umat get_fixedNumber_neighbors(const arma::mat&  x, int  number=6)	{
  // Get the fixed number of neighbors for each spot.
  int i,j, jj, N = x.n_rows;
  int number4 = 4* number+1;
  arma::sp_umat D(N, N);
  //float dis;
  uvec idx, idx2, idx_all;
  for (j = 0; j < N; ++j)
  {   
    //Rprintf("good R!\n");
    idx_all = sort_index(abs(x(j,0) - x.col(0)));
    idx = idx_all.subvec(0, number4);
    //Rprintf("good R %d !\n",j);
    idx2 = find(idx != j);
    int p = idx2.n_elem;
    vec dis(p);
    
    for (i = 0; i < p; ++i)
    {
      dis(i) = norm(x.row(idx(idx2(i))) - x.row(j), 2);
    }
    uvec idx3 = sort_index(dis);
    for(i=0; i<number; ++i){
      jj = idx(idx2(idx3(i)));
      // D(jj,j) = 1;
      D(jj,j) = 1;
    }
  }
  return D;
}


