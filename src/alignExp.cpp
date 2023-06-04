// Date: 2023-02-09

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>


#define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;


// correct for one gene

vec update_alphaj_alignExp(const field<vec>& Xf, const field<mat>& Hf, const field<mat>& Rf, 
                           const mat& Tm, const field<vec>& Eu, const vec& gammaj, const vec& zetaj,
                           const vec& sigmaj){
  
  int r, M=Xf.n_elem, d = Rf(0).n_cols;
  vec y_right(d);
  mat y_left(d,d, fill::zeros);
  for(r=0; r<M; ++r){
    y_left += Rf(r).t() * Rf(r)/sigmaj(r);
    y_right += Rf(r).t() * (Xf(r)- Hf(r)*gammaj  - Eu(r)- as_scalar(Tm.row(r) * zetaj)) /sigmaj(r);
  }
  
  return inv_sympd(y_left) * y_right;
}

//
vec update_gammaj_alignExp(const field<vec>& Xf, const field<mat>& Hf, const field<mat>& Rf, 
                           const mat& Tm, const field<vec>& Eu, const vec& alphaj, const vec& zetaj,
                           const vec& sigmaj){
  
  int r, M=Xf.n_elem, d = Hf(0).n_cols;
  vec y_right(d);
  mat y_left(d,d, fill::zeros);
  for(r=0; r<M; ++r){
    y_left += Hf(r).t() * Hf(r)/sigmaj(r);
    y_right += Hf(r).t() * (Xf(r)- Rf(r)*alphaj  - Eu(r)- as_scalar(Tm.row(r) * zetaj))/sigmaj(r);
  }
  
  return inv_sympd(y_left) * y_right;
}

vec update_zetaj_alignExp(const field<vec>& Xf, const field<mat>& Hf, const field<mat>& Rf, 
                          const mat& Tm, const field<vec>& Eu, const vec& alphaj, const vec& gammaj,
                          const vec& sigmaj){
  
  int r, M=Xf.n_elem, d = Tm.n_cols, n_r;
  vec y_right(d);
  mat y_left(d,d, fill::zeros);
  for(r=0; r<M; ++r){
    n_r = Xf(r).n_rows;
    y_left += Tm.row(r).t() * Tm.row(r)* n_r /sigmaj(r);
    y_right += Tm.row(r).t() * accu(Xf(r)-Hf(r)*gammaj - Rf(r)*alphaj  - Eu(r))/sigmaj(r);
  }
  
  return inv_sympd(y_left) * y_right;
}

vec update_sigmaj_alignExp(const field<vec>& Xf, const field<mat>& Hf, const field<mat>& Rf,
                           const mat& Tm, const field<vec>& Eu, const field<vec>& Varu,
                           const vec& alphaj, const vec& gammaj,const vec& zetaj){
      int r, M=Xf.n_elem, n_r;
      vec sigmaj(M);
      
      for(r=0; r<M; ++r){
        n_r = Xf(r).n_rows;
        vec  tmpvec = Xf(r)-Hf(r)*gammaj - Rf(r)*alphaj  - Eu(r)- as_scalar(Tm.row(r) * zetaj);
        sigmaj(r) = accu(tmpvec % tmpvec + Varu(r))  / n_r;
      }
      
  return sigmaj;
}

vec update_psij_alignExp(const field<vec>& Muu, const field<vec>& Eu, const field<vec>& Varu){
  
  int r, M=Muu.n_elem, n_r;
  vec psij(M);
  for(r=0; r<M; ++r){
    n_r = Muu(r).n_rows;
    psij(r) = accu((Muu(r)-Eu(r)) % (Muu(r)-Eu(r)) + Varu(r)) / n_r;
  }
  
  return(psij);
}

//
arma::vec get_Vmean_alignExp(const arma::vec& V, const arma::sp_mat& Adj){
  int i, n = V.n_rows;
  vec Uv(n, fill::zeros);
  for (i = 0; i < n; i++)
  {
    arma::vec col(Adj.col(i)); // the class label of neighbors of i-th sample.
    uvec q1 = find(col > 0);
    // cout<<q1.n_rows<<endl;
    if( q1.n_rows>0){
      Uv(i) = mean(V(q1));
    }
    
  }
  return Uv;
}

void E_step_alignExp(const field<vec>& Xf, const field<mat>& Hf, const field<mat>& Rf, 
                     const mat& Tm, const field<sp_mat>& Adjf, field<vec>& Muu,  field<vec>& Eu,  field<vec>& Varu,
                     const vec& alphaj, const vec& gammaj,const vec& zetaj,
                     const vec& sigmaj, const vec& psij){
  
  int r, M=Muu.n_elem, n_r;
  field<vec> Mu_x(Xf);
  for(r=0; r<M; ++r){
    n_r = Xf(r).n_rows;
    Mu_x(r) = Muu(r) + Rf(r) * alphaj + Hf(r) * gammaj + as_scalar(Tm.row(r) * zetaj);
    Eu(r) = Muu(r) + (Xf(r) - Mu_x(r))*psij(r)/(psij(r) + sigmaj(r));
    Varu(r) = as_scalar(psij(r) - psij(r)*psij(r)/(psij(r) + sigmaj(r))) * ones(n_r, 1);
  }
  // ICM-step:
  // ICM-step
  for(r = 0; r<M; ++r){
    Muu(r) =  get_Vmean_alignExp(Eu(r),  Adjf(r)); // update Muv
  }
  
}

// Q-function
double Q_fun(const field<vec>& Xf, const field<mat>& Hf, const field<mat>& Rf, 
             const mat& Tm, const field<vec>& Muu,  const field<vec>& Eu,  const field<vec>& Varu,
             const vec& alphaj, const vec& gammaj,const vec& zetaj,
             const vec& sigmaj, const vec& psij){
  
  int r, M=Xf.n_elem, n_r;
  double logPx = 0.0;
  for(r=0; r<M; ++r){
    n_r = Xf(r).n_rows;
    vec  tmpvec = Xf(r)-Hf(r)*gammaj - Rf(r)*alphaj - Eu(r)- as_scalar(Tm.row(r) * zetaj);
    logPx += -0.5* (n_r * log(sigmaj(r)) + accu(tmpvec % tmpvec + Varu(r))  / sigmaj(r));
  }
  
  double logPu = 0.0;
  double entropy = 0.0;
  for(r=0; r<M; ++r){
    n_r = Xf(r).n_rows;
    logPu += -0.5* (n_r * log(psij(r)) + accu((Muu(r)-Eu(r)) % (Muu(r)-Eu(r)) + Varu(r))  / psij(r));
    entropy += -0.5 * accu(log(Varu(r)));
  }
  
  
  return logPx + logPu - entropy;
}



// [[Rcpp::export]]
Rcpp::List correct_one_gene(const Rcpp::List& Xlist, const Rcpp::List& RList,
                            const Rcpp::List& HList, const arma::mat& Tm, const Rcpp::List& Adjlist, 
                            const arma::vec& sigmaj_int, const arma::vec& psij_int, const arma::vec& alphaj_int, 
                            const arma::vec& gammaj_int, const arma::vec& zetaj_int, const int& maxIter,
                            const double& epsELBO, const bool& verbose){
  // Xlist's component is a vector, denoting one gene's expression.
  int r,M = Xlist.length();
  if(M != Tm.n_rows){
    stop("The length of Xlist is not matched to the nrow of Tm!");
  }
  bool Tm_flag = true;
  if(accu(abs(Tm))<1e-10){
    Tm_flag=false;
  }else{
    Tm_flag= true;
  }
  // Initialize the model parameters
  vec sigmaj(sigmaj_int), psij(psij_int), alphaj(alphaj_int), gammaj(gammaj_int), zetaj(zetaj_int);
  
  field<sp_mat> Adjf(M);
  field<mat> Hf(M), Rf(M);
  field<vec> Eu(M), Varu(M), Xf(M); // require to initialize values not only shape!!!
  int n;
  for(r=0; r < M; ++r){ 
    vec Xtmp = Xlist[r]; // enforce to become a matrix.
    Xf(r) = Xtmp;
    mat Htmp = HList[r];
    Hf(r) = Htmp;
    mat Rtmp = RList[r];
    Rf(r) = Rtmp;
    
    sp_mat Adjtmp = Adjlist[r];
    Adjf(r) = Adjtmp;
    //mat tmpMat = Wm_List[r];
    //Wm0.slice(r) = tmpMat;
    
    Eu(r) = zeros(Xf(r).n_rows, 1);
    Varu(r) = ones(Xf(r).n_rows, 1);
  }
  field<vec> Muu(Eu);
  vec elbo_vec(maxIter);
  elbo_vec(0) = INT_MIN;
  int iter;
  
  Rprintf("Finish variable initialization \n");
  
  // begin variational ICM-EM algorithm
  for(iter = 1; iter < maxIter; iter++){
    //E-step+ICM-step
    // Rprintf("Run E-step: \n");
    E_step_alignExp(Xf, Hf, Rf, Tm, Adjf, Muu, Eu,  Varu, 
                    alphaj, gammaj, zetaj, sigmaj, psij);
    // E_step_alignExp(Xf, Hf, Rf, Tm, Muu,   Eu,  Varu,alphaj, gammaj,zetaj, sigmaj, psij)
    
    // Rprintf("Finish E-step: \n");
    double Q0 = Q_fun( Xf, Hf, Rf, Tm, Muu, Eu, Varu, alphaj, gammaj, zetaj, sigmaj, psij); 
    //M-step
    // update alphaj
    // Rprintf("update alphaj: \n");
    alphaj = update_alphaj_alignExp( Xf, Hf, Rf, Tm, Eu, gammaj, zetaj, sigmaj);
    // double Q1=   Q_fun( Xf, Hf, Rf, Tm, Muu, Eu, Varu, alphaj, gammaj, zetaj, sigmaj, psij); 
    // Rprintf("dQ_alphaj = %4f \n", Q1-Q0);
    // update gammaj
    // Rprintf("update gammaj: \n");
    gammaj = update_gammaj_alignExp( Xf, Hf, Rf, Tm, Eu, alphaj, zetaj, sigmaj); 
    // double Q2=   Q_fun( Xf, Hf, Rf, Tm, Muu, Eu, Varu, alphaj, gammaj, zetaj, sigmaj, psij); 
    // Rprintf("dQ_gammaj = %4f \n", Q2-Q1);
    // update zetaj
    // Rprintf("update zetaj: \n");
    if(Tm_flag){
      zetaj = update_zetaj_alignExp( Xf, Hf, Rf, Tm, Eu, alphaj, gammaj, sigmaj); 
    }
    // double Q3=   Q_fun( Xf, Hf, Rf, Tm, Muu, Eu, Varu, alphaj, gammaj, zetaj, sigmaj, psij); 
    // Rprintf("dQ_zetaj = %4f \n", Q3-Q2);
    // update sigmaj
    // Rprintf("update sigmaj: \n");
    sigmaj = update_sigmaj_alignExp( Xf, Hf, Rf, Tm, Eu, Varu, alphaj, gammaj, zetaj); 
    // double Q4=   Q_fun( Xf, Hf, Rf, Tm, Muu, Eu, Varu, alphaj, gammaj, zetaj, sigmaj, psij); 
    // Rprintf("dQ_sigmaj = %4f \n", Q4-Q3);
    // update psij
    // Rprintf("update psij: \n");
    psij = update_psij_alignExp( Muu, Eu, Varu); 
    // double Q5 =   Q_fun( Xf, Hf, Rf, Tm, Muu, Eu, Varu, alphaj, gammaj, zetaj, sigmaj, psij); 
    // Rprintf("dQ_psij = %4f \n", Q5-Q4);
    
    elbo_vec(iter) = Q0;
    // output algorithm info.
    if(verbose){
      Rprintf("iter = %d, ELBO= %4f, dELBO=%4f \n", 
              iter +1, elbo_vec(iter), (elbo_vec(iter)  - elbo_vec(iter-1))/ abs(elbo_vec(iter-1)));
    }
    if(abs((elbo_vec(iter)  - elbo_vec(iter-1))/ elbo_vec(iter-1)) < epsELBO) break;
    
    
  }
  
  // output return value
  List resList = List::create(
    Rcpp::Named("gammaj") = gammaj,
    Rcpp::Named("alphaj") = alphaj,
    Rcpp::Named("zetaj") = zetaj,
    Rcpp::Named("ELBO") = elbo_vec(iter-1),
    Rcpp::Named("ELBO_seq") = elbo_vec.subvec(0, iter-1)
  );
  return(resList);
  
  
}


