// Date: 2023-02-09

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>



using namespace Rcpp;
using namespace arma;
using namespace std;


// correct for one gene

mat update_alpha_alignExps(const field<mat>& Xf, const field<mat>& Hf, const field<mat>& Rf, 
                           const mat& Tm, const field<mat>& Eu, const mat& gamma, const mat& zeta,
                           const mat& sigma){
  
  int r, n_r, M=Xf.n_elem, d = Rf(0).n_cols, p=Xf(0).n_cols;
  cube y_right(d, p, M);
  cube y_left(d,d, M, fill::zeros);
  mat alpha(d, p, fill::zeros);
  for(r=0; r<M; ++r){
    // Rprintf("update_alpha: r = %d\n", r);
    n_r = Rf(r).n_rows;
    y_left.slice(r) = Rf(r).t() * Rf(r);
    y_right.slice(r) = Rf(r).t() * (Xf(r)- Hf(r)*gamma  - Eu(r)-repmat(Tm.row(r) * zeta, n_r, 1));
  }
  int j;
  mat z_left(d,d, fill::zeros);
  vec z_right(d, fill::zeros);
  for(j = 0; j<p; ++j){
    z_left = zeros(d,d);
    z_right = zeros(d,1);
    //Rprintf("update_alpha: j = %d\n", j);
    for(r=0; r<M; ++r){
      
      z_left+= y_left.slice(r) / sigma(r,j);
      z_right+= y_right.slice(r).col(j) / sigma(r,j);
    }
    alpha.col(j) = inv_sympd(z_left) * z_right;
  }
  
  
  return alpha;
}

//
mat update_gamma_alignExps(const field<mat>& Xf, const field<mat>& Hf, const field<mat>& Rf, 
                             const mat& Tm, const field<mat>& Eu, const mat& alpha, const mat& zeta,
                             const mat& sigma){
    
    ///////////////////
    int r, n_r, M=Xf.n_elem, d = Hf(0).n_cols, p=Xf(0).n_cols;
    cube y_right(d, p, M);
    cube y_left(d,d, M, fill::zeros);
    mat gamma(d, p, fill::zeros);
    for(r=0; r<M; ++r){
      //Rprintf("update_gamma: r = %d\n", r);
      n_r = Rf(r).n_rows;
      y_left.slice(r) = Hf(r).t() * Hf(r);
      y_right.slice(r) = Hf(r).t() * (Xf(r)- Rf(r)*alpha  - Eu(r)-repmat(Tm.row(r) * zeta, n_r, 1));
    }
    int j;
    mat z_left(d,d, fill::zeros);
    vec z_right(d, fill::zeros);
    for(j = 0; j<p; ++j){
      z_left = zeros(d,d);
      z_right = zeros(d,1);
      //Rprintf("update_gamma: jr = %d\n", j);
      for(r=0; r<M; ++r){
        z_left+= y_left.slice(r) / sigma(r,j);
        z_right+= y_right.slice(r).col(j) / sigma(r,j);
      }
      gamma.col(j) = inv_sympd(z_left) * z_right;
    }
    
    return gamma;
}

mat update_zeta_alignExps(const field<mat>& Xf, const field<mat>& Hf, const field<mat>& Rf, 
                          const mat& Tm, const field<mat>& Eu, const mat& alpha, const mat& gamma,
                          const mat& sigma){
  
  // int r, M=Xf.n_elem, d = Tm.n_cols, n_r;
  // vec y_right(d);
  // mat y_left(d,d, fill::zeros);
  // for(r=0; r<M; ++r){
  //   n_r = Xf(r).n_rows;
  //   y_left += Tm.row(r).t() * Tm.row(r)* n_r /sigma(r);
  //   y_right += Tm.row(r).t() * accu(Xf(r)-Hf(r)*gamma - Rf(r)*alpha  - Eu(r))/sigma(r);
  // }
  
  
  ///////////////////
  int r, n_r, M=Xf.n_elem, d = Tm.n_cols, p=Xf(0).n_cols;
  cube y_right(d, p, M);
  cube y_left(d,d, M, fill::zeros);
  mat zeta(d, p, fill::zeros);
  for(r=0; r<M; ++r){
    n_r = Rf(r).n_rows;
    y_left.slice(r) = Tm.row(r).t() * Tm.row(r)* n_r;
    y_right.slice(r) =Tm.row(r).t()  * sum(Xf(r)- Hf(r)*gamma - Rf(r)*alpha  - Eu(r)); // d*1 x 1*p
  }
  int j;
  mat z_left(d,d, fill::zeros);
  vec z_right(d, fill::zeros);
  for(j = 0; j<p; ++j){
    z_left = zeros(d,d);
    z_right = zeros(d,1);
    for(r=0; r<M; ++r){
      z_left+= y_left.slice(r) / sigma(r,j);
      z_right+= y_right.slice(r).col(j) / sigma(r,j);
    }
    zeta.col(j) = inv_sympd(z_left) * z_right;
  }
  
  return zeta;
  
}

mat update_sigma_alignExps(const field<mat>& Xf, const field<mat>& Hf, const field<mat>& Rf,
                           const mat& Tm, const field<mat>& Eu, const field<mat>& Varu,
                           const mat& alpha, const mat& gamma,const mat& zeta){
  int r, M=Xf.n_elem, n_r, p = Xf(0).n_cols;
  mat sigma(M, p);
  
  for(r=0; r<M; ++r){
    n_r = Xf(r).n_rows;
    mat  tmpmat = Xf(r)-Hf(r)*gamma - Rf(r)*alpha  - Eu(r)- repmat(Tm.row(r) * zeta, n_r, 1); // n_r * p
    sigma.row(r) = sum(tmpmat % tmpmat + Varu(r))  / n_r;
  }
  
  return sigma;
}

mat update_psi_alignExps(const field<mat>& Muu, const field<mat>& Eu, const field<mat>& Varu){
  
  int r, M=Muu.n_elem, n_r, p = Muu(0).n_cols;
  mat psi(M, p);
  for(r=0; r<M; ++r){
    n_r = Muu(r).n_rows;
    psi.row(r) = sum((Muu(r)-Eu(r)) % (Muu(r)-Eu(r)) + Varu(r)) / n_r;
  }
  
  return(psi);
}

//
arma::mat get_Vmean_alignExps(const arma::mat& V, const arma::sp_mat& Adj){
    int i, n = V.n_rows, q= V.n_cols;
    vec m(n);
    mat Uv(n, q, fill::zeros);
    for (i = 0; i < n; i++)
    {
      arma::vec col(Adj.col(i)); // the class label of neighbors of i-th sample.
      uvec q1 = find(col > 0);
      // cout<<q1.n_rows<<endl;
      if( q1.n_rows>0){
        Uv.row(i) = mean(V.rows(q1));
      }
      
    }
    return Uv;
}

void E_step_alignExps(const field<mat>& Xf, const field<mat>& Hf, const field<mat>& Rf, 
                     const mat& Tm, const field<sp_mat>& Adjf, field<mat>& Muu,  field<mat>& Eu,  field<mat>& Varu,
                     const mat& alpha, const mat& gamma,const mat& zeta,
                     const mat& sigma, const mat& psi){
  // psi(M*p); sigma(M*p)
  
  int r, M=Muu.n_elem, n_r;
  field<mat> Mu_x(Xf); // n_r* p
  for(r=0; r<M; ++r){
    
    n_r = Xf(r).n_rows;
    Mu_x(r) = Muu(r) + Rf(r) * alpha + Hf(r) * gamma + repmat(Tm.row(r) * zeta, n_r, 1); // 1*d x d * p
    
    Eu(r) = Muu(r) + (Xf(r) - Mu_x(r))% repmat(psi.row(r)/(psi.row(r) + sigma.row(r)), n_r, 1);
    Varu(r) = repmat(psi.row(r) - psi.row(r) % psi.row(r)/(psi.row(r) + sigma.row(r)), n_r, 1);
    
  }
  // ICM-step:
    // ICM-step
  for(r = 0; r<M; ++r){
    Muu(r) =  get_Vmean_alignExps(Eu(r),  Adjf(r)); // update Muv
  }
  
}

// Q-function
double Q_fun_expsj(const field<vec>& Xf, const field<mat>& Hf, const field<mat>& Rf, 
             const mat& Tm, const field<vec>& Muu,  const field<vec>& Eu,  const field<vec>& Varu,
             const vec& alpha, const vec& gamma,const vec& zeta,
             const vec& sigma, const vec& psi){
  
  int r, M=Xf.n_elem, n_r;
  double logPx = 0.0;
  for(r=0; r<M; ++r){
    n_r = Xf(r).n_rows;
    vec  tmpvec = Xf(r)-Hf(r)*gamma - Rf(r)*alpha - Eu(r)- as_scalar(Tm.row(r) * zeta);
    logPx += -0.5* (n_r * log(sigma(r)) + accu(tmpvec % tmpvec + Varu(r))  / sigma(r));
  }
  
  double logPu = 0.0;
  double entropy = 0.0;
  for(r=0; r<M; ++r){
    n_r = Xf(r).n_rows;
    logPu += -0.5* (n_r * log(psi(r)) + accu((Muu(r)-Eu(r)) % (Muu(r)-Eu(r)) + Varu(r))  / psi(r));
    entropy += -0.5 * accu(log(Varu(r)));
  }
  
  
  return logPx + logPu - entropy;
}

//
double Q_fun_exps(const field<mat>& Xf, const field<mat>& Hf, const field<mat>& Rf, 
              const mat& Tm, const field<mat>& Muu,  const field<mat>& Eu,  const field<mat>& Varu,
              const mat& alpha, const mat& gamma,const mat& zeta,
              const mat& sigma, const mat& psi){
  
  int r, M=Xf.n_elem, n_r;
  double logPx = 0.0;
  for(r=0; r<M; ++r){
    n_r = Xf(r).n_rows;
    mat  tmpmat = Xf(r)-Hf(r)*gamma - Rf(r)*alpha - Eu(r)- repmat(Tm.row(r) * zeta,n_r, 1); // n_r * p
    logPx += -0.5* accu(n_r * log(sigma.row(r)) + sum(tmpmat % tmpmat + Varu(r))  / sigma.row(r));
  }
  
  double logPu = 0.0;
  double entropy = 0.0;
  for(r=0; r<M; ++r){
    n_r = Xf(r).n_rows;
    logPu += -0.5* accu(n_r * log(psi.row(r)) + sum((Muu(r)-Eu(r)) % (Muu(r)-Eu(r)) + Varu(r))  / psi.row(r));
    entropy += -0.5 * accu(log(Varu(r)));
  }
  
  
  return logPx + logPu - entropy;
}

// simultaneous correction for a group of genes
// [[Rcpp::export]]
Rcpp::List correct_genes(const Rcpp::List& XList, const Rcpp::List& RList,
                         const Rcpp::List& HList, const arma::mat& Tm, const Rcpp::List& Adjlist, 
                         const arma::mat& sigma_int, const arma::mat& psi_int, const arma::mat& alpha_int, 
                         const arma::mat& gamma_int, const arma::mat& zeta_int, const int& maxIter,
                         const double& epsELBO, const bool& verbose){
  // Xlist's component is a vector, denoting one gene's expression.
  int r,M = XList.length();
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
  mat sigma(sigma_int), psi(psi_int), alpha(alpha_int), gamma(gamma_int), zeta(zeta_int);
  
  field<sp_mat> Adjf(M);
  field<mat> Hf(M), Rf(M);
  field<mat> Eu(M), Varu(M), Xf(M); // require to initialize values not only shape!!!
  int  p=1;
  for(r=0; r < M; ++r){ 
    mat Xtmp = XList[r]; // enforce to become a matrix.
    Xf(r) = Xtmp;
    if(r == 0){
      p = Xf(0).n_cols;
    }
    
    mat Htmp = HList[r];
    Hf(r) = Htmp;
    mat Rtmp = RList[r];
    Rf(r) = Rtmp;
    
    sp_mat Adjtmp = Adjlist[r];
    Adjf(r) = Adjtmp;
    
    
    Eu(r) = zeros(Xf(r).n_rows, p);
    Varu(r) = ones(Xf(r).n_rows, p);
  }
  field<mat> Muu(Eu);
  vec elbo_vec(maxIter);
  elbo_vec(0) = -1e15;
  int iter;
  
  Rprintf("Finish variable initialization \n");
  double Q0 = 0.0;
  // begin variational ICM-EM algorithm
  for(iter = 1; iter < maxIter; iter++){
    //E-step+ICM-step
    //Rprintf("Run E-step: \n");
    E_step_alignExps(Xf, Hf, Rf, Tm, Adjf, Muu, Eu,  Varu, 
                    alpha, gamma, zeta, sigma, psi);
    // E_step_alignExps(Xf, Hf, Rf, Tm, Muu,   Eu,  Varu,alpha, gamma,zeta, sigma, psi)
    
    //Rprintf("Finish E-step: \n");
    // double Q0 = Q_fun_expsj(Xfj, Hf, Rf, Tm, Muuj, Euj, Varuj, alphaj, gammaj, zetaj, sigmaj, psij); 
    Q0 = Q_fun_exps(Xf, Hf, Rf, Tm, Muu, Eu, Varu, alpha, gamma, zeta, sigma, psi); 
    //M-step
    // update alpha
    //Rprintf("update alpha: \n");
    alpha = update_alpha_alignExps( Xf, Hf, Rf, Tm, Eu, gamma, zeta, sigma);
    // double Q1=   Q_fun_exps( Xf, Hf, Rf, Tm, Muu, Eu, Varu, alpha, gamma, zeta, sigma, psi);
    // Rprintf("dQ_alpha = %4f \n", Q1-Q0);
    // update gamma
    //Rprintf("update gamma: \n");
    gamma = update_gamma_alignExps( Xf, Hf, Rf, Tm, Eu, alpha, zeta, sigma); 
    // double Q2=   Q_fun_exps( Xf, Hf, Rf, Tm, Muu, Eu, Varu, alpha, gamma, zeta, sigma, psi);
    // Rprintf("dQ_gamma = %4f \n", Q2-Q1);
    // update zeta
    //Rprintf("update zeta: \n");
    if(Tm_flag){
      zeta = update_zeta_alignExps( Xf, Hf, Rf, Tm, Eu, alpha, gamma, sigma); 
    }
    // double Q3=   Q_fun_exps( Xf, Hf, Rf, Tm, Muu, Eu, Varu, alpha, gamma, zeta, sigma, psi);
    // Rprintf("dQ_zeta = %4f \n", Q3-Q2);
    // update sigma
    //Rprintf("update sigma: \n");
    sigma = update_sigma_alignExps( Xf, Hf, Rf, Tm, Eu, Varu, alpha, gamma, zeta); 
    // double Q4=   Q_fun_exps( Xf, Hf, Rf, Tm, Muu, Eu, Varu, alpha, gamma, zeta, sigma, psi);
    // Rprintf("dQ_sigma = %4f \n", Q4-Q3);
    // update psi
    //Rprintf("update psi: \n");
    psi = update_psi_alignExps( Muu, Eu, Varu); 
    // double Q5 =   Q_fun_exps( Xf, Hf, Rf, Tm, Muu, Eu, Varu, alpha, gamma, zeta, sigma, psi);
    // Rprintf("dQ_psi = %4f \n", Q5-Q4);
    
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
    Rcpp::Named("gamma") = gamma,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("zeta") = zeta,
    Rcpp::Named("ELBO") = elbo_vec(iter-1),
    Rcpp::Named("ELBO_seq") = elbo_vec.subvec(0, iter-1)
  );
  return(resList);
  
  
}



