// Date: 2022-09-08



#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
//#include<ctime>
#include <stdio.h>
#include <math.h>
//#include <thread>
//#include <mutex>


#define INT_MIN (-INT_MAX - 1)

using namespace std;
using namespace Rcpp;
using namespace arma;



/*
 * Auxiliary
 */
//' @keywords internal
//' @noRd
//' 



// diag((X - Nu) * diagmat(Lam_vec0).i() * (X-Nu)); sum(log(Lam_vec0))
void multi_det_Sk_variantCpp(const arma::mat& X, const arma::vec& Lam_vec0,  
                     const arma::rowvec Nu, 
                     double& logdSk, arma::vec& mSk){
  //int p = X.n_cols;
  int n = X.n_rows;
  
  logdSk = accu(log(Lam_vec0));
  mat tmp = (X - repmat(Nu, n, 1));
  mat X_tk = tmp.each_row()%trans(1/sqrt(Lam_vec0));
  tmp.reset();
  mSk = sum(X_tk % X_tk, 1);
}



// diag(W0^t* Cki * W0)
vec decomp(const mat& Cki, const mat& W0){
  vec s, tmp1;
  mat U, V, WC12;
  svd(U, s, V, Cki);
  WC12 = W0 * (U * diagmat(sqrt(s))); // p * q
  tmp1 = sum(WC12 % WC12, 1);
  return tmp1;
}  

// 
cube get_CubeInv(const cube& VarvI){
  cube Varv = VarvI;
  
  int i, n = VarvI.n_slices;
  for(i = 0; i < n; ++i){
    Varv.slice(i) = VarvI.slice(i).i();
  }
  
  return(Varv);
}
//
// update alpha and beta
vec update_alpha1(const arma::field<sp_mat>& Xf, const vec& tvec, const arma::field<mat>& Ezv, 
                  const mat& W, const mat& Lam, const vec& beta){
  
  int r, n,  q= Ezv(0).n_cols, M = Xf.n_elem, p = Xf(0).n_cols;
  mat WtLW(q,q, fill::zeros);
  rowvec nLamI(p, fill::zeros), tmpVec(q, fill::zeros);
  for(r=0; r < M; ++r){
    n = Xf(r).n_rows;
    nLamI += n / Lam.row(r);
    tmpVec += sum(Xf(r)- (Ezv(r) + tvec(r)* repmat(beta.t(),n, 1))* W.t()) * diagmat(1.0 /Lam.row(r)) * W;
    
  }
  WtLW = (W.t() % repmat(nLamI, q, 1)) * W;
  
  vec alpha = WtLW.i() * tmpVec.t();
  
  return alpha;
}

vec update_beta1(const arma::field<sp_mat>& Xf, const vec& tvec, const arma::field<mat>& Ezv, 
                 const mat& W, const mat& Lam, const vec& alpha){
  
  int r, n,  q= Ezv(0).n_cols, M = Xf.n_elem, p = Xf(0).n_cols;
  mat WtLW(q,q, fill::zeros);
  rowvec nLamI(p, fill::zeros), tmpVec(q, fill::zeros);
  for(r=0; r < M; ++r){
    n = Xf(r).n_rows;
    nLamI += n * tvec(r)*tvec(r)  / Lam.row(r);
    tmpVec += tvec(r) *  sum(Xf(r)- (Ezv(r) +repmat(alpha.t(),n, 1))* W.t())* diagmat(1/Lam.row(r)) * W ;
    
  }
  WtLW = (W.t() % repmat(nLamI, q, 1)) * W;
  
  vec beta = WtLW.i() * tmpVec.t();
  
  return beta;
}
// update Mu0
mat update_Mu1(const field<mat>& Ez){
  
  int r,  M=Ez.n_elem, q = Ez(0).n_cols;
  mat Mu(M, q);
  for (r=0; r <M; ++r){
    
    Mu.row(r) = mean(Ez(r));
  }
  return(Mu);
}

// update Sigma0
cube update_Sigma1(const mat& Mu, const field<mat>& Ez, const field<mat>& Varz, const bool& Sigma_diag){
  int r, n, M=Ez.n_elem, q = Ez(0).n_cols;
  cube Sigma0(q, q, M);
  mat Smat;
  for(r=0; r < M; ++r){
    Smat = zeros(q, q);
    n = Ez(r).n_rows;
    Smat = (Ez(r).each_row()- Mu.row(r)).t() * (Ez(r).each_row()- Mu.row(r)) /n +   Varz(r); 
    if(Sigma_diag){
      Sigma0.slice(r) = diagmat(Smat) ;
    }else{
      Sigma0.slice(r) = Smat;
    }
    
  }
  return(Sigma0);
}

// update Psi, the covariance of Vmi;
cube update_Psi1(const field<mat>& Ev, const field<mat>& Muv, const field<mat>& Varv, const bool& Psi_diag){
  
  int m,  n, M = Ev.n_elem, q = Ev(0).n_cols;
  cube Psi(q, q, M, fill::zeros);
  for(m = 0; m < M; ++m){
    n = Ev(m).n_rows;
    Psi.slice(m) = (Ev(m) - Muv(m)).t() * (Ev(m) - Muv(m))/n +  Varv(m);
    if(Psi_diag){
      Psi.slice(m) = diagmat(Psi.slice(m));
    }
  }
  
  return(Psi);
}


// 
mat update_nv1(const field<sp_mat>& Xf,  const mat&W,
               const field<mat>& Ev){
  
  int  m, M=Xf.n_elem, p = Xf(0).n_cols;
  mat nv(M, p, fill::zeros);
  for(m=0; m<M; ++m){
    
    nv.row(m) = mean(Xf(m) -  Ev(m) * W.t());
  }
  
  return nv;
}


//
mat update_W1(const field<sp_mat>& Xf, const mat& nv, const mat& Lam, 
              const field<mat>& Ev, const field<mat>& Varv){ 
 
  
  int  r, j, n,  q= Ev(0).n_cols, p = Xf(0).n_cols, M=Xf.n_elem;
  mat  A_w(p,q, fill::zeros);
  
  for(r=0; r<M; ++r){
    n = Xf(r).n_rows;
    mat tmpMat = (Xf(r)- repmat(nv.row(r), n, 1));
    A_w += trans(tmpMat.each_row() % (1.0/Lam.row(r))) * Ev(r);
  }
  
  
  cube B_arr(q,q,M, fill::zeros);
  for(r=0; r<M; ++r){
    n = Xf(r).n_rows;
    B_arr.slice(r) += Ev(r).t()*  Ev(r) +
      n* Varv(r); // There is an omitted mulitiplicative object n that causes the unconvergence of algorithm!!
    
    
  }
  mat tmpMat2;
  mat W(p, q, fill::zeros);
  for(j=0; j<p; ++j){
    tmpMat2 = zeros(q,q);
    for(r=0; r<M; ++r){
      tmpMat2 += 1.0/Lam(r,j)*B_arr.slice(r);
      
    }
    W.row(j) = A_w.row(j)* tmpMat2.i();
  }
  
  return W;
}

// update Lambda0
mat update_Lam1(const field<sp_mat>& Xf, const mat& nv,  const mat& W, 
              const field<mat>& Ev,const field<mat>& Varv,  const bool& homo = false){
  
  int r, n,  p = Xf(0).n_cols, M=Xf.n_elem;
  vec Lsum(p,fill::zeros);
  mat Lam(M, p), tmpX;
  for(r=0; r< M; ++r){
    n = Xf(r).n_rows;
    Lsum = zeros(p,1);
    tmpX = (Xf(r) - repmat(nv.row(r), n, 1)- Ev(r) * W.t());
    Lsum = trans(sum(tmpX % tmpX ))  + n* decomp(Varv(r), W);

    // Lsum +=  n*decomp(Varzv(r), W); // fast SVD method

    if(homo){
      Lam.row(r) = mean(Lsum)* ones(1, p) / (n*1.0);
    }else{
      Lam.row(r) = Lsum.t()/(n*1.0);
    }
  }
  // replace little values to increase the stability of the algorithm.
  uvec id1 =  find(Lam <1e-7);
  Lam(id1) = 1e-7*ones(id1.n_elem, 1);
  
  return Lam; 
}

//
arma::mat get_Vmean(const arma::mat& V, const arma::sp_mat& Adj){
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

//  diag((V-Muk)*SrkI(V-Muk))
void multi_det_Sk_embedCpp(const arma::mat& V,  const arma::mat& SrkI, 
                           const arma::rowvec Muk, 
                           double& logdSk, arma::vec& mSk){
  //int p = X.n_cols;
  int n = V.n_rows;
  // int p = X.n_cols;
  // // mSk = zeros(n);
  // S2k = zeros(p);
  // dSk = 0;
  
  mat WC12,  tmp2;
  vec tmp1, s, tmp3;
  mat U, V1, X_tk;
  
  svd(U, s, V1, SrkI);
  
  logdSk = accu(log(s)); // -log(|Srk|)
  X_tk = (V - repmat(Muk, n, 1)) * (U * diagmat(sqrt(s)));  // 
  mSk = sum(X_tk % X_tk, 1);
}


//
double Q_fun(const field<sp_mat>& Xf, const mat& nv0,  
             const field<mat>& Ev, const arma::field<mat>& Varv, 
             const field<mat>& Muv,  const mat& W0,
             const mat& Lam0, const cube& Psi0){
  
  double Q = 0; // tmp_scalar1 =0; // tmp_scalar2 = 0, tmp_scalar3= 0;
  int r,   M=Xf.n_elem,d= Ev(0).n_cols;
  
  
 
  // P(X)
  int nm;
  double logdSk=0.0;
  double logPx = 0.0;
  for(r=0; r< M; ++r){
    vec mSk;
    nm = Xf(r).n_rows;
    multi_det_Sk_variantCpp(Xf(r)-Ev(r)*W0.t(), Lam0.row(r).t(), nv0.row(r), // Use SVD to speed up
                            logdSk, mSk);
    //
    logPx += -0.5* nm* logdSk  -  0.5 * accu(mSk) - 0.5 * nm* trace(W0.t()* (W0.each_col() % (1.0/Lam0.row(r).t()))*Varv(r));
      
    // Rprintf("tmp_scalar1= %2f, logPost_z= %2f \n", tmp_scalar1, logPx);
  }
  
  // P(V)
  double logPost_v = 0.0;
  for(r=0; r< M; ++r){
      
      vec mSk;
      nm = Xf(r).n_rows;
      multi_det_Sk_embedCpp(Ev(r)-Muv(r), Psi0.slice(r).i(), zeros(1,d), // Use SVD to speed up: log P(V|y)
                            logdSk, mSk);
      logPost_v += 0.5* nm* logdSk  - 0.5 * accu(mSk) -  0.5 * nm* trace(Psi0.slice(r).i() * Varv(r));
    }
  
  Q = logPx + logPost_v;
  return Q;
}


//
double logLikem(const mat& Xm, const mat& bSm, const rowvec& vm, const mat& Wm, const rowvec& Mum,
                const mat& W, const mat& Muv){
  
  int i, n = Xm.n_rows, p = Xm.n_cols;
  
  double dbSm = log(det(bSm)), loglike=0.0;
  for(i = 0; i<n; ++i){
    loglike -= p* log(2* datum::pi) - dbSm + as_scalar((Xm.row(i) - vm-Muv.row(i)*W.t() + Mum* Wm.t()) * bSm * 
      trans(Xm.row(i) - vm-(Muv.row(i))*W.t() + Mum* Wm.t())) ;
  }
  return 0.5*loglike;
}



// Run ICM-step and E-step in a function for efficient computation.
void runICM_sp1(const arma::field<sp_mat>& Xf, const arma::mat& nv0,
                 const arma::mat& W0,  const arma::mat& Lam0,  const cube& Psi0, 
                 const arma::field<sp_mat>& Adjf,  
                 double& Q1,  field<mat>& Muv,
                  field<mat>& Ev, field<mat>& Varv){
  
 
  // basic info.
  int m, M = Xf.n_elem, q = W0.n_cols;
  // int p = Xf(0).n_cols;
  int n;
  
  // --------------- cached objects used for prediction in ICM step and posterior in E-step ----------
  field<mat> WtSW(M), XSW(M); // q*q, n*q
  
  Q1 = 0.0;
  // evaluate energy of x, Ux
  mat Cr_AI, Cr_S, WtLrW, Lam0_rW0, WtAmIW;
  mat  tmpMat;
  double tmp_value;
  // Wmt *Sm* W, Wt*Sm*W, Wmt *Sm*Wm, Wmt*Sm*X, Wt*Sm*X
  //Rprintf("cache some objects! \n");
  for(m = 0; m < M; ++m){
    
    Muv(m) = get_Vmean(Ev(m),  Adjf(m)); // update Muv
    // evaluate energy of x, Ux
    //Rprintf("m = %d ! \n", m);
    n = Xf(m).n_rows; // tmp object
    Lam0_rW0 = trans(repmat(1.0/ Lam0.row(m), q, 1)) % W0;
    WtLrW = W0.t() * Lam0_rW0; // O(p^2 q)
    Cr_AI = inv(Psi0.slice(m).i() +  WtLrW);
    
    // Rprintf("m = %d ! \n", m);
    WtAmIW = WtLrW - WtLrW*Cr_AI * WtLrW;
    // x_mi - W*(alpha + beta t_m)
    // mat XAW =  (Xf(m).each_row() - nv0.row(m)) *  Lam0_rW0 *(eye(q,q) - Cr_AI* WtLrW);
    mat XAW =  (Xf(m) - repmat(nv0.row(m), n, 1)) *  Lam0_rW0 *(eye(q,q) - Cr_AI* WtLrW);
    //Rprintf("m = %d ! \n", m);
    
    // save some objects for posterior of v
    Ev(m) = (XAW - Muv(m) * WtAmIW) *Psi0.slice(m) + Muv(m); // 
   tmpMat = Psi0.slice(m) - Psi0.slice(m)* WtAmIW *Psi0.slice(m);
   
    Varv(m) =  tmpMat;
   //Rprintf("detm = %4f\n", det(tmpMat));
   tmp_value = det(tmpMat);
   Q1 += 0.5*n*log(tmp_value +(tmp_value==0));
    
  }
  
}  



// Assume E(X) = 0.
// [[Rcpp::export]]
Rcpp:: List profast_g_cpp(const Rcpp::List& Xlist, const Rcpp::List& Adjlist, const arma::mat& nu_int,
                                const arma::mat& W_int, 
                                const arma::mat& Lam_int, const arma::cube& Psi_int, const Rcpp::List& EvList,
                                const int& maxIter, const double& epsLogLik, const bool& verbose,
                                const bool& homo = false, const bool& Psi_diag=false){
  
  // homo denotes error covariance is a scalar * identity matrix.
  //bool use_v = false;
  // basic info
  int r,M = Xlist.length(); // get the number of data source
  int q= W_int.n_cols;
  // Rprintf("Step into cpp file \n");
  // transfer list to field.
  field<sp_mat> Xf(M);
  field<sp_mat> Adjf(M);
  // Posterior of y, z, v, and u
  field<mat> Ev(M), Varv(M); // require to intialize values not only shape!!!
  int n;
  for(r=0; r < M; ++r){ 
    sp_mat Xtmp = Xlist[r]; // enforce to become a matrix.
    Xf(r) = Xtmp;
    sp_mat Adjtmp = Adjlist[r];
    Adjf(r) = Adjtmp;
    
    n = Xf(r).n_rows; // tmp object
    
    mat vmat = EvList[r];
    Ev(r) =  vmat;
    Varv(r) = eye(q,q);
  }
  
  //Rprintf("Step into cpp file 2 \n");
  
  // intialize the model parameters
  mat W0(W_int), Lam0(Lam_int), nv0(nu_int);
  cube Psi0(Psi_int);
 
 // Rprintf("Step into cpp file 3 \n");
  vec loglik(maxIter);
  loglik(0) = INT_MIN;
  // vec Qvec(loglik);
  
  // pseudo obserbed loglikelihood.
  double Q1;
  int iter;
  
  // create posterior expectation and covariance of z, v.
  field<mat> Muv(Ev); // Vf = Ev
  
  
  Rprintf("Finish variable intialization \n");
  
  double Q2 = 0;
  // begin ICM-EM algorithm
  for(iter = 1; iter < maxIter; iter++){
    
    // predict V.
    Rprintf("Satrt ICM and E-step! \n");
    runICM_sp1(Xf, nv0,W0, Lam0,  Psi0, Adjf,
               Q1,  Muv, Ev, Varv);
    
    
    
    
    Rprintf("Finish ICM and E-step! \n");
    
    
    Q2 = Q_fun( Xf, nv0, Ev, Varv, Muv, W0, Lam0, Psi0);
    //Rprintf("Q2= %4f \n", Q2);
    
    loglik(iter) = Q2+ Q1; // obtain the pseudo observed log-likelihood.
    
    // update Psi, the covariance of v
    //Rprintf("Run VB-M-step Psi! \n");
    Psi0 = update_Psi1(Ev, Muv, Varv, Psi_diag);
    // double Q3 =  Q_fun( Xf, nv0, Ev, Varv, Muv, W0, Lam0, Psi0);
    // Rprintf("dQ_Psi= %4f \n", Q3-Q2);
    
    
    // update num
    nv0 = update_nv1(Xf, W0, Ev);
    // double Q4 =  Q_fun( Xf, nv0, Ev, Varv, Muv, W0, Lam0, Psi0);
    // Rprintf("dQ_nv= %4f \n", Q4-Q3);
    
    
    W0 = update_W1(Xf, nv0, Lam0, Ev, Varv);
    // double Q5 =  Q_fun( Xf, nv0, Ev, Varv, Muv, W0, Lam0, Psi0);
    // Rprintf("dQ_W= %4f \n", Q5-Q4);
    
    
    
    // update  Lambda
    //Rprintf("Run VB-M-step Lam0! \n");
    Lam0 = update_Lam1(Xf, nv0,  W0, Ev,Varv, homo);
    // double Q6 = Q_fun( Xf, nv0, Ev, Varv, Muv, W0, Lam0, Psi0);
    // Rprintf("dQ_Lam= %4f \n", Q6-Q5);
    
    
    
    // output algorithm info.
    if(verbose){
      Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", 
              iter +1, loglik(iter), (loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
    }
    if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < epsLogLik) break;
    // if(abs(Q  - tmp_Q) < epsLogLik) break;
    
  }
  
  
  
  // output return value
  List resList = List::create(
    Rcpp::Named("hV") = Ev,
    Rcpp::Named("nu") = nv0,
    Rcpp::Named("Psi") = Psi0,
    Rcpp::Named("W") = W0,
    Rcpp::Named("Lam") = Lam0,
    Rcpp::Named("loglik") = loglik(iter-1),
    Rcpp::Named("dLogLik") = loglik(iter-1)  - loglik(iter-2),
    Rcpp::Named("loglik_seq") = loglik.subvec(0, iter-1)
  );
  return(resList);
}





// update Psi, the covariance of Vmi;
cube update_count_Psi(const field<mat>& Ev, const field<mat>& Muv, const field<mat>& Varv, const bool& Psi_diag){
  
  int m,  n, M = Ev.n_elem, q = Ev(0).n_cols;
  cube Psi(q, q, M, fill::zeros);
  for(m = 0; m < M; ++m){
    n = Ev(m).n_rows;
    Psi.slice(m) = (Ev(m) - Muv(m)).t() * (Ev(m) - Muv(m))/n +  Varv(m);
    if(Psi_diag){
      Psi.slice(m) = diagmat(Psi.slice(m));
    }
  }
  
  
  return(Psi);
}

// 
mat update_count_nv(const field<mat>& Ef, const mat&W, const field<mat>& Ev){
  
  int  m, M=Ef.n_elem, p = Ef(0).n_cols;
  mat nv(M, p, fill::zeros);
  for(m=0; m<M; ++m){
    
    nv.row(m) = mean(Ef(m) - Ev(m) * W.t());
  }
  
  return nv;
}

//
mat update_count_W(const field<mat>& Xf, const mat& nv, const mat& Lam, 
                    const field<mat>& Ev, const field<mat>& Varv){ 
  // Ci is the posterior variance matrix of Z
  
  int  r, j, n,  q= Ev(0).n_cols, p = Xf(0).n_cols, M=Xf.n_elem;
  mat  A_w(p,q, fill::zeros);
  
  for(r=0; r<M; ++r){
    n = Xf(r).n_rows;
    mat tmpMat = (Xf(r)- repmat(nv.row(r), n, 1));
    A_w += trans(tmpMat.each_row() % (1.0/Lam.row(r))) * Ev(r);
  }
  
  
  cube B_arr(q,q,M, fill::zeros);
  for(r=0; r<M; ++r){
    n = Xf(r).n_rows;
    B_arr.slice(r) = Ev(r).t()*  Ev(r) + n* Varv(r); // There is an omitted mulitiplicative object n that causes the unconvergence of algorithm!!
  }
  mat tmpMat2;
  mat W(p, q, fill::zeros);
  for(j=0; j<p; ++j){
    tmpMat2 = zeros(q,q);
    for(r=0; r<M; ++r){
      tmpMat2 += 1.0/Lam(r,j)*B_arr.slice(r);
      
    }
    W.row(j) = A_w.row(j)* tmpMat2.i();
  }
  
  return W;
}

// update Lambda0
mat update_count_Lam(const field<mat>& Xf, const field<mat>& Vare, const mat& nv,  const mat& W, 
                     const field<mat>& Ev,const field<mat>& Varv,  const bool& homo = false){
  
  int r, n,  p = Xf(0).n_cols, M=Xf.n_elem;
  vec Lsum(p,fill::zeros);
  mat Lam(M, p), tmpX;
  for(r=0; r< M; ++r){
    n = Xf(r).n_rows;
    Lsum = zeros(p,1);
    tmpX = (Xf(r) - repmat(nv.row(r), n, 1)- Ev(r) * W.t());
    Lsum = trans(sum(tmpX % tmpX + Vare(r))) + n* decomp(Varv(r), W);
    
    
    if(homo){
      Lam.row(r) = mean(Lsum)* ones(1, p) / (n*1.0);
    }else{
      Lam.row(r) = Lsum.t()/(n*1.0);
    }
  }
  // replace little values to increase the stability of the algorithm.
  uvec id1 =  find(Lam <1e-7);
  Lam(id1) = 1e-7*ones(id1.n_elem, 1);
  
  
  return Lam; 
}




double calELBO(const arma::field<sp_mat>& Xf, const field<vec>& Af, const arma::mat& nv0, 
               const arma::mat& W0,  const arma::mat& Lam0, 
               const cube& Psi0, const field<mat>& Muv, const field<mat>& Ee, const field<mat>& Vare,
               const field<mat>& Ev, const field<mat>& Varv){
  int m, M=Xf.n_elem, p = Xf(0).n_cols, n_m, q= W0.n_cols;
  //  P(X|E)
  double logPX_E = 0.0;
  for(m = 0; m<M; ++m){
    logPX_E += accu(Xf(m) % Ee(m) - repmat(Af(m), 1, p) % exp(Ee(m)+ Vare(m)*0.5));
  }
  
  // P(E| Z, V)
  double logPE_ZV=0.0;
  mat tmp_mat1, tmp_mat2,WdL;
  for(m=0; m<M; ++m){
    n_m = Ee(m).n_rows;
    WdL = W0 / repmat(sqrt(Lam0.row(m).t()), 1, q);
    tmp_mat1 = (Ee(m) - repmat(nv0.row(m), n_m, 1) - Ev(m)*W0.t());
    logPE_ZV += accu(tmp_mat1% tmp_mat1  / repmat(Lam0.row(m), n_m,1)) + n_m *accu(log(Lam0.row(m)));
    logPE_ZV +=  accu(decomp(n_m*Varv(m), WdL)) + accu(Vare(m)/ repmat(Lam0.row(m), n_m, 1));
  }
  logPE_ZV = -0.5 * logPE_ZV;
  
  
  // P(V)
  double logPV = 0.0;
  for(m = 0; m<M; ++m){
    n_m = Ee(m).n_rows;
    logPV += n_m * log_det_sympd(Psi0.slice(m)) + accu(decomp(Psi0.slice(m).i(),Ev(m)-Muv(m))) + n_m * trace(Psi0.slice(m).i()* Varv(m));
  }
  logPV = -0.5 * logPV;
  
  // Entropy
  double entropy = 0.0;
  for(m=0; m <M; ++m){
    n_m = Ee(m).n_rows;
    entropy += accu(log(Vare(m))) +  n_m* log_det_sympd(Varv(m));
  }
  entropy = 0.5* entropy;
  
  
  return logPX_E + logPE_ZV + logPV + entropy;
}


void runEstep1(const field<sp_mat>& Xf, const field<vec>& Af, const field<sp_mat>& Adjf,  const mat& W0,
               const mat& nv0, const mat& Lam0, const cube& Psi0,
               field<mat>& Muv, field<mat>& Ev, field<mat> Varv, field<mat>& Ee, field<mat>& Vare){
  
  int r,  n, M = Xf.n_elem, p=Xf(0).n_cols, q = W0.n_cols;
  
  // E-step
  // update Ee and Vare
  // double QQ1 =  calELBO(Xf, Af, nv0, W0, Lam0,Psi0,  Muv, Ee, Vare, Ev, Varv);
  for(r = 0; r<M;  ++r){
    n = Xf(r).n_rows;
    mat tmp_mat =  Ev(r) * W0.t() + repmat(nv0.row(r), n, 1);
    Ee(r) = (Xf(r) - repmat(Af(r), 1, p) % exp(Ee(r)) % (1- Ee(r)) + repmat(1/Lam0.row(r), n,1) % tmp_mat) / 
      (repmat(Af(r), 1, p) % exp(Ee(r)) + repmat(1/Lam0.row(r), n,1) );
    Vare(r) = 1.0 / (repmat(Af(r), 1, p) % exp(Ee(r) ) + repmat(1/Lam0.row(r), n,1) );
  }
  // double QQ2 =  calELBO(Xf, Af, nv0, W0, Lam0,Psi0,  Muv, Ee, Vare, Ev, Varv);
  // Rprintf("dQ_Estep_E= %4f \n", QQ2-QQ1);
  
  
  // update Ev and Varv
  for(r = 0; r<M; ++r){
    n = Xf(r).n_rows;
    Varv(r) = inv_sympd(W0.t() * (W0 / repmat(Lam0.row(r).t(), 1, q)) + inv_sympd(Psi0.slice(r)));
    Ev(r) = ((Ee(r) - repmat(nv0.row(r), n, 1) ) * (W0/ repmat(Lam0.row(r).t(), 1, q))+
      Muv(r) * inv_sympd(Psi0.slice(r))) * (Varv(r)); // inv_sympd
  }
  // double QQ4 =  calELBO(Xf, Af, nv0, W0, Lam0,Psi0,  Muv, Ee, Vare, Ev, Varv);
  // Rprintf("dQ_Estep_V= %4f \n", QQ4-QQ2);
  
  // ICM-step
  for(r = 0; r<M; ++r){
    Muv(r) =  get_Vmean(Ev(r),  Adjf(r)); // update Muv
  }
}


// [[Rcpp::export]]
Rcpp::List profast_p_cpp(const Rcpp::List& Xlist, const Rcpp::List& AList, const Rcpp::List& Adjlist, 
                             const arma::mat& nv_int,const arma::mat& W_int,
                             const arma::mat& Lam_int, const arma::cube& Psi_int, 
                              const Rcpp::List& EvList,
                             const int& maxIter, const double& epsELBO, const bool& verbose,
                             const bool& homo = false, 
                             const bool& Psi_diag=true){
  
  // homo denotes error covariance is a scalar * identity matrix.
  
  //Rprintf("It is good \n");
  // basic info
  int r,M = Xlist.length(), p = W_int.n_rows; // get the number of data source
  int q= W_int.n_cols;
  // transfer list to field.
  field<sp_mat> Xf(M);
  field<sp_mat> Adjf(M);
  field<vec> Af(M);
  // Posterior of y, z, v, and u
  //cube Wm0(p, q_batch, M);
  field<mat> Ee(M), Vare(M),  Ev(M), Varv(M); // require to initialize values not only shape!!!
  for(r=0; r < M; ++r){ 
    //Rprintf("It is good 1.5 \n");
    sp_mat Xtmp = Xlist[r]; // enforce to become a matrix.
    Xf(r) = Xtmp;
    vec atmp = AList[r];
    Af(r) = atmp;
    
    sp_mat Adjtmp = Adjlist[r];
    Adjf(r) = Adjtmp;
    
    
    Ee(r) = log(1+Xf(r));
    Vare(r) = ones(Xf(r).n_rows, p);
    
    mat vmat = EvList[r];
    Ev(r) =  vmat;
    Varv(r) = eye(q,q);
  }
  //Rprintf("It is good2 \n");
  // Initialize the model parameters
  mat W0(W_int), Lam0(Lam_int), nv0(nv_int);
  cube Psi0(Psi_int);
  
  
  vec elbo_vec(maxIter);
  elbo_vec(0) = INT_MIN;
  int iter;
  
  // create posterior expectation and covariance of z, v.
  field<mat> Muv(Ev); // Vf = Ev
  
  Rprintf("Finish variable initialization \n");
  // Rprintf("It is good3 \n");
  // begin variational ICM-EM algorithm
  for(iter = 1; iter < maxIter; iter++){
    
    double Q1 = calELBO(Xf, Af, nv0, W0, Lam0,Psi0,  Muv, Ee, Vare, Ev, Varv);
    elbo_vec(iter) = Q1;
    Rprintf("Satrt ICM step! \n"); // update the variational parameters and V, Muv;
    runEstep1(Xf, Af, Adjf, W0,  nv0, Lam0, Psi0, 
               Muv, Ev, Varv, Ee, Vare);
    
    
    
    
    // update Psi, the covariance of v
    //Rprintf("Run VB-M-step Psi! \n");
    Psi0 = update_count_Psi(Ev, Muv, Varv, Psi_diag);
    // double Q5 =  calELBO(Xf, Af, nv0, W0, Lam0,Psi0,  Muv, Ee, Vare, Ev, Varv);
    // Rprintf("dQ_Psi= %4f \n", Q5-Q1);
    
    
    // update num
    nv0 = update_count_nv(Ee, W0,  Ev);
    // double Q6 =   calELBO(Xf, Af, nv0, W0, Lam0,Psi0,  Muv, Ee, Vare, Ev, Varv);
    // Rprintf("dQ_nv= %4f \n", Q6-Q5);
    
    // update W
    
    W0 = update_count_W(Ee, nv0,  Lam0, Ev, Varv);
    // double Q7 =  calELBO(Xf, Af, nv0, W0, Lam0,Psi0,  Muv, Ee, Vare, Ev, Varv);
    // Rprintf("dQ_W= %4f \n", Q7-Q6);
    
    
    // update  Lambda
    //Rprintf("Run VB-M-step Lam0! \n");
    Lam0 = update_count_Lam(Ee, Vare, nv0,W0, Ev,Varv, homo);
    // double Q8 =  calELBO(Xf, Af, nv0, W0, Lam0,Psi0,  Muv, Ee, Vare, Ev, Varv);
    // Rprintf("dQ_Lam= %4f \n", Q8-Q7);
    
    
    // output algorithm info.
    if(verbose){
      Rprintf("iter = %d, ELBO= %4f, dELBO=%4f \n", 
              iter +1, elbo_vec(iter), (elbo_vec(iter)  - elbo_vec(iter-1))/ abs(elbo_vec(iter-1)));
    }
    if(abs((elbo_vec(iter)  - elbo_vec(iter-1))/ elbo_vec(iter-1)) < epsELBO) break;
    
    
    
  }
  
  
  // output return value
  List resList = List::create(
    Rcpp::Named("hV") = Ev,
    Rcpp::Named("nu") = nv0,
    Rcpp::Named("Psi") = Psi0,
    Rcpp::Named("W") = W0,
    Rcpp::Named("Lam") = Lam0,
    Rcpp::Named("ELBO") = elbo_vec(iter-1),
    Rcpp::Named("ELBO_seq") = elbo_vec.subvec(0, iter-1)
  );
  return(resList);
  
}