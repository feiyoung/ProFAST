// Rewrite iSC-MEB
// Date: 2023-02-14


#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
#include <stdio.h>
#include <math.h>
#include <thread>
#include <mutex>

#include "DR_removebatch.h"

#define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;



// The following functions are moved to util.cpp as common functions to be called.
sp_mat get_spNbs_embed(ivec y, const sp_mat& Adj) {   // ivec是索引型向量
  // row is for pixel.
  //output a sparse matrix, i-th row contains labels of neighbor_i.
  // Make const iterator
  arma::sp_mat::const_iterator start = Adj.begin(); //构造一个sp_mat的常数迭代器,常数迭代器只可读，不可写，实现对矩阵按照列对每个非零元素进行访问。
  //arma::sp_mat::const_iterator end   = Adj.end();
  
  // Calculate number of nonzero points
  //int n = std::distance(start, end);
  int n = Adj.n_nonzero; // 计算Adj的所有非零元的个数
  //cout << "n=" << n << endl;
  //cout << "n=" << Adj.n_nonzero << endl;
  
  sp_mat spNbs(y.n_elem, y.n_elem);    // neiborhood state matrix, matched with Adj.
  
  
  arma::sp_mat::const_iterator it = start; // Note spNbs is not a symmetric matrix, the nonzero in i-th row is the class label of sample i.
  for(int i = 0; i < n; ++i)
  {
    //temp(0) = it.row();
    //temp(1) = it.col();
    spNbs(it.row(), it.col()) = y(it.col()); // it只自加非零元个数次，得到每个i对应的邻居的状态
    ++it; // increment
  }
  
  return spNbs.t(); // return the class label of neighbor matrix, i-th column is the neighbor label of sample i
}



arma::mat calYenergy2D_sp_embed(const arma::ivec& y, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
  // Calculate the value of energy function of y, which is equal to negative logliklihood up to a constant
  int n = y.n_rows;
  arma::sp_mat spNbs_t = get_spNbs_embed(y, Adj); // transform spNbs to iterate by column.
  arma::mat Uy(n, K);
  double n_sameS;
  int i, k, nn;
  for (k = 0; k < K; k++)
  {
    for (i = 0; i < n; i++)
    {
      arma::sp_mat col(spNbs_t.col(i)); // the class label of neighbors of i-th sample.
      n_sameS = 0;
      
      nn = col.n_nonzero; // the number of neighbors of i-th sample
      for (arma::sp_mat::iterator j = col.begin(); j != col.end(); ++j) {
        n_sameS += ((*j) == (k+1));
        
      }
      Uy(i, k) = alpha(k) + beta * (nn - n_sameS); // /2
      
      
    }
  }
  
  arma::mat C_mat = normalise(exp(-Uy), 1, 1); // pseudo likelihood of Y.
  Uy = -log(C_mat); // normalized Uy, this is the energy of y.
  return Uy;
  
}

//
double obj_beta_embed(const field<ivec>& yf, const field<mat>& Rf,
                      const arma::field<sp_mat>& Adjf, int K, const arma::vec alpha, const double beta)	{
  int r, M = yf.n_elem;
  double objval = 0;
  for(r=0; r< M; ++r){
    mat Uy1 = calYenergy2D_sp_embed(yf(r), Adjf(r), K, alpha, beta); // Uy was normalized, so there is no need to normalized Uy.
    objval += -accu(Rf(r) % Uy1);
  }
  
  return objval;
}

//
double objr_beta_embed(const ivec& y, const mat& R,
                       const sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
  
  double objval = 0;
  mat Uy1 = calYenergy2D_sp_embed(y, Adj, K, alpha, beta); // Uy was normalized, so there is no need to normalized Uy.
  objval = -accu(R % Uy1);
  
  
  return objval;
}
// diag(W0^t* Cki * W0)
vec decomp_embed(const mat& Cki, const mat& W0){
  vec s, tmp1;
  mat U, V, WC12;
  svd(U, s, V, Cki);
  WC12 = W0 * (U * diagmat(sqrt(s))); // p * q
  tmp1 = sum(WC12 % WC12, 1);
  return tmp1;
}




void multi_det_Sk_embedCpp2(const arma::mat& V,  const arma::mat& SrkI, 
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



// update Mu0
mat update_embed_Mu1(const arma::field<mat>& Rf, const arma::field<cube>& Ez, 
                     const mat& Nmat){
  
  int r, k, K= Rf(0).n_cols, q= Ez(0).n_cols, r_max = Rf.n_elem;
  mat Mu(q, K, fill::zeros);
  vec b_muk;
  
  for(k=0;k < K ; ++k){
    b_muk = zeros(q,1);
    for(r=0; r< r_max; ++r){
      b_muk += trans(Ez(r).slice(k)) * Rf(r).col(k);
      
    }
    Mu.col(k) = b_muk/ accu(Nmat.col(k));
  }
  return Mu.t();
}

// update Sigma0
cube update_embed_Sigma1(const field<mat>& Rf, const field<cube>& Ez, const field<cube>& Varz,
                         const mat& Mu, const mat& Nmat,
                         const bool& homoClust=false, const bool& Sigma_diag=false){
  
  // homoCluste implies whether Sigma_k = Sigma_l, for all k,l.
  int r,k, K= Mu.n_rows, q=Mu.n_cols, r_max=Rf.n_elem;
  cube Sigma0(q,q , K);
  mat Smat(q, q, fill::zeros);
  // Rprintf("r_max = %d, K= %d, q= %d", r_max, K, q);
  if(homoClust){// Sigma is not related to  k.
    Smat = zeros<mat>(q,q);
    for(k = 0; k<K; ++k){  
      for(r=0; r< r_max; ++r){
        Smat += (trans(Ez(r).slice(k) - repmat(Mu.row(k), Rf(r).n_rows, 1)) % trans(repmat(Rf(r).col(k), 1, q))) * 
          (Ez(r).slice(k) - repmat(Mu.row(k), Rf(r).n_rows, 1)) + Nmat(r,k)*  (Varz(r).slice(k)); //yangyi n-by-n matrix to n-by-q matrix
      }
    }
    for(k = 0; k<K; ++k){ 
      if(Sigma_diag){
        Sigma0.slice(k) = diagmat(Smat) / accu(Nmat);
      }else{
        Sigma0.slice(k) = Smat /  accu(Nmat);
      }
    }
  }else{// Sigma is  related to  k.
    for(k = 0; k<K; ++k){  
      Smat = zeros<mat>(q,q);
      for(r=0; r< r_max; ++r){
        Smat += (trans(Ez(r).slice(k) - repmat(Mu.row(k), Rf(r).n_rows, 1)) % trans(repmat(Rf(r).col(k), 1, q))) * 
          (Ez(r).slice(k) - repmat(Mu.row(k), Rf(r).n_rows, 1)) + Nmat(r,k)*  Varz(r).slice(k); //yangyi n-by-n matrix to n-by-q matrix
        
      }
      if(Sigma_diag){
        Sigma0.slice(k) = diagmat(Smat) / accu(Nmat.col(k));
      }else{
        Sigma0.slice(k) = Smat / accu(Nmat.col(k));
      }
      // Sigma0.slice(k).print();
    }
  }   
  
  return Sigma0;
}

// function2
cube update_embed_Sigma2(const field<mat>& Rf, const field<cube>& Ez, const field<cube>& Varz,
                         const mat& Mu, const mat& Nmat,
                         const bool& homoClust=false, const bool& Sigma_diag=false){
  
  // homoCluste implies whether Sigma_k = Sigma_l, for all k,l.
  int i, r,k,n, K= Mu.n_rows, q=Mu.n_cols, r_max=Rf.n_elem;
  cube Sigma0(q,q , K);
  mat Smat(q, q, fill::zeros);
  // Rprintf("r_max = %d, K= %d, q= %d", r_max, K, q);
  for(k = 0; k<K; ++k){  
    Smat = zeros<mat>(q,q);
    for(r=0; r< r_max; ++r){
      n = Rf(r).n_rows;
      for(i=0; i<n; ++i){
        Smat += Rf(r)(i,k) * (trans((Ez(r).slice(k).row(i) - Mu.row(k)))*
          (Ez(r).slice(k).row(i) - Mu.row(k)) + Varz(r).slice(k)); 
      }
      
      
    }
    if(Sigma_diag){
      Sigma0.slice(k) = diagmat(Smat) / accu(Nmat.col(k));
    }else{
      Sigma0.slice(k) = Smat / accu(Nmat.col(k));
    }
    //Sigma0.slice(k).print();
  }
  
  
  return Sigma0;
}





// update Lambda0
mat update_embed_Lam1(const field<mat>& Vf, const field<mat>& Rf, const mat& Nmat,
                      const field<cube>& Ezv, const field<cube>& Varzv, const bool& homo = false){
  
  int r, k, K= Rf(0).n_cols, q = Vf(0).n_cols, r_max=Rf.n_elem;
  vec Lsum(q,fill::zeros);
  mat Lam(r_max, q);
  for(r=0; r< r_max; ++r){
    Lsum = zeros<vec>(q,1);
    mat tmpXk;
    for(k=0; k<K; ++k){
      tmpXk = Vf(r) - Ezv(r).slice(k);
      Lsum += trans(sum(tmpXk % tmpXk % repmat(Rf(r).col(k), 1, q))) + Nmat(r,k) * diagvec(Varzv(r).slice(k));
    }
    if(homo){
      Lam.row(r) = mean(Lsum)* ones<mat>(1, q) / (Vf(r).n_rows*1.0);
    }else{
      Lam.row(r) = Lsum.t()/(Vf(r).n_rows*1.0);
    }
    
  }
  // replace little values to increase the stability of the algorithm.
  uvec id1 =  find(Lam <1e-7);
  Lam(id1) = 1e-7*ones<vec>(id1.n_elem, 1);
  return Lam; 
}

// update Psi
cube update_embed_Psi(const field<mat>& Muf, const field<mat>& Rf, const field<cube>& Ef,
                      const field<cube>& Varf, const mat& Nmat){
  
  int r, k, n, M= Muf.n_elem, q = Muf(0).n_cols, K= Rf(0).n_cols;
  cube Psi0(q,q, M, fill::zeros);
  mat Psi_new;
  for(r=0; r < M; ++r){
    
    Psi_new = zeros(q,q);
    // update Psi0
    
    for(k=0; k< K; ++k){
      Psi_new += trans(Muf(r)- Ef(r).slice(k)) * ((Muf(r)- Ef(r).slice(k)) % repmat(Rf(r).col(k), 1, q)) +
        Nmat(r,k) * Varf(r).slice(k);
    }
    n = Rf(r).n_rows;
    Psi0.slice(r) = Psi_new/ n;
  }
  return Psi0;
}

//update beta
vec update_embed_beta(const arma::field<ivec>& yf, const field<mat>& Rf, const field<sp_mat>& Adjf,
                      const vec& beta_grid, const vec& alpha0, const bool& mix_prop_heter=true){
  
  int r, k,  M = Rf.n_elem, K = Rf(0).n_cols;
  // update beta: grid search.
  int ng_beta = beta_grid.n_elem;
  vec objBetaVec(ng_beta), beta0(M);
  if(mix_prop_heter){
    for(r=0; r < M; ++r){
      objBetaVec = zeros(ng_beta,1);
      for(k = 0; k < ng_beta; ++k){ // each sample has a smoothing parameter
        objBetaVec(k) =  objr_beta_embed(yf(r), Rf(r), Adjf(r), K, alpha0, beta_grid(k)); 
      }
      beta0(r) = beta_grid(index_max(objBetaVec));
    }
  }else{
    for(k=0; k < ng_beta; ++k){ // all samples have a same smoothing parameter.
      objBetaVec(k) = obj_beta_embed(yf, Rf, Adjf, K, alpha0, beta_grid(k));
    }
    beta0 = ones(M, 1) * beta_grid(index_max(objBetaVec));
  }
  
  return beta0;
}

//
arma::mat get_Vmean_embed(const arma::mat& V, const arma::sp_mat& Adj){
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



// //
// double Q_fun(const field<mat>& Vf, const field<mat>& Rf,
//              const field<cube>& Ez, const arma::field<cube>& Varz, 
//              const field<cube>& Ev, const arma::field<cube>& Varv, 
//              const field<cube>& Ezv, const arma::field<cube>& Varzv, 
//              const field<mat>& Muv,
//              const mat& Mu0, const vec& alpha0,const vec& beta0, const cube& Sigma0, 
//              const mat& Lam0, const cube& Psi0, const bool& Sp_embed){
//   
//   double Q = 0, tmp_scalar1 =0, tmp_scalar2 = 0, tmp_scalar3= 0, tmp_scalar4=0;
//   int i, k, r, K= Rf(0).n_cols,  M=Vf.n_elem,  q = Ez(0).n_cols;
//   // bool Temporal_embed = false;
//   
//   //Rprintf("good varzv: \n");
//   // P(X)
//   for(r=0; r< M; ++r){
//     //Rprintf("good Px: \n");
//     for(k=0; k<K; ++k){
//       for(i = 0; i < Vf(r).n_rows; i++){
//         tmp_scalar1 +=  Rf(r)(i,k) *( 0.5* accu(log(Lam0.row(r))) +
//           0.5* accu((Vf(r).row(i) - Ezv(r).slice(k).row(i))%
//           (Vf(r).row(i) -  Ezv(r).slice(k).row(i))/ Lam0.row(r)) + 
//           0.5* trace(diagmat(1.0/Lam0.row(r))*Varzv(r).slice(k))) + 0.5*q* log(2* datum::pi);
//         
//       }
//     }
//     
//   }
//   
//   // P(Z)
//   int nm;
//   for(r=0; r< M; ++r){
//     //Rprintf("good Pz: \n");
//     nm = Vf(r).n_rows;
//     for(k=0; k<K; ++k){
//       for(i = 0; i < nm; i++){
//         //mat tmpMat = ;
//         tmp_scalar2 += Rf(r)(i,k) *(0.5* log(det(Sigma0.slice(k)))+ 0.5* trace(Sigma0.slice(k).i() *
//           Varz(r).slice(k))+ 0.5* as_scalar((Ez(r).slice(k).row(i) - Mu0.row(k))*
//           Sigma0.slice(k).i()*(Ez(r).slice(k).row(i) - Mu0.row(k)).t() )) + 0.5*q* log(2* datum::pi);
//       }
//     }
//   }
//   
//   
//   // P(V)
//   if(Sp_embed){
//     for(r=0; r< M; ++r){
//       //Rprintf("good Pv: \n");
//       for(k=0; k<K; ++k){
//         for(i = 0; i < Vf(r).n_rows; i++){
//           tmp_scalar3 +=  Rf(r)(i,k) *(0.5* log(det(Psi0.slice(r))) + 0.5* as_scalar((Ev(r).slice(k).row(i) - Muv(r).row(i))*
//             Psi0.slice(r).i()*trans(Ev(r).slice(k).row(i) - Muv(r).row(i)) ) +
//             0.5*trace(Psi0.slice(r).i()* Varv(r).slice(k)));
//           
//         }
//       }
//       
//     }
//     
//   }
//   
//   
//   
//   Q = - (tmp_scalar1 + tmp_scalar2 + tmp_scalar3);
//   
//   
//   return Q;
// }




// Run ICM-step and partial E-step in a function for efficient computation.
void runICM_sp_embed(const arma::field<mat>& Vf, arma::field<ivec>& yf, arma::field<mat>& Ff, 
                     const arma::mat& Lam0,  const arma::mat& Mu0,
                     const arma::cube& Sigma0, cube& Psi0, 
                     const arma::field<sp_mat>& Adjf, const arma::vec& alpha0,
                     const vec& beta0, int maxIter_ICM,
                     double& loglik,  field<mat>& Muf, field<mat>& Rf, 
                     field<cube>& Ez, field<cube>& Varz,  field<cube>& Ef, field<cube>& Varf, 
                     field<cube>& Ezf, field<cube>& Varzf,
                     const bool& Sp_embed){
  
  
  // basic info.
  int r, r_max = Vf.n_elem, K = Mu0.n_rows;
  int i, iter, k, n;
  
  // two cached objects used for parameters update.
  field<mat> Ux(r_max);
  double  logdSk;
  // Efaluate energy of x, Ux
  arma::mat SrkI, VSrk;
  
  for(r = 0; r < r_max; ++r){
    // Efaluate energy of x, Ux
    n = Vf(r).n_rows; // tmp object
    Ux(r) = zeros(n, K);
    //Rprintf("Efaluate posterior! \n");
    vec mSk;
    for (k = 0; k < K; k++)	{
      //Rprintf("k= %d \n", k);
      SrkI = inv_sympd(Psi0.slice(r) + diagmat(Lam0.row(r)) + Sigma0.slice(k));
      
      
      VSrk = (Vf(r) -Muf(r)-repmat(Mu0.row(k), n, 1))* SrkI; // O(nq^2)
      
      
      Ez(r).slice(k) = VSrk*Sigma0.slice(k) + repmat(Mu0.row(k), n, 1); // complexity: O(nq^2)
      Varz(r).slice(k) = Sigma0.slice(k) - Sigma0.slice(k)*SrkI*Sigma0.slice(k); // O(q^3)
      
      // save some objects for posterior of v
      Ef(r).slice(k) = VSrk *Psi0.slice(r) + Muf(r); // 
      Varf(r).slice(k) = Psi0.slice(r) - Psi0.slice(r)* SrkI *Psi0.slice(r);
      
      Ezf(r) = Ez(r) + Ef(r); // use the addition of two cubes.
      Varzf(r).slice(k) = Varz(r).slice(k)+ Varf(r).slice(k) - Psi0.slice(r)*SrkI*Sigma0.slice(k) -
        Sigma0.slice(k)*SrkI*Psi0.slice(r);
      
      // Efaluate energy of x, Ux for updating y and caculating responsibility R
      multi_det_Sk_embedCpp2(Vf(r) - Muf(r), SrkI, Mu0.row(k), // Use SVD to speed up.
                             logdSk, mSk);
      
      Ux(r).col(k) = -0.5*logdSk  + 0.5 * mSk; // calculate energy by column.
      
    }
  }
  
  // Estimate Y by ICM
  arma::vec Energy(maxIter_ICM);
  //--------------------------------------------------------------------------------	
  //--------- ICM algrithm to predict V------------------------------------
  field<mat> U(r_max);
  // Rprintf("predict Y and V! \n");
  for(r= 0; r< r_max; r++){
    
    Energy(0) = INFINITY;
    arma::mat Uy1, U1;
    vec U1min;
    uvec y1_u;
    n = Vf(r).n_rows; // tmp object
    //Rprintf("predict Y and V, r = %d ! \n", r);
    // ivec y1;
    for (iter = 1; iter < maxIter_ICM;  ++iter ) {
      //Rprintf("predict Y and V,  iter = %d ! \n", iter);
      Uy1 = calYenergy2D_sp_embed(yf(r), Adjf(r), K, alpha0, beta0(r));
      //Rprintf("predict Y and V,  iter = %d ! \n", iter);
      U1 = Uy1 + Ux(r); // log likelihood of (x, y).
      U1min = min(U1, 1);
      y1_u = index_min(U1, 1);
      yf(r) = conv_to< ivec >::from(y1_u) + 1;
      // Rprintf("predict Y and V,  iter = %d ! \n", iter);
      if(Sp_embed){
        Muf(r) = get_Vmean_embed(Ff(r),  Adjf(r));
        
        // for(i=0; i<n; i++){  // O(nq^2*maxIter_ICM)
        //   k = y1_u(i);
        //   Ff(r).row(i)= ((Vf(r).row(i)- Mu0.row(k))* (Sigma0.slice(k)+ diagmat(Lam0.row(r))).i() + Muf(r).row(i)* Psi0.slice(r).i())*
        //        inv_sympd((Sigma0.slice(k)+ diagmat(Lam0.row(r))).i() + Psi0.slice(r).i());
        // 
        // }
        
        for(k = 0; k<K; ++k){
          uvec index_k = find(y1_u == k);
          int nk = index_k.n_elem;
          // Rprintf("k= %d,  nk = %d ! \n", k, nk);
          if(nk > 0){// if the number of spots whose cluster is k is greater than 0
            Ff(r).rows(index_k) = ((Vf(r).rows(index_k)- repmat(Mu0.row(k), nk,1))* (Sigma0.slice(k)+ diagmat(Lam0.row(r))).i()  + Muf(r).rows(index_k)* Psi0.slice(r).i()) * 
              inv_sympd((Sigma0.slice(k)+ diagmat(Lam0.row(r))).i() + Psi0.slice(r).i());
            
          }
          
        }
        
        // centering V1 for identifiability.
        // if(r == 0) Vf(r) =   Vf(r)- repmat(mean(Vf(r)), Vf(r).n_rows, 1);
        // Since energy_V is computationally high, we do not caculate it.
        // Energy(iter) = energy_V(X, V, W0, Lam_vec0, Muv, Mu0, Sigma0,Psi0,y, Cki) + sum(Umin); // 
        
      }
      
      // Rprintf("predict Y and V! \n");
      Energy(iter) = sum(U1min);
      if (Energy(iter) - Energy(iter - 1) > 1e-5) {
        // cout << "diff Energy = " << Energy(iter) - Energy(iter - 1)  << endl;
        // Rprintf("diff Energy = %4f \n", Energy(iter) - Energy(iter - 1));
        break;
      }
      
      if(!Sp_embed){
        if (Energy(iter-1) - Energy(iter) < 1e-5)
        {
          Rprintf("ICM Converged at Iteration = %d \n", iter);
          break;
        }
      }
    }
    U(r) = U1;
    
  }
  
  // calculate R and pseudo observed loglikelihood
  loglik=0;
  for(r=0; r< r_max; ++r){
    vec maxA1 = max(-U(r), 1);
    //cout<<"good3"<<endl;
    mat negU = (-U(r) - repmat(maxA1, 1, K));
    vec loglik_more_vec = sum(exp(negU),1);
    loglik += sum(log(loglik_more_vec) + maxA1);  
    Rf(r) = exp(negU) / repmat(loglik_more_vec, 1, K);
  }
  
}  


Obj_SCMEBTwo SepSpatialClusterCpp_GivenK(
    const arma::field<arma::mat>& V_init, 
    const arma::field<arma::sp_mat>& Adjf, 
    const arma::ivec& y_init,
    const arma::mat& Mu_int, 
    const arma::cube& Sigma_int,
    const arma::cube& Psi_int, 
    const vec& beta_init, 
    const arma::vec& beta_grid,
    const int& maxIter_ICM,
    const int& maxIter, 
    const double& epsLogLik, 
    const bool& verbose,
    const bool& homo, 
    const bool&  Sigma_equal,
    const bool& Sigma_diag, 
    const bool& Sp_embed){
  // homo denotes error covariance is a scalar * identity matrix.
  // basic info
  int r, M = V_init.n_elem; // get the number of data source
  int n, q = Mu_int.n_cols, K = Mu_int.n_rows;
  bool mix_prop_heter = true;
  
  vec alpha0(K, fill::zeros);
  vec beta0(beta_init);
  arma::field<arma::mat> Vf(M);
  arma::field<arma::ivec> yf(M);
  int n0 = 0;
  int n1 = -1;
  for(r=0; r < M; ++r){ 
    Vf(r) = V_init(r);
    n1 = n1 + Vf(r).n_rows;
    yf(r) = y_init.subvec(n0, n1);
    n0 = n0 + Vf(r).n_rows;
  }
  
  // Initialize the model parameters
  mat Mu0(Mu_int);
  cube Sigma0(Sigma_int), Psi0(Psi_int);
  if(!Sp_embed) Psi0 = zeros<cube>(q,q, M);
  
  mat Lam0 = zeros(M,q);
  
  vec loglik(maxIter);
  loglik(0) = INT_MIN;
  // vec Qvec(loglik);
  
  // pseudo obserbed loglikelihood.
  double loglikVal;
  int iter;
  
  if (verbose) Rprintf("Finish variable initialization \n");
  // create posterior expectation and covariance of z+v+u.
  field<mat> Rf(M),  Ff(M), Muf(M);
  // Posterior of z, v and z+v
  field<cube> Ezf(M), Varzf(M), Ez(M), Varz(M), Ef(M),  Varf(M);
  for(r = 0; r < M; ++r){// initialize the shape of cube in these fields.
    n = Vf(r).n_rows; // tmp object
    Muf(r) = zeros(n, q);
    Ff(r) = Muf(r);
    // Rprintf("nrow(Muv(0)) = %d, nrow(Vf(0))=%d\n", Muv(r).n_rows, Vf(r).n_rows);
    Ezf(r) = zeros(n, q, K);
    Ez(r) = Ezf(r);
    Ef(r) = Ezf(r);
    Varzf(r) = zeros(q,q,K);
    Varz(r) = Varzf(r);
    Varf(r) = Varzf(r);
  }
  
  mat Nmat(M, K);
  // begin ICM-EM algorithm
  for(iter = 1; iter < maxIter; iter++){
    
    // predict V and y.
    //Rprintf("Satrt ICM step! \n");
    runICM_sp_embed(Vf, yf, Ff, Lam0, Mu0, Sigma0, Psi0, 
                    Adjf, alpha0, beta0, maxIter_ICM, 
                    loglikVal,  Muf, Rf, Ez, Varz, Ef, Varf, Ezf, Varzf, Sp_embed);
    loglik(iter) = loglikVal; // obtain the pseudo observed log-likelihood.
    //Rprintf("Finish ICM and E-step! \n");
    
    
    //double Q1 = Q_fun(Xf, tvec, Ez, Varz, Ef, Varf, Ezf, Varzf, Muv, W0, Mu0, beta0, Sigma0, Lam0, Psi0, Sp_embed);
    //double elbo = ELBOcpp(Xf, tvec, W0, alpha0, beta0, Sigma0, Lam0, Psi0, Sp_embed);
    // loglik(iter) = Q1;
    //Qvec(iter) = Q1;
    
    // compute N
    for(r=0; r< M; ++r){
      Nmat.row(r) = sum(Rf(r));
    }
    
    // double Q2 = Q_fun(Vf,Rf, Ez, Varz, Ef, Varf, Ezf, Varzf, Muf, Mu0,alpha0, beta0, Sigma0, Lam0, Psi0,Sp_embed);
    
    
    // update Mu0
    Mu0 = update_embed_Mu1(Rf, Ez, Nmat);
    
    // double Q3= Q_fun(Vf,Rf, Ez, Varz, Ef, Varf, Ezf, Varzf, Muf, Mu0,alpha0, beta0, Sigma0, Lam0, Psi0,Sp_embed);
    // Rprintf("dQ_Mu= %4f \n", Q3-Q2);
    // update Sigma0
    //Rprintf("Run M-step Sigma0! \n");
    
    Sigma0 = update_embed_Sigma1(Rf, Ez, Varz, Mu0, Nmat, Sigma_equal, Sigma_diag); // 
    // double Q4 = Q_fun(Vf,Rf, Ez, Varz, Ef, Varf, Ezf, Varzf, Muf, Mu0, alpha0, beta0, Sigma0, Lam0, Psi0,Sp_embed);
    // Rprintf("dQ_Sigma= %4f \n", Q4-Q3);
    // double Q5 =  Q_fun(Vf,Rf, Ez, Varz, Ef, Varf, Ezf, Varzf, Muf, Mu0,alpha0, beta0, Sigma0, Lam0, Psi0,Sp_embed);
    // Rprintf("dQ_Lam= %4f \n", Q5-Q4);
    
    // update Psi, the covariance of v
    //Rprintf("Run VB-M-step Psi! \n");
    Psi0 = update_embed_Psi(Muf, Rf, Ef, Varf, Nmat);
    // double Q6 = Q_fun(Vf,Rf, Ez, Varz, Ef, Varf, Ezf, Varzf, Muf, Mu0,alpha0, beta0, Sigma0, Lam0, Psi0,Sp_embed);
    // Rprintf("dQ_Psi= %4f \n", Q6-Q5);
    
    // update beta
    beta0 = update_embed_beta(yf, Rf, Adjf, beta_grid, alpha0, mix_prop_heter);
    // double Q7 = Q_fun(Vf,Rf, Ez, Varz, Ef, Varf, Ezf, Varzf, Muf, Mu0,alpha0, beta0, Sigma0, Lam0, Psi0,Sp_embed);
    // Rprintf("dQ_beta= %4f \n", Q7-Q6);
    
    // calculate loglikelihood
    if(loglik(iter)  - loglik(iter-1)   < -1e-7){
      // perror("The likelihood failed to increase!");
      //break;
    }
    
    //finish1 = clock();
    //cout << finish1 - start1 << "/" << CLOCKS_PER_SEC << " (s) " << endl;
    // output algorithm info.
    if(verbose){
      Rprintf("K = %d, iter = %d, loglik= %4f, dloglik=%4f \n", 
              K, iter +1, loglik(iter), abs(loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
    }
    if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < epsLogLik) break;
    // if(abs(Q  - tmp_Q) < epsLogLik) break;
    
  }
  
  field<mat> Ezz(M);
  for(r=0; r< M; ++r){
    Ezz(r) = zeros<mat>(Rf(r).n_rows, q);
    for(int k=0; k<K; k++){
      
      Ezz(r) +=   Ez(r).slice(k) % repmat(Rf(r).col(k), 1, q);
    }
  }
  
  Obj_SCMEBTwo output;
  output.cluster    = yf;
  output.hZ         = Ezz;
  output.hV         = Vf;
  output.Rf         = Rf;
  output.Mu         = Mu0;
  output.Sigma      = Sigma0;
  output.Psi        = Psi0;
  output.beta       = beta0;
  output.loglik     = loglik(iter-1);
  output.dLogLik    = loglik(iter-1)  - loglik(iter-2);
  output.loglik_seq = loglik.subvec(0, iter-1);
  
  return(output);
}


// [[Rcpp::export]]
Rcpp::List iSCMEBCpp(
    const Rcpp::List& vList, 
    const Rcpp::List& Adjlist, 
    const Rcpp::List& yList_int,
    const Rcpp::List& Mu_int, 
    const Rcpp::List& Sigma_int,
    const arma::cube& Psi_int, 
    const double& beta_int, 
    const arma::vec& beta_grid,
    const int& maxIter_ICM,
    const int& maxIter, 
    const double& epsLogLik, 
    const bool& verbose,
    const bool& homo, 
    const bool&  Sigma_equal,
    const bool& Sigma_diag, 
    const bool& Sp_embed, 
    const arma::uword& maxK,
    const arma::uword& minK,
    const int& coreNum){
  // homo denotes error covariance is a scalar * identity matrix.
  int r, M = vList.length(); // get the number of data source
  field<mat> Vf(M);
  field<sp_mat> Adjf(M);
  vec  beta0(M);
  // transfer list to field
  for(r=0; r < M; ++r){ 
    mat Xtmp = vList[r]; // enforce to become a matrix.
    Vf(r) = Xtmp;
    sp_mat Adjtmp = Adjlist[r];
    Adjf(r) = Adjtmp;
    
    beta0(r) = beta_int;
  }
  
  // transfer list
  uword lengthK = maxK - minK + 1;
  arma::field<arma::mat> Mu0(lengthK);
  arma::field<arma::cube> Sigma0(lengthK);
  arma::field<arma::ivec> yf(lengthK);
  
  for (uword k = 0; k < lengthK; k++) {
    mat Mutmp = Mu_int[k]; // enforce to become a matrix.
    Mu0(k) = Mutmp;
    
    cube Sigmatmp = Sigma_int[k];
    Sigma0(k) = Sigmatmp;
    
    ivec ytmp = yList_int[k];
    yf(k) = ytmp;
  }
  
  if(lengthK==1){
    
    struct Obj_SCMEBTwo output[1];
    output[0] = SepSpatialClusterCpp_GivenK(
      Vf, Adjf, yf(0), Mu0(0), Sigma0(0), Psi_int, beta0, beta_grid, maxIter_ICM, 
      maxIter, epsLogLik, verbose, homo, Sigma_equal, Sigma_diag, Sp_embed); 
    
    
    List Obj_SCMEBTwo_Rcpp(maxK-minK+1);
    
    for (uword k = 0; k < maxK - minK + 1; k++){
      // output return value
      Obj_SCMEBTwo_Rcpp[k] = List::create(
        Rcpp::Named("cluster")    = output[k].cluster,
        Rcpp::Named("hZ")         = output[k].hZ,
        Rcpp::Named("hV")         = output[k].hV,
        Rcpp::Named("Rf")         = output[k].Rf,
        Rcpp::Named("Mu")         = output[k].Mu,
        Rcpp::Named("Sigma")      = output[k].Sigma,
        Rcpp::Named("Psi")        = output[k].Psi,
        Rcpp::Named("beta")       = output[k].beta,
        Rcpp::Named("loglik")     = output[k].loglik,
        Rcpp::Named("dLogLik")    = output[k].dLogLik,
        Rcpp::Named("loglik_seq") = output[k].loglik_seq);
    }
    return(Obj_SCMEBTwo_Rcpp);
  }
  
  //set parallel structure object
  par_SCMEBTwo parObj(
      Vf, Adjf, yf, Mu0, Sigma0, Psi_int, beta0, beta_grid, maxIter_ICM, maxIter, 
      epsLogLik, verbose, homo, Sigma_equal, Sigma_diag, Sp_embed, maxK, minK);
  
  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  
  for (int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&par_SCMEBTwo::update_by_thread_SCMEBTwo, &parObj, i_thread);
  }
  for (int i = 0; i < n_thread; i++){
    threads[i].join();
  }
  
  List Obj_SCMEBTwo_Rcpp(maxK-minK+1);
  
  for (uword k = 0; k < maxK - minK + 1; k++){
    // output return value
    Obj_SCMEBTwo_Rcpp[k] = List::create(
      Rcpp::Named("cluster")    = parObj.output[k].cluster,
      Rcpp::Named("hZ")         = parObj.output[k].hZ,
      Rcpp::Named("hV")         = parObj.output[k].hV,
      Rcpp::Named("Rf")         = parObj.output[k].Rf,
      Rcpp::Named("Mu")         = parObj.output[k].Mu,
      Rcpp::Named("Sigma")      = parObj.output[k].Sigma,
      Rcpp::Named("Psi")        = parObj.output[k].Psi,
      Rcpp::Named("beta")       = parObj.output[k].beta,
      Rcpp::Named("loglik")     = parObj.output[k].loglik,
      Rcpp::Named("dLogLik")    = parObj.output[k].dLogLik,
      Rcpp::Named("loglik_seq") = parObj.output[k].loglik_seq);
  }
  return(Obj_SCMEBTwo_Rcpp);
}


// =============================================================================
// mt_parallel_job.cpp
// =============================================================================
void par_SCMEBTwo::loop_by_K_SCMEBTwo(int g){ 
  arma::ivec y0     = y_init(g);
  arma::mat  Mu1    = Mu0(g);
  arma::cube Sigma1 = Sigma0(g);
  
  output[g] = SepSpatialClusterCpp_GivenK(
    V_init, Adjf, y0, Mu1, Sigma1, Pi0, beta0, beta_grid, maxIter_ICM, 
    maxIter, epsLogLik, verbose, homo, Sigma_equal, Sigma_diag, Sp_embed); 
  
  // reset to free memory
  y0.reset();
  Mu1.reset();
  Sigma1.reset();
}

std::mutex _mtx22;
int par_SCMEBTwo::next_SCMEBTwo() {
  std::lock_guard<std::mutex> lockGuard(_mtx22);
  if (current_idx >= maxK - minK + 1) {
    return -1;
  }
  current_idx++;
  return current_idx - 1;
}

void par_SCMEBTwo::update_by_thread_SCMEBTwo(int thread_id) {
  while (true){
    int idx = next_SCMEBTwo();
    if (idx == -1) {
      break;
    }
    loop_by_K_SCMEBTwo(idx);
  }
}


// =============================================================================
// getNB_fast.cpp
// =============================================================================
// [[Rcpp::export]]
arma::sp_umat getneighborhood_fast(const arma::mat x, double radius)	{
  int N = x.n_rows;
  arma::sp_umat D(N, N);
  double dis;
  uvec idx, idx2;
  for (int j = 0; j < N-1; ++j)
  {    
    idx = find(abs(x(j,0) - x.col(0))<radius); 
    idx2 = find(idx>j);
    int p = idx2.n_elem;
    for (int i = 0; i < p; ++i)
    {
      dis = norm(x.row(idx(idx2(i))) - x.row(j), 2);
      if (dis < radius){
        D(idx(idx2(i)),j) = 1;
        D(j,idx(idx2(i))) = 1;
      }
    }
  }
  return D;
}

// =============================================================================
// wpca.cpp
// =============================================================================
//[[Rcpp::export]]   
Rcpp::List wpcaCpp(const arma::mat&X, const int& nPCs, const bool& weighted=true) {
  arma::mat U, V;
  arma::vec s;
  arma::mat PCs, loadings;
  svd_econ(U, s, V, X);
  PCs = U.cols(0, nPCs-1) *diagmat(s.subvec(0, nPCs-1));
  loadings = V.cols(0, nPCs-1);
  arma::mat dX = PCs * loadings.t() - X;
  arma::rowvec Lam_vec = mean(dX % dX);
  if(weighted){
    svd_econ(U, s, V, X*diagmat(1.0/ sqrt(Lam_vec)));
    // vec s2 =  s % s; // s; //
    arma::mat loadings_unscale = diagmat(sqrt(Lam_vec)) * V.cols(0, nPCs-1);
    arma::mat  V1;
    arma::vec s1;
    svd_econ(loadings, s1, V1, loadings_unscale);
    PCs = U.cols(0, nPCs-1) * diagmat(s.subvec(0, nPCs-1)) * V1 * diagmat(s1);
    dX = PCs * loadings.t() - X;
    Lam_vec = mean(dX % dX);
  }
  List output = List::create(
    Rcpp::Named("PCs") = PCs,
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("Lam_vec") = Lam_vec);
  
  return output;
}
