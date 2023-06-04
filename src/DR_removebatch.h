#ifndef SCMEBTwo_h
#define SCMEBTwo_h

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

struct Obj_SCMEBTwo{
    arma::field<ivec> cluster;
    arma::field<arma::mat> hZ;
    arma::field<arma::mat> hV;
    arma::field<arma::mat> Rf;
    arma::mat Mu;
    arma::cube Sigma;
    arma::cube Psi;
    vec beta;
    double loglik;
    double dLogLik;
    vec loglik_seq;
};

Obj_SCMEBTwo SepSpatialClusterCpp_GivenK(
  const arma::field<arma::mat>& V_init, 
  const arma::field<arma::sp_mat>& Adjf, 
  const arma::field<arma::ivec>& y_init,
  const arma::mat& Mu_int, 
  const arma::cube& Sigma_int,
  const arma::cube& Psi_int, 
  const vec& beta_init, 
  const arma::vec& beta_grid,
  const int& maxIter_ICM,
  const int& maxIter, 
  const double& epsLogLik, 
  const bool& verbose,
  const bool& homo = false, 
  const bool&  Sigma_equal=false,
  const bool& Sigma_diag=false, 
  const bool& Sp_embed=false
);

class par_SCMEBTwo{
public:
	arma::field<arma::mat> V_init;
    arma::field<arma::sp_mat> Adjf;
    arma::field<arma::ivec> y_init;
    arma::field<arma::mat> Mu0;
    arma::field<arma::cube> Sigma0;
    arma::cube Pi0;
    vec beta0;
    vec beta_grid;
    int maxIter_ICM;
    int maxIter;
    double epsLogLik;
    bool verbose;
    bool homo;
    bool Sigma_equal;
    bool Sigma_diag;
    bool Sp_embed;
    arma::uword maxK;
    arma::uword minK;
    int g;    
    uword current_idx = 0;
    struct Obj_SCMEBTwo output[50];

    par_SCMEBTwo(
        const arma::field<arma::mat>& V_init, 
        const arma::field<arma::sp_mat>& Adjf, 
        const arma::field<arma::ivec>& y_init,
        const field<arma::mat>& Mu0, 
        const field<arma::cube>& Sigma0,
        const arma::cube& Pi0, 
        const vec& beta_init, 
        const arma::vec& beta_grid,
        const int& maxIter_ICM,
        const int& maxIter, 
        const double& epsLogLik, 
        const bool& verbose,
        const bool& homo, 
        const bool& Sigma_equal,
        const bool& Sigma_diag, 
        const bool& Sp_embed,
        const arma::uword& maxK,
        const arma::uword& minK
    ) {
        this->V_init      = V_init;
        this->Adjf        = Adjf;
        this->y_init      = y_init;
        this->Mu0         = Mu0;
        this->Sigma0      = Sigma0;
        this->Pi0         = Pi0;
        this->beta0       = beta_init;
        this->beta_grid   = beta_grid;
        this->maxIter_ICM = maxIter_ICM;
        this->maxIter     = maxIter;
        this->epsLogLik   = epsLogLik;
        this->verbose     = verbose;
        this->homo        = homo;
        this->Sigma_equal = Sigma_equal;
        this->Sigma_diag  = Sigma_diag;
        this->Sp_embed    = Sp_embed;
        this->maxK = maxK;
        this->minK = minK;
    }

	void loop_by_K_SCMEBTwo(int g);
	void update_by_thread_SCMEBTwo(int thread_id);
	int next_SCMEBTwo();
};

#endif
