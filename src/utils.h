#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <random>
#include <chrono>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

// Aux functions
void reportProgress(int iter, int nBurnin, int nSims, int nReport, int save_idx = 0,
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now());
arma::mat dist1(const arma::mat& coords);
arma::mat dist2(const arma::mat& coords);
arma::mat dist_mat(const arma::mat& A, const arma::mat& B);
double dtnorm(double x, double mu, double sigma, double a, double b);
double rtnorm(double mu, double sigma, double a, double b);
double psi1(double x, double alpha, double lambda);
double psi2(double x, double alpha, double lambda);
arma::vec rgig(const int N, const double a, const arma::vec b, const double nu);
arma::vec rig(const int N, const arma::vec& mu, const double lambda, 
  const bool parallel = false, const int nThreads = 0);
arma::vec ralRcpp(const arma::vec sigma, const double tau);
arma::vec RandomMultiNormalC(const arma::mat& Q, const arma::vec& b);
arma::mat ffbs(const arma::mat& Y, const double rho, const arma::mat& Sigma_w, 
               const std::vector<arma::vec>& D_list);
arma::mat dotW(const arma::mat& R, const arma::vec& Wtls,
               const std::vector<arma::uvec>& indW, const double rho, 
               const int T, const arma::ivec& L, const double onemrho2);
    
#endif