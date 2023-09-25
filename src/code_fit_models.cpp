#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double dtnorm(double x, double mu, double sigma, double a, double b) { // log
  return R::dnorm(x, mu, sigma, 1) - std::log(R::pnorm(b, mu, sigma, 1, 0) - R::pnorm(a, mu, sigma, 1, 0));
}

// [[Rcpp::export]]
double rtnorm(double mu, double sigma, double a, double b) {
  double alpha = (a - mu) / sigma;
  double beta  = (b - mu) / sigma;
  double pAlpha = R::pnorm(alpha, 0, 1, 1, 0);
  double pBeta  = R::pnorm(beta, 0, 1, 1, 0);
  return R::qnorm(pAlpha + R::runif(0, 1) * (pBeta - pAlpha), 0, 1, 1, 0) * sigma + mu;
}

// [[Rcpp::export]]
double psi1(double x, double alpha, double lambda){
  Rcpp::NumericVector coshx = { x };
  coshx = Rcpp::cosh(coshx);
  return - alpha * (coshx(0) - 1) - lambda * (exp(x) - x - 1);
}
// [[Rcpp::export]]
double psi2(double x, double alpha, double lambda){
  Rcpp::NumericVector sinhx = { x };
  sinhx = Rcpp::sinh(sinhx);
  return - alpha * sinhx(0) - lambda * (std::exp(x) - 1);
}
// [[Rcpp::export]]
arma::vec rgig(
    const int N, 
    const double a, 
    const arma::vec b, 
    const double nu) {
  
  arma::vec X(N);
  static const double cosh1 = 1.543081;
  
  double lambda = nu;
  arma::vec omega = arma::sqrt(a * b);
  
  double alpha;
  double t, s;
  double tp, sp, q;
  double eta, zeta, theta, xi;
  double p, r;
  double chi;
  double U, V, W;
  
  double aux;
  
  Rcpp::NumericVector v(2);
  
  for (int i = 0; i < N; ++i) {
    alpha = std::sqrt(omega(i) * omega(i) + lambda * lambda) - lambda;
    aux = - psi1(1, alpha, lambda);
    if (aux > 2) {
      t = std::sqrt(2 / (alpha + lambda));
    } else if (aux < 0.5) {
      t = std::log(4 / (alpha + 2 * lambda));
    } else {
      t = 1;
    }
    aux = - psi1(-1, alpha, lambda);
    if (aux > 2) {
      s = std::sqrt(4 / (alpha * cosh1 + lambda));
    } else if (aux < 0.5) {
      aux = 1 / alpha;
      v = { 1 / lambda, log(1 + aux + std::sqrt(aux * aux + 2 * aux)) };
      s = Rcpp::min(v);
    } else {
      s = 1;
    }
    
    eta   = - psi1( t, alpha, lambda);
    zeta  = - psi2( t, alpha, lambda);
    theta = - psi1(-s, alpha, lambda);
    xi    =   psi2(-s, alpha, lambda);
    
    p = 1 / xi;
    r = 1 / zeta;
    
    tp = t - r * eta;
    sp = s - p * theta;
    q  = tp + sp;
    
    do {
      U = R::runif(0, 1);
      V = R::runif(0, 1);
      W = R::runif(0, 1);
      if (U < (q / (p + q + r))) {
        X(i) = - sp + q * V;
      } else if (U < ((q + r) / (p + q + r))) {
        X(i) = tp + r * std::log(1 / V);
      } else {
        X(i) = - sp - p * std::log(1 / V);
      }
      if (X(i) > tp) {
        chi = exp(- eta - zeta * (X(i) - t));
      } else if (X(i) < (- sp)) {
        chi = exp(- theta + xi * (X(i) + s));
      } else {
        chi = 1;
      }
    } while ((W * chi) > std::exp(psi1(X(i), alpha, lambda)));
      
    aux = lambda / omega(i);
    X(i) = (aux + std::sqrt(1 + aux * aux)) * std::exp(X(i));
  }
  
  return X % arma::sqrt(b / a);
}

// [[Rcpp::export]]
arma::vec rig(
    const int N, 
    const arma::vec mu, const double lambda) {
  
  arma::vec X(N);
  
  double nu, y, z;
  double muy, mu2lambda;
  
  for (int i = 0; i < N; ++i) {
    nu = R::rnorm(0, 1);
    y  = nu * nu;
    muy       = mu(i) * y;
    mu2lambda = mu(i) / (2 * lambda);
    X(i) = mu(i) + 
      muy * mu2lambda - 
      mu2lambda * std::sqrt(4 * muy * lambda + muy * muy);
    z = R::runif(0, 1);
    if (z > mu(i) / (mu(i) + X(i))) {
      X(i) = mu(i) * mu(i) / X(i);
    }
  }
  
  return X;
}

//// [[Rcpp::export]]
//arma::vec priorSpatialGP(
//  arma::vec v, 
//  double mean, double prec, double decay,
//  double Rsum, arma::mat Rinv
//  ) {
//  
//  const double na = 0;
//  const double nb = 0.000001;
//  
//  double ZtRZ = arma::as_scalar(Z.t() * Rinv * Z);
//  
//  // mean
//  V    = 1 / (prec * Rsum + nb); 
//  mean = V * (prec * arma::as_scalar(onen.t() * Rinv * v) + nb * na);
//  mean = R::rnorm(mean, sqrt(V));
//  
//  // prec
//  Z = v - mean;
//  prec = R::rgamma(n / 2 + ga, 1 / (ZtRZ / 2 + gb));
//  
//  // decay
//  decay_aux   = rtnorm(decay, sd, Umin, Umax);
//  Rinv_aux    = arma::inv_sympd(exp(- decay_aux * dist));
//  Rlogdet_aux = arma::log_det_sympd(Rinv_aux);
//  ZtRZ_aux    = arma::as_scalar(Z.t() * Rinv_aux * Z);  
//  A = (Rlogdet_aux - prec * ZtRZ_aux) / 2 - 
//    (Rlogdet - prec * ZtRZ) / 2 + 
//    dtnorm(decay, decay_aux, sd, Umin, Umax) - 
//   dtnorm(decay_aux, decay, sd, Umin, Umax);
//  if (log(R::runif(0, 1)) <= A) {
//    decay = decay_aux;
//    Rinv = Rinv_aux;
//    Rlogdet = Rlogdet_aux;
//    Rsum = arma::accu(Rinv);
//    ZtRZ = ZtRZ_aux;
//  }
//  
//  return { mean, prec, decay };
//}

// [[Rcpp::export]]
arma::mat iidMeanRcpp(
  const arma::vec Y, const arma::mat X, 
  const int N, const int k,
  arma::vec beta, double prec, 
  arma::mat keep,
  const int nSims, const int nThin, const int nBurnin, const int nReport) {
  
  // prior values
  const double na = 0;
  const double nb = 0.000001;
  const double ga = 0.0001;
  const double gb = 0.0001;
  
  // constants
  const arma::vec XtY = X.t() * Y;
  const arma::mat XtX = X.t() * X;
  
  // auxiliary
  arma::vec Z(N);
  
  // prior Mean M and Precision P of beta, PM = P * M
  const arma::mat P  = nb * arma::mat(k, k, arma::fill::eye);
  const arma::vec PM = P  * arma::vec(k, arma::fill::value(na));
  
  // full posterior Variance of beta
  arma::mat V(k, k);
  // full posterior parameters A B of prec
  double A = N / 2 + ga;
  double B;
  
  // iterations
  for (int b = 1 - nBurnin; b <= nSims; ++b) {

    if (b % nReport == 0) {
      Rcpp::Rcout << "Current iteration : " << b << "\n";
    }
      
    // beta
    V    = arma::inv_sympd(prec * XtX + P);
    beta = V * (prec * XtY + PM);
    beta += arma::chol(V, "lower") * arma::randn(k);
    
    // prec
    Z = Y - X * beta;
    B = arma::accu(Z % Z) / 2 + gb;
    prec = R::rgamma(A, 1 / B);
    
    // keep
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, k - 1)) = beta.t();
      keep(b / nThin - 1, k) = prec;
    }
  }
  
  return keep;
}

// [[Rcpp::export]]
arma::mat iidQuantileRcpp(
    double tau,
    const arma::vec Y, const arma::mat X, 
    const int N, const int k,
    arma::vec beta, double prec, 
    arma::mat keep,
    const int nSims, const int nThin, const int nBurnin, const int nReport) {
  
  // prior values
  const double na = 0;
  const double nb = 0.000001;
  const double ga = 0.0001;
  const double gb = 0.0001;
  
  // constants
  const double c1 = (1 - 2 * tau) / (tau * (1 - tau));
  const double c2 = tau * (1 - tau) / 2;
  const double c3 = 2 + c1 * c1 * c2;
  const double c4 = sqrt(c3 / c2);
  
  // latent
  arma::vec xi(N, arma::fill::ones);
  
  // auxiliary
  arma::vec Z = Y - X * beta - c1 * xi;
  arma::vec c2xiInv(N);
  arma::vec Xaux(k);
  
  // prior Mean M and Precision P of beta, PM = P * M
  const arma::mat P  = nb * arma::mat(k, k, arma::fill::eye);
  const arma::vec PM = P  * arma::vec(k, arma::fill::value(na));
  
  // full posterior Variance of beta
  arma::mat V(k, k);
  // full posterior parameters A B of prec
  double A = 3 * N / 2 + ga;
  double B;
  
  // iterations
  for (int b = 1 - nBurnin; b <= nSims; ++b) {
    
    if (b % nReport == 0) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    // xi
    Z += c1 * xi;
    xi = 1 / rig(N, c4 / abs(Z), prec * c3);
    c2xiInv = c2 / xi;
    Z -= c1 * xi;
    
    // beta
    Z += X * beta;
    V    = P;
    beta = PM;
    for (int i = 0; i < N; ++i) {
      Xaux = prec * c2xiInv(i) * X.row(i).t();
      V += Xaux * X.row(i);
      beta += Xaux * Z(i);
    }
    V    = arma::inv_sympd(V);
    beta = V * beta + arma::chol(V, "lower") * arma::randn(k);
    Z -= X * beta;
    
    // prec
    B = arma::accu(xi + c2xiInv % Z % Z / 2) + gb;
    prec = R::rgamma(A, 1 / B);
    
    // keep
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, k - 1)) = beta.t();
      keep(b / nThin - 1, k) = prec;
    }
  }
  
  return keep;
}

// [[Rcpp::export]]
arma::mat arMeanRcpp(
    const arma::vec Y, const arma::mat X, 
    const int N, const int k,
    arma::vec beta, double rho, double prec, 
    arma::mat keep,
    const int nSims, const int nThin, const int nBurnin, const int nReport) {
  
  // prior values
  const double na = 0;
  const double nb = 0.000001;
  const double ga = 0.0001;
  const double gb = 0.0001;
  
  // constants
  const arma::vec XtY = X.t() * Y;
  const arma::mat XtX = X.t() * X;
  
  // auxiliary
  arma::vec Z(N);
  
  // prior Mean M and Precision P of beta, PM = P * M
  const arma::mat P  = nb * arma::mat(k, k, arma::fill::eye);
  const arma::vec PM = P  * arma::vec(k, arma::fill::value(na));
  
  // full posterior Variance of beta
  arma::mat V(k, k);
  // full posterior parameters A B of prec
  double A = N / 2 + ga;
  double B;
  
  // iterations
  for (int b = 1 - nBurnin; b <= nSims; ++b) {
    
    if (b % nReport == 0) {
      Rcpp::Rcout << "Current iteration : " << b << "\n";
    }
    
    // beta
    V    = arma::inv_sympd(prec * XtX + P);
    beta = V * (prec * XtY + PM);
    beta += arma::chol(V, "lower") * arma::randn(k);
    
    // rho
    
    
    // prec
    Z = Y - X * beta;
    B = arma::accu(Z % Z) / 2 + gb;
    prec = R::rgamma(A, 1 / B);
    
    // keep
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, k - 1)) = beta.t();
      keep(b / nThin - 1, k) = prec;
      keep(b / nThin - 1, k + 1) = rho;
    }
  }
  
  return keep;
}

// [[Rcpp::export]]
double spTMeanRcpp() {

    return 1;
}

// [[Rcpp::export]]
double spTQuantileRcpp() {
  
  return 1;
}