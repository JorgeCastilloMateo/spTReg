#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// @description Computes the distance matrix.
// @param coords \eqn{n \times 2} coordinates matrix.
// [[Rcpp::export]]
arma::mat dist1(arma::mat coords) {
  int n = coords.n_rows;
  arma::mat D(n, n, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      D(i, j) = arma::norm(coords.row(i) - coords.row(j));
      D(j, i) = D(i, j);
    }
  }
  return D;
}

// @description Computes the squared distance matrix.
// @param coords \eqn{n \times m} coordinates matrix.
// [[Rcpp::export]]
arma::mat dist2(arma::mat coords) {
  int n = coords.n_rows;
  int m = coords.n_cols;
  arma::mat D(n, n, arma::fill::zeros);
  arma::rowvec v(m);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      v = coords.row(i) - coords.row(j);
      D(i, j) = arma::dot(v, v);
      D(j, i) = D(i, j);
    }
  }
  return D;
}

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
  const arma::vec mu, 
  const double lambda) {
  
  arma::vec X(N);
  double V, W, C, P1, Y;
  
  for (int i = 0; i < N; ++i) {
    V = R::rchisq(1);
    W = mu(i) * V;
    C = mu(i) / (2 * lambda);
    X(i) = mu(i) + C * 
      (W - std::sqrt(W * (4 * lambda + W)));
    P1 = mu(i) / (mu(i) + X(i));
    Y = R::runif(0, 1);
    if (Y > P1) {
      X(i) = mu(i) * mu(i) / X(i);
    }
  }
  
  return X;
}

// [[Rcpp::export]]
arma::vec ral(
  const int N, 
  const double theta, 
  const double w2Inv,
  const double sigma) {
  
  arma::vec _logU = -log(arma::randu(N));
  arma::vec X = (theta * _logU + sqrt(_logU / w2Inv) % arma::randn(N)) * sigma;
  
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

// iid mean model
// [[Rcpp::export]]
arma::mat iidMeanRcpp(
  const arma::vec Y, // data
  const arma::mat X,
  const arma::vec M, // prior values
  const arma::mat P,
  const double ga,
  const double gb,
  arma::vec beta,    // initial values
  double prec,
  const int N,       // constants
  const int k,
  arma::mat keep,    // replicates
  const int nSims,   // MCMC numbers
  const int nThin,
  const int nBurnin,
  const int nReport
) {
  
  // constants
  const arma::vec XtY = X.t() * Y;
  const arma::mat XtX = X.t() * X;
  
  // aux
  arma::vec e(N);
  
  // prior Mean M and Precision P of beta, PM = P * M
  const arma::vec PM = P * M;
  
  // full posterior Variance of beta
  arma::mat V(k, k);
  // full posterior parameters A B of prec
  double A = N / 2 + ga;
  double B;
  
  // iterations
  for (int b = 1 - nBurnin; b <= nSims; ++b) {

    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "Current iteration : " << b << "\n";
    }
      
    // beta
    V    = arma::inv_sympd(prec * XtX + P);
    beta = V * (prec * XtY + PM);
    beta += arma::chol(V, "lower") * arma::randn(k);
    
    // prec
    e = Y - X * beta;
    B = arma::accu(e % e) / 2 + gb;
    prec = R::rgamma(A, 1 / B);
    
    // keep
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, k - 1)) = beta.t();
      keep(b / nThin - 1, k) = prec;
      //keep(b / nThin - 1, arma::span(k + 1, k + N)) = (X * beta).t();
        //(X * beta + arma::randn(N) / sqrt(prec)).t();
      //keep(b / nThin - 1, arma::span(k + N + 1, k + 2 * N)) = 
      //  Y.t() - keep(b / nThin - 1, arma::span(k + 1, k + N));
    }
  }
  
  return keep;
}

// iid quantile model
// [[Rcpp::export]]
arma::mat iidQuantileRcpp(
  const double tau,  // quantile level
  const arma::vec Y, // data
  const arma::mat X,
  const arma::vec M, // prior values
  const arma::mat P,
  const double ga,
  const double gb,
  arma::vec beta,    // initial values
  double prec,
  const int N,       // constants
  const int k,
  arma::mat keep,    // replicates
  const int nSims,   // MCMC numbers
  const int nThin,
  const int nBurnin,
  const int nReport
) {
  
  // constants
  const double c1 = (1 - 2 * tau) / (tau * (1 - tau));
  const double c2 = tau * (1 - tau) / 2;
  const double c3 = 2 + c1 * c1 * c2;
  const double c4 = sqrt(c3 / c2);
  
  // latent
  arma::vec xi(N, arma::fill::ones);
  
  // aux
  arma::vec e = Y - X * beta - c1 * xi;
  arma::vec c2xiInv(N);
  arma::vec Xaux(k);
  
  // prior Mean M and Precision P of beta, PM = P * M
  const arma::vec PM = P * M;
  
  // full posterior Variance of beta
  arma::mat V(k, k);
  // full posterior parameters A B of prec
  double A = 3 * N / 2 + ga;
  double B;
  
  // iterations
  for (int b = 1 - nBurnin; b <= nSims; ++b) {
    
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "Current iteration : " << b << "\n";
    }
    
    // xi
    e += c1 * xi;
    xi = 1 / rig(N, c4 / abs(e), prec * c3);
    c2xiInv = c2 / xi;
    e -= c1 * xi;
    
    // beta
    e += X * beta;
    V    = P;
    beta = PM;
    for (int i = 0; i < N; ++i) {
      Xaux = prec * c2xiInv(i) * X.row(i).t();
      V    += Xaux * X.row(i);
      beta += Xaux * e(i);
    }
    V    = arma::inv_sympd(V);
    beta = V * beta + arma::chol(V, "lower") * arma::randn(k);
    e -= X * beta;
    
    // prec
    B = arma::accu(xi + c2xiInv % e % e / 2) + gb;
    prec = R::rgamma(A, 1 / B);
    
    // keep
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, k - 1)) = beta.t();
      keep(b / nThin - 1, k) = prec;
      //keep(b / nThin - 1, arma::span(k + 1, k + N)) = (X * beta).t();
        //(X * beta + ral(N, c1, c2, 1 / prec)).t();
      //keep(b / nThin - 1, arma::span(k + N + 1, k + 2 * N)) = 
      //  Y.t() - keep(b / nThin - 1, arma::span(k + 1, k + N));
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

// spatial mean model
// [[Rcpp::export]]
arma::mat spMeanRcpp(
    const arma::vec Y,    // data
    const arma::mat X,
    const arma::mat V,  
    const arma::mat dist, 
    const arma::vec M,    // prior values
    const arma::mat P,
    const double da,    
    const double db,    
    const double ga,
    const double gb,
    const double na,    
    const double nb,    
    arma::vec beta,       // initial values
    arma::mat betas,    
    double prec,
    arma::mat hp,       
    const int N,          // constants
    const int n,        
    const int k,
    const int r,        
    const arma::uvec s,  
    arma::mat keep,       // replicates
    const int nSims,      // MCMC numbers
    const int nThin,
    const int nBurnin,
    const int nReport
) {
  
  // constants
  const arma::mat Xt  = X.t();
  const arma::mat XtX = Xt * X;
  
  // aux
  arma::uvec ind(1);
  arma::vec Xb = X * beta;
  arma::vec e  = Y - Xb;
  for (int m = 0; m < r; ++m) {
    ind = m;
    e -= V.col(m) % betas(s, ind);
  }
  
  // aux GP
  arma::vec onen(n, arma::fill::ones);
  int count0, count1; 
  int Nn = N / n;
  
  hp.row(0).zeros(); //temporal
  
  // aux decay
  arma::vec accept(r, arma::fill::zeros);
  int total = 0;
  arma::vec ratio(r);
  arma::vec sd(r, arma::fill::ones);
  arma::vec lsd(r, arma::fill::zeros); // log(sd);
  
  arma::cube R(n, n, r);
  arma::vec Rlogdet(r);
  arma::vec oneRone(r);
  
  for (int m = 0; m < r; ++m) {
    R.slice(m) = arma::inv_sympd(exp(- hp(2, m) * dist));
    Rlogdet(m) = arma::log_det_sympd(R.slice(m));
    oneRone(m) = arma::accu(R.slice(m));
  }
  
  double decay_aux;
  double ldecay_aux = 0;
  arma::vec ldecay = log(hp.row(2).t());
  arma::mat R_aux(n, n);
  double Rlogdet_aux;
  
  arma::vec vn(n);
  double vtRv, vtRv_aux;
  double alpha = 0;
  
  // prior Mean M and Precision P of beta, PM = P * M
  const arma::vec PM = P * M;
  
  // full posterior Variance of beta
  arma::mat Vk(k, k);
  // full posterior Variance of betas
  arma::mat Vn(n, n);
  // full posterior parameters chi delta of mu
  double chi;
  double delta;
  // full posterior parameters A B of prec, and C - of prec (betas)
  double A = N / 2 + ga;
  double B;
  double C = n / 2 + ga;
  
  // iterations
  for (int b = 1 - nBurnin; b <= nSims; ++b) {
    
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "Current iteration : " << b << "\n";
    }
    
    // beta
    e += Xb;
    Vk   = arma::inv_sympd(prec * XtX + P);
    beta = Vk * (prec * Xt * e + PM);
    beta += arma::chol(Vk, "lower") * arma::randn(k);
    Xb = X * beta;
    e -= Xb;
    
    // betas m = 1,...,r
    for (int m = 0; m < r; ++m) {
      ind = m;
      e += V.col(m) % betas(s, ind);
      Vn = hp(1, m) * R.slice(m);
      betas.col(m) = Vn * onen * hp(0, m);
      count1 = 0;
      for (int i = 0; i < n; ++i) {
        count0 = count1;
        count1 += Nn;
        Vn(i, i)    += prec * arma::accu(V(arma::span(count0, count1 - 1), m) % V(arma::span(count0, count1 - 1), m));
        betas(i, m) += prec * arma::accu(e(arma::span(count0, count1 - 1))    % V(arma::span(count0, count1 - 1), m));
      }
      Vn = arma::inv_sympd(Vn);
      betas.col(m) = Vn * betas.col(m) + arma::chol(Vn, "lower") * arma::randn(n);
      betas.col(m) -= Vn * onen * (arma::accu(betas.col(m)) / arma::accu(Vn));
      e -= V.col(m) % betas(s, ind);
      
      // mu 
      //delta = 1 / (oneRone(m) * hp(1, m) + nb);
      //chi   = arma::as_scalar(onen.t() * R.slice(m) * betas.col(m)) * hp(1, m) + na * nb;
      //hp(0, m) = R::rnorm(delta * chi, sqrt(delta));
      
      // decay
      ldecay_aux  = R::rnorm(ldecay(m), sd(m));
      decay_aux   = exp(ldecay_aux);
      R_aux       = arma::inv_sympd(exp(- decay_aux * dist));
      Rlogdet_aux = arma::log_det_sympd(R_aux);
      vn       = betas.col(m) - hp(0, m);
      vtRv_aux = arma::as_scalar(vn.t() * R_aux * vn);
      vtRv     = arma::as_scalar(vn.t() * R.slice(m) * vn);
      alpha = 
        (Rlogdet_aux - hp(1, m) * vtRv_aux) / 2 + 
        da * ldecay_aux - db * decay_aux - 
        ((Rlogdet(m) - hp(1, m) * vtRv) / 2 +
        da * ldecay(m) - db * hp(2, m));
      if (log(R::runif(0, 1)) < alpha) {
        ++accept(m);
        hp(2, m) = decay_aux;
        ldecay(m) = ldecay_aux;
        R.slice(m) = R_aux;
        Rlogdet(m) = Rlogdet_aux;
        oneRone(m) = arma::accu(R);
        vtRv = vtRv_aux;
      }
      
      // prec
      hp(1, m) = R::rgamma(C, 
        1 / (vtRv / 2 + gb));
    }
  
    // tune sd of the proposal for decay
    if (b == 0) {
      accept.zeros();
      total = 0;
    } else if ((b < 1) && (++total % 25 == 0)) {
      ratio = accept / total;
      for (int m = 0; m < r; ++m) {
        if (ratio(m) > 0.33) {
          lsd(m) += 1 / sqrt(total / 25);
        } else {
          lsd(m) -= 1 / sqrt(total / 25);
        }
        sd(m) = exp(lsd(m));
      }
    }

    // prec
    B = arma::accu(e % e) / 2 + gb;
    prec = R::rgamma(A, 1 / B);
    
    // keep
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, k - 1)) = beta.t();
      keep(b / nThin - 1, k) = prec;
      for (int m = 0; m < r; ++m) {
        keep(b / nThin - 1, arma::span(k + 1 + m * (n + 3), k + n + m * (n + 3))) = betas.col(m).t();
        keep(b / nThin - 1, arma::span(k + 1 + n + m * (n + 3), k + (m + 1) * (n + 3))) = hp.col(m).t();
      }
    }
  }
  
  return keep;
}

// spatial quantile model
// [[Rcpp::export]]
arma::mat spQuantileRcpp(
    const double tau,     // quantile level
    const arma::vec Y,    // data
    const arma::mat X,
    const arma::mat V,  
    const arma::mat dist, 
    const arma::vec M,    // prior values
    const arma::mat P,
    const double da,    
    const double db,    
    const double ga,
    const double gb,
    const double na,    
    const double nb,    
    arma::vec beta,       // initial values
    arma::mat betas,    
    double prec,
    arma::mat hp,       
    const int N,          // constants
    const int n,        
    const int k,
    const int r,        
    const arma::uvec s,  
    arma::mat keep,       // replicates
    const int nSims,      // MCMC numbers
    const int nThin,
    const int nBurnin,
    const int nReport
) {
  
  // constants
  const double c1 = (1 - 2 * tau) / (tau * (1 - tau));
  const double c2 = tau * (1 - tau) / 2;
  const double c3 = 2 + c1 * c1 * c2;
  const double c4 = sqrt(c3 / c2);
  
  // latent
  arma::vec xi(N, arma::fill::ones);
  
  // aux
  arma::uvec ind(1);
  arma::vec Xb = X * beta;
  arma::vec e  = Y - Xb - c1 * xi;
  for (int m = 0; m < r; ++m) {
    ind = m;
    e -= V.col(m) % betas(s, ind);
  }
  
  arma::vec c2xiInv(N);
  arma::vec Xaux(k);
  
  // aux GP
  arma::vec onen(n, arma::fill::ones);
  int count0, count1; 
  int Nn = N / n;
  
  hp.row(0).zeros(); //temporal
  
  // aux decay
  arma::vec accept(r, arma::fill::zeros);
  int total = 0;
  arma::vec ratio(r);
  arma::vec sd(r, arma::fill::ones);
  arma::vec lsd(r, arma::fill::zeros); // log(sd);
  
  arma::cube R(n, n, r);
  arma::vec Rlogdet(r);
  arma::vec oneRone(r);
  
  for (int m = 0; m < r; ++m) {
    R.slice(m) = arma::inv_sympd(exp(- hp(2, m) * dist));
    Rlogdet(m) = arma::log_det_sympd(R.slice(m));
    oneRone(m) = arma::accu(R.slice(m));
  }
  
  double decay_aux;
  double ldecay_aux = 0;
  arma::vec ldecay = log(hp.row(2).t());
  arma::mat R_aux(n, n);
  double Rlogdet_aux;
  
  arma::vec vn(n);
  double vtRv, vtRv_aux;
  double alpha = 0;
  
  // prior Mean M and Precision P of beta, PM = P * M
  const arma::vec PM = P * M;
  
  // full posterior Variance of beta
  arma::mat Vk(k, k);
  // full posterior Variance of betas
  arma::mat Vn(n, n);
  // full posterior parameters chi delta of mu
  double chi;
  double delta;
  // full posterior parameters A B of prec, and C - of prec (betas)
  double A = 3 * N / 2 + ga;
  double B;
  double C = n / 2 + ga;
  
  // iterations
  for (int b = 1 - nBurnin; b <= nSims; ++b) {
    
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "Current iteration : " << b << "\n";
    }
    
    // xi
    e += c1 * xi;
    xi = 1 / rig(N, c4 / abs(e), prec * c3);
    c2xiInv = c2 / xi;
    e -= c1 * xi;
    
    // beta
    e += Xb;
    Vk   = P;
    beta = PM;
    for (int i = 0; i < N; ++i) {
      Xaux = prec * c2xiInv(i) * X.row(i).t();
      Vk   += Xaux * X.row(i);
      beta += Xaux * e(i);
    }
    Vk   = arma::inv_sympd(Vk);
    beta = Vk * beta + arma::chol(Vk, "lower") * arma::randn(k);
    Xb = X * beta;
    e -= Xb;
    
    // betas m = 1,...,r
    for (int m = 0; m < r; ++m) {
      ind = m;
      e += V.col(m) % betas(s, ind);
      Vn = hp(1, m) * R.slice(m);
      betas.col(m) = Vn * onen * hp(0, m);
      count1 = 0;
      for (int i = 0; i < n; ++i) {
        count0 = count1;
        count1 += Nn;
        Vn(i, i)    += prec * arma::accu(V(arma::span(count0, count1 - 1), m) % 
                                         V(arma::span(count0, count1 - 1), m) %
                                   c2xiInv(arma::span(count0, count1 - 1)));
        betas(i, m) += prec * arma::accu(e(arma::span(count0, count1 - 1))    % 
                                         V(arma::span(count0, count1 - 1), m) %
                                   c2xiInv(arma::span(count0, count1 - 1)));
      }
      Vn = arma::inv_sympd(Vn);
      betas.col(m) = Vn * betas.col(m) + arma::chol(Vn, "lower") * arma::randn(n);
      betas.col(m) -= Vn * onen * (arma::accu(betas.col(m)) / arma::accu(Vn));
      e -= V.col(m) % betas(s, ind);
      
      // mu 
      //delta = 1 / (oneRone(m) * hp(1, m) + nb);
      //chi   = arma::as_scalar(onen.t() * R.slice(m) * betas.col(m)) * hp(1, m) + na * nb;
      //hp(0, m) = R::rnorm(delta * chi, sqrt(delta));
      
      // decay
      ldecay_aux  = R::rnorm(ldecay(m), sd(m));
      decay_aux   = exp(ldecay_aux);
      R_aux       = arma::inv_sympd(exp(- decay_aux * dist));
      Rlogdet_aux = arma::log_det_sympd(R_aux);
      vn       = betas.col(m) - hp(0, m);
      vtRv_aux = arma::as_scalar(vn.t() * R_aux * vn);
      vtRv     = arma::as_scalar(vn.t() * R.slice(m) * vn);
      alpha = 
        (Rlogdet_aux - hp(1, m) * vtRv_aux) / 2 + 
        da * ldecay_aux - db * decay_aux - 
        ((Rlogdet(m) - hp(1, m) * vtRv) / 2 +
        da * ldecay(m) - db * hp(2, m));
      if (log(R::runif(0, 1)) < alpha) {
        ++accept(m);
        hp(2, m) = decay_aux;
        ldecay(m) = ldecay_aux;
        R.slice(m) = R_aux;
        Rlogdet(m) = Rlogdet_aux;
        oneRone(m) = arma::accu(R);
        vtRv = vtRv_aux;
      }
      
      // prec
      hp(1, m) = R::rgamma(C, 
         1 / (vtRv / 2 + gb));
    }
    
    // tune sd of the proposal for decay
    if (b == 0) {
      accept.zeros();
      total = 0;
    } else if ((b < 1) && (++total % 25 == 0)) {
      ratio = accept / total;
      for (int m = 0; m < r; ++m) {
        if (ratio(m) > 0.33) {
          lsd(m) += 1 / sqrt(total / 25);
        } else {
          lsd(m) -= 1 / sqrt(total / 25);
        }
        sd(m) = exp(lsd(m));
      }
    }
    
    // prec
    B = arma::accu(xi + c2xiInv % e % e / 2) + gb;
    prec = R::rgamma(A, 1 / B);
    
    // keep
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, k - 1)) = beta.t();
      keep(b / nThin - 1, k) = prec;
      for (int m = 0; m < r; ++m) {
        keep(b / nThin - 1, arma::span(k + 1 + m * (n + 3), k + n + m * (n + 3))) = betas.col(m).t();
        keep(b / nThin - 1, arma::span(k + 1 + n + m * (n + 3), k + (m + 1) * (n + 3))) = hp.col(m).t();
      }
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