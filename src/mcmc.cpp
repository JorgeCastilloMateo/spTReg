#include "utils.h"

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// iid mean model
// [[Rcpp::export]]
Rcpp::List iidMeanRcpp(
  arma::vec Y,        // data
  const arma::mat& X,
  const arma::vec& M, // priors
  const arma::mat& P,
  const double ga,
  const double gb,
  arma::vec beta,     // initial
  double prec,
  const int N,        // constants
  const int p,
  const int nSims,    // MCMC numbers
  const int nThin,
  const int nBurnin,
  const int nReport
) {
  
  // missing index
  int nKeep = nSims / nThin;
  const arma::uvec missing_idx = arma::find_nonfinite(Y);
  const arma::uword missing_n = missing_idx.n_elem;
  arma::mat keep_Y(nKeep, missing_n);
  const double mean_Y = arma::mean(Y.elem(arma::find_finite(Y)));
  Y.elem(missing_idx).fill(mean_Y);
  
  // constants
  arma::vec XtY = X.t() * Y;
  const arma::mat XtX = X.t() * X;
  const arma::vec PM = P * M;
  
  // residual
  arma::vec e(N);
  // full posterior Precision rhs of beta
  arma::mat Q(p, p);
  arma::vec b(p);
  // full posterior parameters A B of prec
  double A = 0.5 * N + ga;
  double B;
  
  // save
  int save_idx = 0;
  int nCols = p + 1;
  arma::mat keep(nKeep, nCols);

  // time
  auto start_time = std::chrono::steady_clock::now();
  
  // iterations
  for (int iter = 1 - nBurnin; iter <= nSims; ++iter) {

    reportProgress(iter, nBurnin, nSims, nReport, 
                   (iter > 0 ? iter / nThin : 0), start_time);
      
    // beta
    Q = P  + prec * XtX;
    b = PM + prec * XtY;
    beta = RandomMultiNormalC(Q, b);
    
    // prec
    e = Y - X * beta;
    B = 0.5 * arma::dot(e, e) + gb;
    prec = R::rgamma(A, 1.0 / B);
    
    // Y missing
    if (missing_n > 0) {
      XtY -= X.rows(missing_idx).t() * Y.elem(missing_idx);
      Y.elem(missing_idx) = X.rows(missing_idx) * beta + 
        arma::randn(missing_n) / std::sqrt(prec);
      XtY += X.rows(missing_idx).t() * Y.elem(missing_idx);
      
      if (iter > 0 && iter % nThin == 0) {
        keep_Y.row(save_idx) = Y.elem(missing_idx).t();
      }
    }
    
    // keep
    if (iter > 0 && iter % nThin == 0) {
      keep(save_idx, arma::span(0, p - 1)) = beta.t();
      keep(save_idx, p) = prec;
      ++save_idx;
    }
  }
  
  if (missing_n == 0) {
    return Rcpp::List::create(
      Rcpp::Named("params") = keep
    );
  } else {
    return Rcpp::List::create(
      Rcpp::Named("params") = keep,
      Rcpp::Named("missing") = keep_Y
    );
  }
}

// iid quantile model
// [[Rcpp::export]]
Rcpp::List iidQuantileRcpp(
  const double tau,   // quantile level
  arma::vec Y,        // data
  const arma::mat& X,
  const arma::vec& M, // priors
  const arma::mat& P,
  const double ga,
  const double gb,
  arma::vec beta,     // initial
  double prec,
  const int N,        // constants
  const int p,
  const int nSims,    // MCMC numbers
  const int nThin,
  const int nBurnin,
  const int nReport,
  const bool parallel,// parallel
  const int nThreads
) {
  
  // missing index
  int nKeep = nSims / nThin;
  const arma::uvec missing_idx = arma::find_nonfinite(Y);
  const arma::uword missing_n = missing_idx.n_elem;
  arma::vec mm(missing_n), ss(missing_n);
  arma::mat keep_Y(nKeep, missing_n);
  const double mean_Y = arma::mean(Y.elem(arma::find_finite(Y)));
  Y.elem(missing_idx).fill(mean_Y);
  
  // constants
  const double c1 = (1.0 - 2.0 * tau) / (tau * (1.0 - tau));
  const double c2 = 0.5 * tau * (1.0 - tau);
  const double c3 = 2.0 + c1 * c1 * c2;
  const double c4 = sqrt(c3 / c2);
  const arma::vec PM = P * M;
  
  // aux
  arma::vec c2dxi(N);
  arma::vec Xaux(p);
  // latent
  arma::vec xi(N, arma::fill::ones);
  // residual
  arma::vec Xb = X * beta;
  arma::vec e = Y - Xb - c1 * xi;
  // full posterior Precision rhs of beta
  arma::mat Q(p, p);
  arma::vec b(p);
  // full posterior parameters A B of prec
  double A = 1.5 * N + ga;
  double B;
  
  // save
  int save_idx = 0;
  int nCols = p + 1;
  arma::mat keep(nKeep, nCols);
  
  // time
  auto start_time = std::chrono::steady_clock::now();
  
  // iterations
  for (int iter = 1 - nBurnin; iter <= nSims; ++iter) {
    
    reportProgress(iter, nBurnin, nSims, nReport, 
                   (iter > 0 ? iter / nThin : 0), start_time);
    
    // xi
    e += c1 * xi;
    xi = 1.0 / rig(N, c4 / arma::abs(e), prec * c3, parallel, nThreads);
    c2dxi = c2 / xi;
    e -= c1 * xi;
    
    // beta
    e += Xb;
    Q = P;
    b = PM;
    for (int i = 0; i < N; ++i) {
      Xaux = prec * c2dxi(i) * X.row(i).t();
      Q += Xaux * X.row(i);
      b += Xaux * e(i);
    }
    beta = RandomMultiNormalC(Q, b);
    Xb = X * beta;
    e -= Xb;
    
    // prec
    B = arma::accu(xi + 0.5 * c2dxi % arma::square(e)) + gb;
    prec = R::rgamma(A, 1.0 / B);
    
    // Y missing
    if (missing_n > 0) {
      e.elem(missing_idx) -= Y.elem(missing_idx);
      mm = X.rows(missing_idx) * beta + c1 * xi.elem(missing_idx);
      ss = arma::sqrt(xi.elem(missing_idx) / (prec * c2));
      Y.elem(missing_idx) = mm + ss % arma::randn(missing_n);
      e.elem(missing_idx) += Y.elem(missing_idx);
      
      if (iter > 0 && iter % nThin == 0) {
        keep_Y.row(save_idx) = Y.elem(missing_idx).t();
      }
    }
    
    // keep
    if (iter > 0 && iter % nThin == 0) {
      keep(save_idx, arma::span(0, p - 1)) = beta.t();
      keep(save_idx, p) = prec;
      ++save_idx;
    }
  }
  
  if (missing_n == 0) {
    return Rcpp::List::create(
      Rcpp::Named("params") = keep
    );
  } else {
    return Rcpp::List::create(
      Rcpp::Named("params") = keep,
      Rcpp::Named("missing") = keep_Y
    );
  }
}

// spatial mean model
// [[Rcpp::export]]
Rcpp::List spMeanRcpp(
    arma::vec Y,           // data
    const arma::mat& X,
    const arma::mat& U,
    const arma::mat& V,
    const arma::mat& dist, 
    const arma::vec& M,    // priors
    const arma::mat& P,
    const double da,    
    const double db,    
    const double ga,
    const double gb,
    const double na,    
    const double nb,    
    arma::vec beta,        // initial
    arma::mat gamma,
    arma::mat alpha,
    double prec,
    arma::mat hp,       
    const int N,           // constants
    const int n,
    const int p,
    const int q,
    const int r,        
    const arma::uvec& s,
    const int nSims,       // MCMC numbers
    const int nThin,
    const int nBurnin,
    const int nReport
) {
  
  // missing index
  int nKeep = nSims / nThin;
  const arma::uvec missing_idx = arma::find_nonfinite(Y);
  const arma::uword missing_n = missing_idx.n_elem;
  arma::mat keep_Y(nKeep, missing_n);
  const double mean_Y = arma::mean(Y.elem(arma::find_finite(Y)));
  Y.elem(missing_idx).fill(mean_Y);
  
  // constants
  const arma::mat Xt  = X.t();
  const arma::mat XtX = Xt * X;
  const arma::vec PM = P * M;
  
  // residual
  arma::vec Xb = X * beta;
  arma::vec e = Y - Xb;
  arma::vec alpha_m(n);
  arma::vec V_m(N);
  arma::uvec s_missing(missing_n);
  for (int m = 0; m < r; ++m) {
    alpha_m = alpha.col(m);
    V_m = V.col(m);
    e -= V_m % alpha_m.elem(s);
  }

  // aux GP
  arma::vec onen(n, arma::fill::ones);
  int Ndn = N / n;
  
  std::vector<arma::uvec> s_group(n);
  for (int i = 0; i < n; ++i) {
    s_group[i] = arma::regspace<arma::uvec>(i*Ndn, i*Ndn + Ndn - 1);
  }
  
  arma::vec V_block(Ndn), e_block(Ndn);
  
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
  double ALPHA = 0;
  
  // full posterior Precision rhs of beta
  arma::mat Qp(p, p);
  arma::vec bp(p);
  // full posterior Precision rhs of alpha
  arma::mat Qn(n, n);
  arma::vec bn(n);
  // full posterior parameters chi delta of mu
  //double chi;
  //double delta;
  // full posterior parameters A B of prec, and C D of prec (alpha)
  double A = 0.5 * N + ga;
  double B;
  double C = 0.5 * n + ga;
  double D;
  
  // save
  int save_idx = 0;
  int nCols = p + 1 + r * (n + 3);
  arma::mat keep(nKeep, nCols);

  // time
  auto start_time = std::chrono::steady_clock::now();
  
  // iterations
  for (int iter = 1 - nBurnin; iter <= nSims; ++iter) {
    
    reportProgress(iter, nBurnin, nSims, nReport, 
                   (iter > 0 ? iter / nThin : 0), start_time);
    
    // beta
    e += Xb;
    Qp = P  + prec * XtX;
    bp = PM + prec * Xt * e;
    beta = RandomMultiNormalC(Qp, bp);
    Xb = X * beta;
    e -= Xb;
    
    // gamma t = 0,1,...,T-1,T
    //if (q > 0) {
    // t = 0
    
    // t = 1,...,T-1
    //for (int t = 1; t < T; ++t) {
      
    //}
    // t = T
    
    //}
    
    // alpha m = 1,...,r
    if (r > 0) {
    for (int m = 0; m < r; ++m) {
      alpha_m = alpha.col(m);
      V_m = V.col(m);
      e += V_m % alpha_m.elem(s);
      Qn = hp(1, m) * R.slice(m);
      bn = Qn * (onen * hp(0, m));
      for (int i = 0; i < n; ++i) {
        V_block = V_m.elem(s_group[i]);
        e_block = e.elem(s_group[i]);
        Qn(i,i) += prec * arma::dot(V_block, V_block);
        bn(i)   += prec * arma::dot(V_block, e_block);
      }
      alpha_m = RandomMultiNormalC(Qn, bn);
      //alpha_m -= Vn * onen * (arma::accu(alpha_m) / arma::accu(Vn));
      alpha.col(m) = alpha_m;
      e -= V_m % alpha_m.elem(s);
      
      // mu 
      //delta = 1 / (oneRone(m) * hp(1, m) + nb);
      //chi   = arma::as_scalar(onen.t() * R.slice(m) * alpha.col(m)) * hp(1, m) + na * nb;
      //hp(0, m) = R::rnorm(delta * chi, sqrt(delta));
      
      // decay
      ldecay_aux  = R::rnorm(ldecay(m), sd(m));
      decay_aux   = exp(ldecay_aux);
      R_aux       = arma::inv_sympd(exp(- decay_aux * dist));
      Rlogdet_aux = arma::log_det_sympd(R_aux);
      vn       = alpha.col(m) - hp(0, m);
      vtRv_aux = arma::as_scalar(vn.t() * R_aux * vn);
      vtRv     = arma::as_scalar(vn.t() * R.slice(m) * vn);
      ALPHA = 
        (Rlogdet_aux - hp(1, m) * vtRv_aux) / 2 + 
        da * ldecay_aux - db * decay_aux - 
        ((Rlogdet(m) - hp(1, m) * vtRv) / 2 +
        da * ldecay(m) - db * hp(2, m));
      if (log(R::runif(0, 1)) < ALPHA) {
        ++accept(m);
        hp(2, m) = decay_aux;
        ldecay(m) = ldecay_aux;
        R.slice(m) = R_aux;
        Rlogdet(m) = Rlogdet_aux;
        oneRone(m) = arma::accu(R);
        vtRv = vtRv_aux;
      }
      
      // prec
      D = 0.5 * vtRv + gb;
      hp(1, m) = R::rgamma(C, 1.0 / D);
    }
    }
  
    // tune sd of the proposal for decay
    if (iter == 0) {
      accept.zeros();
      total = 0;
    } else if ((iter < 1) && (++total % 25 == 0)) {
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

    // w t=1,...,T l=1,...,L
    //w
    //ffbs(
    //  const arma::mat& Y,       // n x T : observations
    //  const double rho,         // AR parameter
    //  const arma::mat& Sigma_w, // state conditional covariance
    //  const Rcpp::List& D_list  // list of length T with vec of diag matrices
    //)
    
    // prec
    B = 0.5 * arma::dot(e, e) + gb;
    prec = R::rgamma(A, 1.0 / B);
    
    // Y missing
    if (missing_n > 0) {
      e.elem(missing_idx) -= Y.elem(missing_idx);
      Y.elem(missing_idx) = X.rows(missing_idx) * beta + 
        arma::randn(missing_n) / std::sqrt(prec);
      if (r > 0) {
        for (int m = 0; m < r; ++m) {
          alpha_m = alpha.col(m);
          V_m = V.col(m);
          s_missing = s.elem(missing_idx);
          Y.elem(missing_idx) += V_m.elem(missing_idx) % alpha_m.elem(s_missing);
        }
      }
      e.elem(missing_idx) += Y.elem(missing_idx);
      
      if (iter > 0 && iter % nThin == 0) {
        keep_Y.row(save_idx) = Y.elem(missing_idx).t();
      }
    }
    
    // keep
    if (iter > 0 && iter % nThin == 0) {
      keep(save_idx, arma::span(0, p - 1)) = beta.t();
      keep(save_idx, p) = prec;
      if (r > 0) {
      for (int m = 0; m < r; ++m) {
        keep(save_idx, arma::span(p + 1 + m * (n + 3), p + n + m * (n + 3))) = alpha.col(m).t();
        keep(save_idx, arma::span(p + 1 + n + m * (n + 3), p + (m + 1) * (n + 3))) = hp.col(m).t();
      }
      }
      ++save_idx;
    }
  }
  
  if (missing_n == 0) {
    return Rcpp::List::create(
      Rcpp::Named("params") = keep
    );
  } else {
    return Rcpp::List::create(
      Rcpp::Named("params") = keep,
      Rcpp::Named("missing") = keep_Y
    );
  }
}

// spatial quantile model
// [[Rcpp::export]]
Rcpp::List spQuantileRcpp(
    const double tau,      // quantile level
    arma::vec Y,           // data
    const arma::mat& X,
    const arma::mat& V,  
    const arma::mat& dist, 
    const arma::vec& M,    // prior values
    const arma::mat& P,
    const double da,    
    const double db,    
    const double ga,
    const double gb,
    const double na,    
    const double nb,    
    arma::vec beta,        // initial
    arma::mat alpha,    
    double prec,
    arma::mat hp,       
    const int N,          // constants
    const int n,        
    const int p,
    const int r,        
    const arma::uvec& s,  
    const int nSims,      // MCMC numbers
    const int nThin,
    const int nBurnin,
    const int nReport,
    const bool parallel,  // parallel
    const int nThreads
) {
  
  // missing index
  int nKeep = nSims / nThin;
  const arma::uvec missing_idx = arma::find_nonfinite(Y);
  const arma::uword missing_n = missing_idx.n_elem;
  arma::vec mm(missing_n), ss(missing_n);
  arma::mat keep_Y(nKeep, missing_n);
  const double mean_Y = arma::mean(Y.elem(arma::find_finite(Y)));
  Y.elem(missing_idx).fill(mean_Y);
  
  // constants
  const double c1 = (1 - 2 * tau) / (tau * (1 - tau));
  const double c2 = tau * (1 - tau) / 2;
  const double c3 = 2 + c1 * c1 * c2;
  const double c4 = sqrt(c3 / c2);
  const arma::vec PM = P * M;
  
  // aux
  arma::vec c2dxi(N);
  arma::vec Xaux(p);
  // latent
  arma::vec xi(N, arma::fill::ones);
  
  // residual
  arma::vec Xb = X * beta;
  arma::vec e  = Y - Xb - c1 * xi;
  arma::vec alpha_m(n);
  arma::vec V_m(N);
  arma::uvec s_missing(missing_n);
  for (int m = 0; m < r; ++m) {
    alpha_m = alpha.col(m);
    V_m = V.col(m);
    e -= V_m % alpha_m.elem(s);
  }
  
  // aux GP
  arma::vec onen(n, arma::fill::ones);
  int Ndn = N / n;
  
  std::vector<arma::uvec> s_group(n);
  for (int i = 0; i < n; ++i) {
    s_group[i] = arma::regspace<arma::uvec>(i*Ndn, i*Ndn + Ndn - 1);
  }
  
  arma::vec V_block(Ndn), e_block(Ndn);
  arma::vec c2dxi_block(Ndn), V_c2dxi_block(Ndn);
  
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
  double ALPHA = 0;
  
  // full posterior Precision rhs of beta
  arma::mat Qp(p, p);
  arma::vec bp(p);
  // full posterior Precision rhs of alpha
  arma::mat Qn(n, n);
  arma::vec bn(n);
  // full posterior parameters chi delta of mu
  //double chi;
  //double delta;
  // full posterior parameters A B of prec, and C D of prec (alpha)
  double A = 1.5 * N + ga;
  double B;
  double C = 0.5 * n + ga;
  double D;
  
  // save
  int save_idx = 0;
  int nCols = p + 1 + r * (n + 3);
  arma::mat keep(nKeep, nCols);
  
  // time
  auto start_time = std::chrono::steady_clock::now();
  
  // iterations
  for (int iter = 1 - nBurnin; iter <= nSims; ++iter) {
    
    reportProgress(iter, nBurnin, nSims, nReport, 
                   (iter > 0 ? iter / nThin : 0), start_time);
    
    // xi
    e += c1 * xi;
    xi = 1.0 / rig(N, c4 / arma::abs(e), prec * c3, parallel, nThreads);
    c2dxi = c2 / xi;
    e -= c1 * xi;
    
    // beta
    e += Xb;
    Qp = P;
    bp = PM;
    for (int i = 0; i < N; ++i) {
      Xaux = prec * c2dxi(i) * X.row(i).t();
      Qp += Xaux * X.row(i);
      bp += Xaux * e(i);
    }
    beta = RandomMultiNormalC(Qp, bp);
    Xb = X * beta;
    e -= Xb;

    // alpha m = 1,...,r
    if (r > 0) {
    for (int m = 0; m < r; ++m) {
      alpha_m = alpha.col(m);
      V_m = V.col(m);
      e += V_m % alpha_m.elem(s);
      Qn = hp(1, m) * R.slice(m);
      bn = Qn * (onen * hp(0, m));
      for (int i = 0; i < n; ++i) {
        V_block = V_m.elem(s_group[i]);
        e_block = e.elem(s_group[i]);
        c2dxi_block = c2dxi.elem(s_group[i]);
        V_c2dxi_block = V_block % c2dxi_block;
        Qn(i,i) += prec * arma::accu(V_c2dxi_block % V_block);
        bn(i)   += prec * arma::accu(V_c2dxi_block % e_block);
      }
      alpha_m = RandomMultiNormalC(Qn, bn);
      //alpha_m -= Vn * onen * (arma::accu(alpha_m) / arma::accu(Vn));
      alpha.col(m) = alpha_m;
      e -= V_m % alpha_m.elem(s);

      // mu 
      //delta = 1 / (oneRone(m) * hp(1, m) + nb);
      //chi   = arma::as_scalar(onen.t() * R.slice(m) * alpha.col(m)) * hp(1, m) + na * nb;
      //hp(0, m) = R::rnorm(delta * chi, sqrt(delta));
      
      // decay
      ldecay_aux  = R::rnorm(ldecay(m), sd(m));
      decay_aux   = exp(ldecay_aux);
      R_aux       = arma::inv_sympd(exp(- decay_aux * dist));
      Rlogdet_aux = arma::log_det_sympd(R_aux);
      vn       = alpha.col(m) - hp(0, m);
      vtRv_aux = arma::as_scalar(vn.t() * R_aux * vn);
      vtRv     = arma::as_scalar(vn.t() * R.slice(m) * vn);
      ALPHA = 
        (Rlogdet_aux - hp(1, m) * vtRv_aux) / 2 + 
        da * ldecay_aux - db * decay_aux - 
        ((Rlogdet(m) - hp(1, m) * vtRv) / 2 +
        da * ldecay(m) - db * hp(2, m));
      if (log(R::runif(0, 1)) < ALPHA) {
        ++accept(m);
        hp(2, m) = decay_aux;
        ldecay(m) = ldecay_aux;
        R.slice(m) = R_aux;
        Rlogdet(m) = Rlogdet_aux;
        oneRone(m) = arma::accu(R);
        vtRv = vtRv_aux;
      }
      
      // prec
      D = 0.5 * vtRv + gb;
      hp(1, m) = R::rgamma(C, 1.0 / D);
    }
    }
    
    // tune sd of the proposal for decay
    if (iter == 0) {
      accept.zeros();
      total = 0;
    } else if ((iter < 1) && (++total % 25 == 0)) {
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
    B = arma::accu(xi + 0.5 * c2dxi % arma::square(e)) + gb;
    prec = R::rgamma(A, 1.0 / B);
    
    // Y missing
    if (missing_n > 0) {
      e.elem(missing_idx) -= Y.elem(missing_idx);
      mm = X.rows(missing_idx) * beta + c1 * xi.elem(missing_idx);
      ss = arma::sqrt(xi.elem(missing_idx) / (prec * c2));
      Y.elem(missing_idx) = mm + ss % arma::randn(missing_n);
      if (r > 0) {
        for (int m = 0; m < r; ++m) {
          alpha_m = alpha.col(m);
          V_m = V.col(m);
          s_missing = s.elem(missing_idx);
          Y.elem(missing_idx) += V_m.elem(missing_idx) % alpha_m.elem(s_missing);
        }
      }
      e.elem(missing_idx) += Y.elem(missing_idx);
      
      if (iter > 0 && iter % nThin == 0) {
        keep_Y.row(save_idx) = Y.elem(missing_idx).t();
      }
    }
    
    // keep
    if (iter > 0 && iter % nThin == 0) {
      keep(save_idx, arma::span(0, p - 1)) = beta.t();
      keep(save_idx, p) = prec;
      if (r > 0) {
        for (int m = 0; m < r; ++m) {
          keep(save_idx, arma::span(p + 1 + m * (n + 3), p + n + m * (n + 3))) = alpha.col(m).t();
          keep(save_idx, arma::span(p + 1 + n + m * (n + 3), p + (m + 1) * (n + 3))) = hp.col(m).t();
        }
      }
      ++save_idx;
    }
  }
  
  if (missing_n == 0) {
    return Rcpp::List::create(
      Rcpp::Named("params") = keep
    );
  } else {
    return Rcpp::List::create(
      Rcpp::Named("params") = keep,
      Rcpp::Named("missing") = keep_Y
    );
  }
}

// [[Rcpp::export]]
double spTMeanRcpp() {

    return 1;
}

// [[Rcpp::export]]
double spTQuantileRcpp() {
  
  return 1;
}

// bayesian kriging
// [[Rcpp::export]]
arma::mat krigeBayesRcpp(
    const arma::mat& w,
    const arma::mat& hp,
    const arma::mat& coords,
    const arma::mat& newcoords,
    const bool parallel = false,
    int nThreads = 0) {
  
  int B  = w.n_rows;
  int n0 = newcoords.n_rows;
  
  arma::mat out(B, n0);
  
  arma::mat d22 = dist_mat(coords, coords);
  arma::mat d11 = dist_mat(newcoords, newcoords);
  arma::mat d21 = dist_mat(coords, newcoords);
  
  if (!parallel) {
    
    int n = coords.n_rows;
    double mu, sigma, decay;
    arma::mat R22(n,n), R11(n0,n0), R21(n,n0); 
    arma::mat C(n0,n0), L(n0,n0), R12R22inv(n0,n); 
    arma::vec wb(n), z(n0);
    
    for (int b = 0; b < B; ++b) {
      mu    = hp(b,0);
      sigma = hp(b,1);
      decay = hp(b,2);
      
      wb = w.row(b).t();
      R22 = arma::exp(-decay * d22);
      R11 = arma::exp(-decay * d11);
      R21 = arma::exp(-decay * d21);
      
      R12R22inv = arma::solve(R22, R21, arma::solve_opts::fast).t();
      C = R11 - R12R22inv * R21;
      L = arma::chol(C, "lower");
      
      z = arma::randn(n0);
      
      out.row(b) = (mu + R12R22inv * (wb - mu) + sigma * (L * z)).t();
    }
    return out;
  }
  
  unsigned long seed = 
    static_cast<unsigned long>(R::runif(0, std::numeric_limits<unsigned long>::max()));
  
#ifdef _OPENMP
  if (nThreads <= 0)
    nThreads = omp_get_max_threads();
  omp_set_num_threads(nThreads);
#endif
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int b = 0; b < B; ++b) {
    
    std::mt19937 rng(seed + 101ULL * b);
    std::normal_distribution<> norm(0.0, 1.0);
    
    double mu    = hp(b,0);
    double sigma = hp(b,1);
    double decay = hp(b,2);
    
    arma::vec wb = w.row(b).t();
    
    arma::mat R22 = arma::exp(-decay * d22);
    arma::mat R11 = arma::exp(-decay * d11);
    arma::mat R21 = arma::exp(-decay * d21);
    
    arma::mat R12R22inv = arma::solve(R22, R21, arma::solve_opts::fast).t();
    arma::mat C = R11 - R12R22inv * R21;
    arma::mat L = arma::chol(C, "lower");
    
    arma::vec z(n0);
    for (int i = 0; i < n0; ++i)
      z(i) = norm(rng);
    
    arma::rowvec outrow =
      (mu + R12R22inv * (wb - mu) + sigma * (L * z)).t();
    
    out.row(b) = outrow;
  }
  
  return out;
}
