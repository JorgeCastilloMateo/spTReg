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
    const std::vector<arma::mat>& X_gamma,
    const std::vector<arma::mat>& X_alpha,
    const arma::mat& dist, 
    const arma::vec& M,    // priors
    const arma::mat& P,
    const arma::vec& prior_sigma,
    const std::vector<arma::vec>& M_beta_gamma,
    const std::vector<arma::mat>& P_beta_gamma,
    const arma::mat& prior_rho_gamma,
    const arma::mat& prior_sigma_gamma,
    const std::vector<arma::vec>& M_beta_alpha,
    const std::vector<arma::mat>& P_beta_alpha,
    const arma::mat& prior_sigma_alpha,
    const arma::mat& prior_phi_alpha,
    const arma::vec& prior_rho_w,
    const arma::vec& prior_sigma_w,
    const arma::vec& prior_phi_w,
    arma::vec beta,        // initial
    double prec,
    arma::mat gamma,    // n x q
    std::vector<arma::vec>& beta_gamma,
    arma::mat hp_gamma, // 2 x q
    arma::mat alpha,    // n x r
    std::vector<arma::vec>& beta_alpha,
    arma::mat hp_alpha, // 2 x r
    arma::vec Wtls,
    arma::vec hp_w,
    const int N,           // constants
    const int n,
    const int T,
    const arma::vec& L,
    const int p,
    const int q,
    const int r,
    const arma::vec& p_gamma,
    const arma::vec& p_alpha,
    const bool wBool,
    const arma::uvec& site,
    const arma::uvec& year,
    const arma::uvec& yday,
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
  std::vector<arma::vec> PM_beta_alpha(r);
  for (int m = 0; m < r; ++m) {
    if (p_alpha(m) > 0)
      PM_beta_alpha[m] = P_beta_alpha[m] * M_beta_alpha[m];
  }
  
  // residual
  arma::vec Xb = X * beta;
  arma::vec e = Y - Xb;
  arma::vec alpha_m(n);
  arma::vec V_m(N);
  arma::uvec site_missing(missing_n);
  for (int m = 0; m < r; ++m) {
    alpha_m = alpha.col(m);
    V_m = V.col(m);
    e -= V_m % alpha_m.elem(site);
  }

  std::vector<arma::vec> Xb_alpha(r);
  for (int m = 0; m < r; ++m) {
    if (p_alpha(m) > 0)
      Xb_alpha[m] = X_alpha[m] * beta_alpha[m];
  }
  
  // aux GP
  std::vector<arma::uvec> site_group(n);
  if (r > 0) {
    for (int i = 0; i < n; ++i) {
      site_group[i] = arma::find(site == i);
    }
  }
  
  int Ndn = N / n; // TL = T * L = N / n
  arma::vec V_block(Ndn), e_block(Ndn);
  
  // aux decay
  int total = 0;
  arma::vec accept(r, arma::fill::zeros);
  arma::vec ratio(r);
  arma::vec sd(r, arma::fill::ones);
  arma::vec lsd(r, arma::fill::zeros); // log(sd);
  
  arma::cube R(n, n, r);
  arma::vec Rlogdet(r);
  std::vector<arma::mat> xR(r);
  
  for (int m = 0; m < r; ++m) {
    R.slice(m) = arma::inv_sympd(exp(- hp_alpha(1, m) * dist));
    Rlogdet(m) = arma::log_det_sympd(R.slice(m));
  }
  
  double decay_aux;
  double ldecay_aux = 0.0;
  arma::vec ldecay = log(hp_alpha.row(1).t());
  arma::mat R_aux(n, n);
  double Rlogdet_aux;
  
  arma::vec vn(n);
  double vtRv, vtRv_aux;
  double ALPHA;
  
  // aux W rho_w sigma_w decay_w
  arma::vec accept_w(2, arma::fill::zeros);
  arma::vec ratio_w(2);
  double sd_w = 1.0;
  double lsd_w = 0.0; // log(sd_w);
  
  arma::mat S_aux(n, n);
  arma::mat S_w(n, n);
  arma::mat R_w(n, n);
  double Rlogdet_w;
  
  if (wBool) {
    S_w = exp(- hp_w(2) * dist);
    R_w = arma::inv_sympd(S_w);
    Rlogdet_w = arma::log_det_sympd(R_w);
  }
  
  double ldecay_w = log(hp_w(2));
  
  arma::vec wn(n);
  arma::mat wtRw(1,1), wtRw_aux(1,1), w1tRw1(1,1);
  
  double rho_aux;
  double onemrho2 = 1 - hp_w(0) * hp_w(0);
  
  double nTd2 = 0.5 * n * T;
  
  int col_idx = 0;
  std::vector<arma::uvec> indW(Ndn);
  if (wBool) {
    for (int t = 0; t < T; ++t) {
      for (int l = 0; l < L(t); ++l) {
        indW[col_idx] = arma::find((year == t) && (yday == l));
        col_idx++;
      }
    }
  }
  
  // full posterior precision rhs of rho
  double Q1, b1;
  // full posterior Precision rhs of beta
  arma::mat Qp(p, p);
  arma::vec bp(p);
  // full posterior Precision rhs of alpha
  arma::mat Qn(n, n);
  arma::vec bn(n);
  // full posterior Precision rhs of beta_alpha
  std::vector<arma::mat> Qp_alpha(r);
  std::vector<arma::vec> bp_alpha(r);
  // full posterior parameters 
  //   A B of prec, C D of prec (alpha), E F of prec (w)
  double A = 0.5 * N + prior_sigma(0);
  double B;
  arma::vec C = 0.5 * n + prior_sigma_alpha.row(0).t();
  double D;
  double E = 0.5 * N + prior_sigma_w(0);
  double F;
  
  // save
  int save_idx = 0;
  int nCols1 = p + 1;
  int nCols2 = r * (n + 2) + arma::accu(p_alpha);
  arma::mat keep(nKeep, nCols1);
  arma::mat keep_sp(nKeep, nCols2);
  arma::mat keep_st(nKeep, N + 3);
  
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
    for (int m = 0; m < r; ++m) {
      alpha_m = alpha.col(m);
      V_m = V.col(m);
      e += V_m % alpha_m.elem(site);
      Qn = hp_alpha(0, m) * R.slice(m);
      bn = Qn * Xb_alpha[m];
      for (int i = 0; i < n; ++i) {
        V_block = V_m.elem(site_group[i]);
        e_block = e.elem(site_group[i]);
        Qn(i,i) += prec * arma::dot(V_block, V_block);
        bn(i)   += prec * arma::dot(V_block, e_block);
      }
      alpha_m = RandomMultiNormalC(Qn, bn);
      //alpha_m -= Vn * onen * (arma::accu(alpha_m) / arma::accu(Vn));
      alpha.col(m) = alpha_m;
      e -= V_m % alpha_m.elem(site);
      
      // mu 
      if (p_alpha(m) > 0) {
        xR[m] = hp_alpha(0, m) * X_alpha[m].t() * R.slice(m);
        Qp_alpha[m] = xR[m] * X_alpha[m] + P_beta_alpha[m];
        bp_alpha[m] = xR[m] * alpha_m + PM_beta_alpha[m];
        beta_alpha[m] = RandomMultiNormalC(Qp_alpha[m], bp_alpha[m]);
        Xb_alpha[m] = X_alpha[m] * beta_alpha[m];
      }
      
      // decay
      ldecay_aux  = R::rnorm(ldecay(m), sd(m));
      decay_aux   = exp(ldecay_aux);
      R_aux       = arma::inv_sympd(exp(- decay_aux * dist));
      Rlogdet_aux = arma::log_det_sympd(R_aux);
      vn       = alpha_m - Xb_alpha[m];
      vtRv_aux = arma::as_scalar(vn.t() * R_aux * vn);
      vtRv     = arma::as_scalar(vn.t() * R.slice(m) * vn);
      ALPHA = 
        0.5 * (Rlogdet_aux - Rlogdet(m) + 
          hp_alpha(0, m) * (vtRv - vtRv_aux)) +
        prior_phi_alpha(0, m) * (ldecay_aux - ldecay(m)) + 
        prior_phi_alpha(1, m) * (hp_alpha(1, m) - decay_aux);
      if (log(R::runif(0, 1)) < ALPHA) {
        ++accept(m);
        hp_alpha(1, m) = decay_aux;
        ldecay(m) = ldecay_aux;
        R.slice(m) = R_aux;
        Rlogdet(m) = Rlogdet_aux;
        vtRv = vtRv_aux;
      }
      
      // prec
      D = 0.5 * vtRv + prior_sigma_alpha(1, m);
      hp_alpha(0, m) = R::rgamma(C(m), 1.0 / D);
    }
    
    // w t=1,...,T l=1,...,L
    if (wBool) {
      e += Wtls;
      col_idx = 0;
      for (int t = 0; t < T; ++t) {
        std::vector<arma::vec> D_list(L(t), arma::vec(n));
        arma::mat et(n, L(t));
        for (int l = 0; l < L(t); ++l) {
          D_list[l].fill(1.0 / prec);
          et.col(l) = e(indW[col_idx]);
          col_idx++;
        }
        et = ffbs(
          et, // n x L_t : observations
          hp_w(0),
          S_w / hp_w(1),
          D_list
        );
        col_idx -= L(t);
        for (int l = 0; l < L(t); ++l) {
          Wtls(indW[col_idx]) = et.col(l);
          col_idx++;
        }
      }
      e -= Wtls;
      
      //// phi_w
      ldecay_aux  = R::rnorm(ldecay_w, sd_w);
      decay_aux   = exp(ldecay_aux);
      S_aux       = exp(- decay_aux * dist);
      R_aux       = arma::inv_sympd(S_aux);
      Rlogdet_aux = arma::log_det_sympd(R_aux);
      
      wtRw_aux.zeros(); 
      wtRw.zeros();
      col_idx = 0;
      for (int t = 0; t < T; ++t) {
        wn = Wtls(indW[col_idx]);
        vn = wn;
        wtRw_aux += onemrho2 * wn.t() * R_aux * wn;
        wtRw     += onemrho2 * wn.t() * R_w * wn;
        col_idx++;
        for (int l = 1; l < L(t); ++l) {
          wn = Wtls(indW[col_idx]) - hp_w(0) * vn; // e_w = ... // w[- l == 1] - rho_w * w[- l == L];
          vn = Wtls(indW[col_idx]);
          wtRw_aux += wn.t() * R_aux * wn;
          wtRw     += wn.t() * R_w * wn;
          col_idx++;
        }
      }
      
      ALPHA = 
        0.5 * (Ndn * (Rlogdet_aux - Rlogdet_w) +
        hp_w(1) * arma::as_scalar(wtRw - wtRw_aux)) + 
        prior_phi_w(0) * (ldecay_aux - ldecay_w) - 
        prior_phi_w(1) * (decay_aux - hp_w(2));
      
      if (log(R::runif(0, 1)) < ALPHA) {
        ++accept_w(1);
        hp_w(2) = decay_aux;
        ldecay_w = ldecay_aux;
        S_w = S_aux;
        R_w = R_aux;
        Rlogdet_w = Rlogdet_aux;
        wtRw = wtRw_aux;
      }
      
      //// prec_w
      F = 0.5 * arma::as_scalar(wtRw) + prior_sigma_w(1);
      hp_w(1) = R::rgamma(E, 1.0 / F);
      
      //// rho_w
      w1tRw1.zeros();
      wtRw.zeros();
      wtRw_aux.zeros();
      col_idx = 0;
      for (int t = 0; t < T; ++t) {
        wn = Wtls(indW[col_idx]);
        w1tRw1 += wn.t() * R_w * wn;
        col_idx++;
        for (int l = 1; l < L(t); ++l) {
          vn = wn;
          wn = Wtls(indW[col_idx]);
          wtRw += vn.t() * R_w * wn;
          wtRw_aux += vn.t() * R_w * vn;
          col_idx++;
        }
      }
      Q1 = hp_w(1) * arma::as_scalar(wtRw_aux) + prior_rho_w(1);
      b1 = hp_w(1) * arma::as_scalar(wtRw) + prior_rho_w(0) * prior_rho_w(1);
      rho_aux = rtnorm(b1 / Q1, 1 / sqrt(Q1), -1.0, 1.0);
      ALPHA = nTd2 * log((1 - rho_aux*rho_aux) / onemrho2) + 
        0.5 * hp_w(1) * (rho_aux*rho_aux - hp_w(0)*hp_w(0)) * arma::as_scalar(w1tRw1);
      
      if (log(R::runif(0, 1)) < ALPHA) {
        ++accept_w(0);
        hp_w(0) = rho_aux;
        onemrho2 = 1 - hp_w(0) * hp_w(0);
      }
    }
    
    // prec
    B = 0.5 * arma::dot(e, e) + prior_sigma(1);
    prec = R::rgamma(A, 1.0 / B);
    
    // tune sd of the proposal for decay's
    if (iter == 0) {
      if (r > 0) accept.zeros();
      if (wBool) accept_w.zeros();
      total = 0;
    } else if ((iter < 1) && (++total % 25 == 0)) {
      if (r > 0) {
        ratio = accept / total;
        for (int m = 0; m < r; ++m) {
          if (ratio(m) > 0.33) {
            lsd(m) += 1 / sqrt(0.04 * total);
          } else {
            lsd(m) -= 1 / sqrt(0.04 * total);
          }
          sd(m) = exp(lsd(m));
        }
      }
      if (wBool) {
        ratio_w = accept_w / total;
        if (ratio_w(1) > 0.33) {
          lsd_w += 1 / sqrt(0.04 * total);
        } else {
          lsd_w -= 1 / sqrt(0.04 * total);
        }
        sd_w = exp(lsd_w);
      }
    }
    
    // Y missing
    if (missing_n > 0) {
      e.elem(missing_idx) -= Y.elem(missing_idx);
      Y.elem(missing_idx) = X.rows(missing_idx) * beta +
        arma::randn(missing_n) / std::sqrt(prec);
      for (int m = 0; m < r; ++m) {
        alpha_m = alpha.col(m);
        V_m = V.col(m);
        site_missing = site.elem(missing_idx);
        Y.elem(missing_idx) += V_m.elem(missing_idx) % alpha_m.elem(site_missing);
      }
      if (wBool)
        Y.elem(missing_idx) += Wtls.elem(missing_idx);
      
      e.elem(missing_idx) += Y.elem(missing_idx);
      
      if (iter > 0 && iter % nThin == 0) {
        keep_Y.row(save_idx) = Y.elem(missing_idx).t();
      }
    }
    
    // keep
    if (iter > 0 && iter % nThin == 0) {
      keep(save_idx, arma::span(0, p - 1)) = beta.t();
      keep(save_idx, p) = prec;
      
      nCols2 = 0;
      for (int m = 0; m < r; ++m) {
        keep_sp(save_idx, arma::span(nCols2, n - 1 + nCols2)) = alpha.col(m).t();
        if (p_alpha(m) > 0)
          keep_sp(save_idx, arma::span(n + nCols2, n + p_alpha(m) - 1 + nCols2)) = beta_alpha[m].t();
        keep_sp(save_idx, arma::span(n + p_alpha(m) + nCols2, n + p_alpha(m) + 1 + nCols2)) = hp_alpha.col(m).t();
        nCols2 += n + p_alpha(m) + 2;
      }
      
      if (wBool) {
        keep_st(save_idx, arma::span(0, N - 1)) = Wtls.t();
        keep_st(save_idx, arma::span(N, N + 2)) = hp_w.t();
      }
      
      ++save_idx;
    }
  }
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("params") = keep);
  
  if (r > 0) out["sp"] = keep_sp;
  if (wBool) out["st"] = keep_st;
  if (missing_n != 0) out["missing"] = keep_Y;
  
  return out;
}

// spatial quantile model
// [[Rcpp::export]]
Rcpp::List spQuantileRcpp(
    const double tau,      // quantile level
    arma::vec Y,           // data
    const arma::mat& X,
    const arma::mat& U,
    const arma::mat& V,
    const std::vector<arma::mat>& X_gamma,
    const std::vector<arma::mat>& X_alpha,
    const arma::mat& dist, 
    const arma::vec& M,    // priors
    const arma::mat& P,
    const arma::vec& prior_sigma,
    const std::vector<arma::vec>& M_beta_gamma,
    const std::vector<arma::mat>& P_beta_gamma,
    const arma::mat& prior_rho_gamma,
    const arma::mat& prior_sigma_gamma,
    const std::vector<arma::vec>& M_beta_alpha,
    const std::vector<arma::mat>& P_beta_alpha,
    const arma::mat& prior_sigma_alpha,
    const arma::mat& prior_phi_alpha,
    const arma::vec& prior_rho_w,
    const arma::vec& prior_sigma_w,
    const arma::vec& prior_phi_w,
    arma::vec beta,        // initial
    double prec,
    arma::mat gamma,    // n x q
    std::vector<arma::vec>& beta_gamma,
    arma::mat hp_gamma, // 2 x q
    arma::mat alpha,    // n x r
    std::vector<arma::vec>& beta_alpha,
    arma::mat hp_alpha, // 2 x r
    arma::vec Wtls,
    arma::vec hp_w,
    const int N,          // constants
    const int n,
    const int T,
    const arma::vec& L,
    const int p,
    const int q,
    const int r,
    const arma::vec& p_gamma,
    const arma::vec& p_alpha,
    const bool wBool,
    const arma::uvec& site,
    const arma::uvec& year,
    const arma::uvec& yday,
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
  std::vector<arma::vec> PM_beta_alpha(r);
  for (int m = 0; m < r; ++m) {
    if (p_alpha(m) > 0)
      PM_beta_alpha[m] = P_beta_alpha[m] * M_beta_alpha[m];
  }
  
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
  arma::uvec site_missing(missing_n);
  for (int m = 0; m < r; ++m) {
    alpha_m = alpha.col(m);
    V_m = V.col(m);
    e -= V_m % alpha_m.elem(site);
  }
  
  std::vector<arma::vec> Xb_alpha(r);
  for (int m = 0; m < r; ++m) {
    if (p_alpha(m) > 0)
      Xb_alpha[m] = X_alpha[m] * beta_alpha[m];
  }
  
  // aux GP
  std::vector<arma::uvec> site_group(n);
  if (r > 0) {
    for (int i = 0; i < n; ++i) {
      site_group[i] = arma::find(site == i);
    }
  }
  
  int Ndn = N / n; // TL = T * L = N / n
  arma::vec V_block(Ndn), e_block(Ndn);
  arma::vec c2dxi_block(Ndn), V_c2dxi_block(Ndn);
  
  // aux decay
  int total = 0;
  arma::vec accept(r, arma::fill::zeros);
  arma::vec ratio(r);
  arma::vec sd(r, arma::fill::ones);
  arma::vec lsd(r, arma::fill::zeros); // log(sd);
  
  arma::cube R(n, n, r);
  arma::vec Rlogdet(r);
  std::vector<arma::mat> xR(r);
  
  for (int m = 0; m < r; ++m) {
    R.slice(m) = arma::inv_sympd(exp(- hp_alpha(1, m) * dist));
    Rlogdet(m) = arma::log_det_sympd(R.slice(m));
  }
  
  double decay_aux;
  double ldecay_aux = 0.0;
  arma::vec ldecay = log(hp_alpha.row(1).t());
  arma::mat R_aux(n, n);
  double Rlogdet_aux;
  
  arma::vec vn(n);
  double vtRv, vtRv_aux;
  double ALPHA;
  
  // aux W rho_w sigma_w decay_w
  arma::vec accept_w(2, arma::fill::zeros);
  arma::vec ratio_w(2);
  double sd_w = 1.0;
  double lsd_w = 0.0; // log(sd_w);

  arma::mat S_aux(n, n);
  arma::mat S_w(n, n);
  arma::mat R_w(n, n);
  double Rlogdet_w;
  
  if (wBool) {
    S_w = exp(- hp_w(2) * dist);
    R_w = arma::inv_sympd(S_w);
    Rlogdet_w = arma::log_det_sympd(R_w);
  }

  double ldecay_w = log(hp_w(2));

  arma::vec wn(n);
  arma::mat wtRw(1,1), wtRw_aux(1,1), w1tRw1(1,1);

  double rho_aux;
  double onemrho2 = 1 - hp_w(0) * hp_w(0);

  double nTd2 = 0.5 * n * T;
  
  int col_idx = 0;
  std::vector<arma::uvec> indW(Ndn);
  if (wBool) {
    for (int t = 0; t < T; ++t) {
      for (int l = 0; l < L(t); ++l) {
        indW[col_idx] = arma::find((year == t) && (yday == l));
        col_idx++;
      }
    }
  }
  
  // full posterior precision rhs of rho
  double Q1, b1;
  // full posterior Precision rhs of beta
  arma::mat Qp(p, p);
  arma::vec bp(p);
  // full posterior Precision rhs of alpha
  arma::mat Qn(n, n);
  arma::vec bn(n);
  // full posterior Precision rhs of beta_alpha
  std::vector<arma::mat> Qp_alpha(r);
  std::vector<arma::vec> bp_alpha(r);
  // full posterior parameters 
  //   A B of prec, C D of prec (alpha), E F of prec (w)
  double A = 1.5 * N + prior_sigma(0);
  double B;
  arma::vec C = 0.5 * n + prior_sigma_alpha.row(0).t();
  double D;
  double E = 0.5 * N + prior_sigma_w(0);
  double F;
  
  // save
  int save_idx = 0;
  int nCols1 = p + 1;
  int nCols2 = r * (n + 2) + arma::accu(p_alpha);
  arma::mat keep(nKeep, nCols1);
  arma::mat keep_sp(nKeep, nCols2);
  arma::mat keep_st(nKeep, N + 3);
  
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
    for (int m = 0; m < r; ++m) {
      alpha_m = alpha.col(m);
      V_m = V.col(m);
      e += V_m % alpha_m.elem(site);
      Qn = hp_alpha(0, m) * R.slice(m);
      bn = Qn * Xb_alpha[m];
      for (int i = 0; i < n; ++i) {
        V_block = V_m.elem(site_group[i]);
        e_block = e.elem(site_group[i]);
        c2dxi_block = c2dxi.elem(site_group[i]);
        V_c2dxi_block = V_block % c2dxi_block;
        Qn(i,i) += prec * arma::accu(V_c2dxi_block % V_block);
        bn(i)   += prec * arma::accu(V_c2dxi_block % e_block);
      }
      alpha_m = RandomMultiNormalC(Qn, bn);
      //alpha_m -= Vn * onen * (arma::accu(alpha_m) / arma::accu(Vn));
      alpha.col(m) = alpha_m;
      e -= V_m % alpha_m.elem(site);

      // mu 
      if (p_alpha(m) > 0) {
        xR[m] = hp_alpha(0, m) * X_alpha[m].t() * R.slice(m);
        Qp_alpha[m] = xR[m] * X_alpha[m] + P_beta_alpha[m];
        bp_alpha[m] = xR[m] * alpha_m + PM_beta_alpha[m];
        beta_alpha[m] = RandomMultiNormalC(Qp_alpha[m], bp_alpha[m]);
        Xb_alpha[m] = X_alpha[m] * beta_alpha[m];
      }
      
      // decay
      ldecay_aux  = R::rnorm(ldecay(m), sd(m));
      decay_aux   = exp(ldecay_aux);
      R_aux       = arma::inv_sympd(exp(- decay_aux * dist));
      Rlogdet_aux = arma::log_det_sympd(R_aux);
      vn       = alpha_m - Xb_alpha[m];
      vtRv_aux = arma::as_scalar(vn.t() * R_aux * vn);
      vtRv     = arma::as_scalar(vn.t() * R.slice(m) * vn);
      ALPHA = 
        0.5 * (Rlogdet_aux - Rlogdet(m) + 
          hp_alpha(0, m) * (vtRv - vtRv_aux)) +
        prior_phi_alpha(0, m) * (ldecay_aux - ldecay(m)) + 
        prior_phi_alpha(1, m) * (hp_alpha(1, m) - decay_aux);
      if (log(R::runif(0, 1)) < ALPHA) {
        ++accept(m);
        hp_alpha(1, m) = decay_aux;
        ldecay(m) = ldecay_aux;
        R.slice(m) = R_aux;
        Rlogdet(m) = Rlogdet_aux;
        vtRv = vtRv_aux;
      }
      
      // prec
      D = 0.5 * vtRv + prior_sigma_alpha(1, m); ;
      hp_alpha(0, m) = R::rgamma(C(m), 1.0 / D);
    }
    
    // w t=1,...,T l=1,...,L
    if (wBool) {
      e += Wtls;
      col_idx = 0;
      for (int t = 0; t < T; ++t) {
        std::vector<arma::vec> D_list(L(t), arma::vec(n));
        arma::mat et(n, L(t));
        for (int l = 0; l < L(t); ++l) {
          D_list[l] = 1.0 / (c2dxi(indW[col_idx]) * prec);
          et.col(l) = e(indW[col_idx]);
          col_idx++;
        }
        et = ffbs(
          et, // n x L_t : observations
          hp_w(0),
          S_w / hp_w(1),
          D_list
        );
        col_idx -= L(t);
        for (int l = 0; l < L(t); ++l) {
          Wtls(indW[col_idx]) = et.col(l);
          col_idx++;
        }
      }
      e -= Wtls;
      
      //// phi_w
      ldecay_aux  = R::rnorm(ldecay_w, sd_w);
      decay_aux   = exp(ldecay_aux);
      S_aux       = exp(- decay_aux * dist);
      R_aux       = arma::inv_sympd(S_aux);
      Rlogdet_aux = arma::log_det_sympd(R_aux);
      
      wtRw_aux.zeros(); 
      wtRw.zeros();
      col_idx = 0;
      for (int t = 0; t < T; ++t) {
        wn = Wtls(indW[col_idx]);
        vn = wn;
        wtRw_aux += onemrho2 * wn.t() * R_aux * wn;
        wtRw     += onemrho2 * wn.t() * R_w * wn;
        col_idx++;
        for (int l = 1; l < L(t); ++l) {
          wn = Wtls(indW[col_idx]) - hp_w(0) * vn; // e_w = ... // w[- l == 1] - rho_w * w[- l == L];
          vn = Wtls(indW[col_idx]);
          wtRw_aux += wn.t() * R_aux * wn;
          wtRw     += wn.t() * R_w * wn;
          col_idx++;
        }
      }
      
      ALPHA = 
        0.5 * (Ndn * (Rlogdet_aux - Rlogdet_w) +
        hp_w(1) * arma::as_scalar(wtRw - wtRw_aux)) + 
        prior_phi_w(0) * (ldecay_aux - ldecay_w) - 
        prior_phi_w(1) * (decay_aux - hp_w(2));
      
      if (log(R::runif(0, 1)) < ALPHA) {
        ++accept_w(1);
        hp_w(2) = decay_aux;
        ldecay_w = ldecay_aux;
        S_w = S_aux;
        R_w = R_aux;
        Rlogdet_w = Rlogdet_aux;
        wtRw = wtRw_aux;
      }
      
      //// prec_w
      F = 0.5 * arma::as_scalar(wtRw) + prior_sigma_w(1);
      hp_w(1) = R::rgamma(E, 1.0 / F);
      
      //// rho_w
      w1tRw1.zeros();
      wtRw.zeros();
      wtRw_aux.zeros();
      col_idx = 0;
      for (int t = 0; t < T; ++t) {
        wn = Wtls(indW[col_idx]);
        w1tRw1 += wn.t() * R_w * wn;
        col_idx++;
        for (int l = 1; l < L(t); ++l) {
          vn = wn;
          wn = Wtls(indW[col_idx]);
          wtRw += vn.t() * R_w * wn;
          wtRw_aux += vn.t() * R_w * vn;
          col_idx++;
        }
      }
      Q1 = hp_w(1) * arma::as_scalar(wtRw_aux) + prior_rho_w(1);
      b1 = hp_w(1) * arma::as_scalar(wtRw) + prior_rho_w(0) * prior_rho_w(1);
      rho_aux = rtnorm(b1 / Q1, 1 / sqrt(Q1), -1.0, 1.0);
      ALPHA = nTd2 * log((1 - rho_aux*rho_aux) / onemrho2) + 
        0.5 * hp_w(1) * (rho_aux*rho_aux - hp_w(0)*hp_w(0)) * arma::as_scalar(w1tRw1);
      
      if (log(R::runif(0, 1)) < ALPHA) {
        ++accept_w(0);
        hp_w(0) = rho_aux;
        onemrho2 = 1 - hp_w(0) * hp_w(0);
      }
    }
    
    // prec
    B = arma::accu(xi + 0.5 * c2dxi % arma::square(e)) + prior_sigma(1);
    prec = R::rgamma(A, 1.0 / B);
    
    // tune sd of the proposal for decays
    if (iter == 0) {
      if (r > 0) accept.zeros();
      if (wBool) accept_w.zeros();
      total = 0;
    } else if ((iter < 1) && (++total % 25 == 0)) {
      if (r > 0) {
        ratio = accept / total;
        for (int m = 0; m < r; ++m) {
          if (ratio(m) > 0.33) {
            lsd(m) += 1 / sqrt(0.04 * total);
          } else {
            lsd(m) -= 1 / sqrt(0.04 * total);
          }
          sd(m) = exp(lsd(m));
        }
      }
      if (wBool) {
        ratio_w = accept_w / total;
        if (ratio_w(1) > 0.33) {
          lsd_w += 1 / sqrt(0.04 * total);
        } else {
          lsd_w -= 1 / sqrt(0.04 * total);
        }
        sd_w = exp(lsd_w);
      }
    }
    
    // Y missing
    if (missing_n > 0) {
      e.elem(missing_idx) -= Y.elem(missing_idx);
      mm = X.rows(missing_idx) * beta + c1 * xi.elem(missing_idx);
      ss = arma::sqrt(xi.elem(missing_idx) / (prec * c2));
      Y.elem(missing_idx) = mm + ss % arma::randn(missing_n);
      for (int m = 0; m < r; ++m) {
        alpha_m = alpha.col(m);
        V_m = V.col(m);
        site_missing = site.elem(missing_idx);
        Y.elem(missing_idx) += V_m.elem(missing_idx) % alpha_m.elem(site_missing);
      }
      if (wBool)
        Y.elem(missing_idx) += Wtls.elem(missing_idx);
      
      e.elem(missing_idx) += Y.elem(missing_idx);
      
      if (iter > 0 && iter % nThin == 0) {
        keep_Y.row(save_idx) = Y.elem(missing_idx).t();
      }
    }
    
    // keep
    if (iter > 0 && iter % nThin == 0) {
      keep(save_idx, arma::span(0, p - 1)) = beta.t();
      keep(save_idx, p) = prec;
      
      nCols2 = 0;
      for (int m = 0; m < r; ++m) {
        keep_sp(save_idx, arma::span(nCols2, n - 1 + nCols2)) = alpha.col(m).t();
        if (p_alpha(m) > 0)
          keep_sp(save_idx, arma::span(n + nCols2, n + p_alpha(m) - 1 + nCols2)) = beta_alpha[m].t();
        keep_sp(save_idx, arma::span(n + p_alpha(m) + nCols2, n + p_alpha(m) + 1 + nCols2)) = hp_alpha.col(m).t();
        nCols2 += n + p_alpha(m) + 2;
      }
      
      if (wBool) {
        keep_st(save_idx, arma::span(0, N - 1)) = Wtls.t();
        keep_st(save_idx, arma::span(N, N + 2)) = hp_w.t();
      }
      
      ++save_idx;
    }
  }
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("params") = keep);

  if (r > 0) out["sp"] = keep_sp;
  if (wBool) out["st"] = keep_st;
  if (missing_n != 0) out["missing"] = keep_Y;

  return out;
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
