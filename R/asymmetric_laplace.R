#' @name Asymmetric-Laplace
#' @title The Asymmetric Laplace Distribution
#' @aliases dal qal pal ral
#' 
#' @description Density, distribution function, quantile function and random
#'   generation for the asymmetric Laplace distribution with location equal 
#'   to \code{mu}, scale equal to \code{sigma}, and asymmetry equal to 
#'   \code{tau}.
#'   
#' @details The asymmetric Laplace distribution with ... has density
#'   \deqn{f(x) = ...}
#'   for \eqn{x ...}.
#'   
#'   \eqn{p(x)} is computed using ... 
#'   
#'   The quantile is defined as ...
#'   
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param size The Poisson binomial distribution has \code{size} times the
#'   vector of probabilities \code{prob}.
#' @param prob Vector with the probabilities of success on each trial.
#' @param log,log.p Logical. If \code{TRUE}, probabilities \eqn{p} are given as 
#'   \eqn{\log(p)}.
#' @param lower.tail Logical. If \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le x)}, otherwise, \eqn{P(X > x)}. 
#' @return 
#'   \code{dal} gives the density, \code{pal} gives the distribution function,
#'   \code{qal} gives the quantile function and \code{ral} generates random 
#'   deviates.
#'   
#'   The length of the result is determined by ...
#'   
#' @author Jorge Castillo-Mateo 
#' @references 
#' Kozumi H, Kobayashi G (2011). 
#' Gibbs sampling methods for Bayesian quantile regression. 
#' \emph{Journal of Statistical Computation and Simulation}, \strong{81}(11), 1565--1578. 
#' \doi{10.1080/00949655.2010.496117}
#' 
#' Yu K, Moyeed RA (2001). 
#' Bayesian quantile regression. 
#' \emph{Statistics & Probability Letters}, \strong{54}(4), 437--447.
#' \doi{10.1016/S0167-7152(01)00124-9}
#' 
#' @rdname dal
#' @export dal
dal <- function(x, mu = 0, sigma = 1, tau = 0.5, log = FALSE) {
  if (log) {
    return(log(tau * (1 - tau) / sigma) - .delta((x - mu) / sigma, tau))
  } else {
    return(tau * (1 - tau) / sigma * exp(-.delta((x - mu) / sigma, tau)))
  }
}

.delta <- function(u, tau) {
  u * (tau - (u < 0))
  #if (u < 0) {
  #  return(u * tau - u)
  #} else {
  #  return(u * tau)
  #}
}

#' @rdname dal
#' @export pal
pal <- function(q, mu = 0, sigma = 1, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  return(1)
}

#' @rdname dal
#' @export qal
qal <- function(p, mu = 0, sigma = 1, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  
  if (log.p) { p <- exp(p) }
  
  if (tau < 0 | tau > 1) { stop("'tau' should take values in [0,1]" )}
  
  return(1)
}

#' @rdname dal
#' @export ral
ral <- function(n, mu = 0, sigma = 1, tau = 0.5) {
  return(1)
}
