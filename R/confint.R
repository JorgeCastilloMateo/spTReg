#' @title Credible intervals for model parameters
#' @aliases confint.iidm confint.spTm
#' 
#' @description Computes credible intervals for one or more parameters in a 
#'   fitted model. Package \strong{spTReg} adds methods for \code{iidm} and 
#'   \code{spTm} fits.
#'   
#' @details 
#'   \code{\link[stats]{confint}} is a generic function in package stats.
#'   
#'   These \code{confint} methods call the \code{\link[stats]{quantile}} 
#'   function, then find the credible intervals through the empirical quantiles
#'   of the posterior samples.
#'  
#' @param object a fitted model object. Methods currently exist for the classes 
#'   \code{"iidm"} and \code{"spTm"}.
#' @param param a specification of which parameters are to be given credible 
#'   intervals, either a vector of numbers or a vector of names. If missing, 
#'   all parameters are considered.
#' @param level the credible level required.
#' @param ... currently no additional arguments.
#' @return A matrix (or vector) with columns giving the mean, lower and upper 
#'   credible limits for each parameter. These will be labelled as Estimate, 
#'   \eqn{(1 - \code{level})/2} and \eqn{1 - (1 - \code{level})/2} in \% (by 
#'   default 2.5\% and 97.5\%).
#'   
#' @author Jorge Castillo-Mateo
#' 
#' @seealso \code{\link{spTm}}, \code{\link{iidm}}
#' 
#' @rdname confint
#' @method confint iidm
#' @export confint.iidm
#' @export 
confint.iidm <- function(object, param, level = 0.95, ...) {
  
  if (missing(param)) {
    params <- object$p.params.samples
  } else {
    params <- object$p.params.samples[, param, drop = FALSE]
  }
  
  alpha <- (1 - level) / 2
  
  m  <- colMeans(params)
  ci <- apply(params, 2, stats::quantile, prob = c(alpha, 1 - alpha))
  
  return(t(rbind("Estimate" = m, ci)))
}

#' @rdname confint
#' @method confint spTm
#' @export confint.spTm
#' @export 
confint.spTm <- function(object, param, level = 0.95, ...) {
  
  if (missing(param)) {
    params <- object$p.params.samples
  } else {
    params <- object$p.params.samples[, param, drop = FALSE]
  }
  
  alpha <- (1 - level) / 2
  
  m  <- colMeans(params)
  ci <- apply(params, 2, stats::quantile, prob = c(alpha, 1 - alpha))
  
  return(t(rbind("Estimate" = m, ci)))
}