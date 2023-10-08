#' @title Credible intervals for model parameters
#' @aliases confint.iidm
#' 
#' @description Computes credible intervals for one or more parameters in a 
#'   fitted model. Package \strong{spTReg} adds methods for \code{iidm} and ...
#'   fits.
#'   
#' @details 
#'   \code{\link[stats]{confint}} is a generic function in package stats.
#'   
#'   These \code{confint} methods call the \code{\link[stats]{quantile}} 
#'   function, then find the credible intervals through the empirical quantiles
#'   of the posterior samples.
#'  
#' @param object a fitted model object. Methods currently exist for the classes 
#'   \code{"iidm"}, ... and ... .
#' @param parm a specification of which parameters are to be given credible 
#'   intervals, either a vector of numbers or a vector of names. If missing, 
#'   all parameters are considered.
#' @param level the credible level required.
#' @param ... additional argument(s) for methods.
#' @return A matrix (or vector) with columns giving the mean, lower and upper 
#'   credible limits for each parameter. These will be labelled as Estimate, 
#'   (1 - level)/2 and 1 - (1 - level)/2 in \% (by default 2.5\% and 97.5\%).
#'   
#' @author Jorge Castillo-Mateo
#' 
#' @rdname confint
#' @method confint iidm
#' @export confint.iidm
#' @export 
confint.iidm <- function(object, parm, level = 0.95, ...) {
  
  if (missing(parm)) {
    params <- object$p.params.samples
  } else {
    params <- object$p.params.samples[, parm]
  }
  
  alpha <- (1 - level) / 2
  
  m  <- colMeans(params)
  ci <- apply(params, 2, stats::quantile, prob = c(alpha, 1 - alpha))
  
  return(t(rbind("Estimate" = m, ci)))
}
