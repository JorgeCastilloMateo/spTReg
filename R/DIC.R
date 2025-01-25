#' @title Deviance information criterion
#' @aliases DIC.default DIC
#' @description Generic function calculating Deviance information criterion.
#'   
#' @details The DIC is calculated according to the formula 
#'   \deqn{\text{DIC} = \overline{D(\bm{\theta})} + p_{D},}
#'   where \eqn{D(\bm{\theta})} is the deviance, and 
#'   \eqn{p_{D} = \overline{D(\bm{\theta})} - D(\bar{\bm{\theta})}} is the 
#'   effective number of parameters. Bars represent the posterior mean.
#' 
#'   When comparing hierarchical models fitted by MCMC to the same 
#'   data, the smaller the DIC, the better the fit.
#'   
#' @param object Object of class \code{"iidm"} or \code{"spTm"}.
#' @param ... currently no additional arguments.
#'  
#' @return A numeric vector with three components corresponding to the DIC, the
#'   mean deviance \eqn{\overline{D(\bm{\theta})}}, and the effective number
#'   of parameters \eqn{p_{D}}.
#'   
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{deviance}}, \code{\link{iidm}}, \code{\link{spTm}}
#' @references 
#' Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A (2002). 
#' Bayesian measures of model complexity and fit (with discussion). 
#' \emph{Journal of the Royal Statistical Society, Series B}, \strong{64}(4), 583--639. 
#' \doi{10.1111/1467-9868.00353}
#' 
#' @export DIC
DIC <- function(object, ...) {
  UseMethod("DIC")
}

#' @rdname DIC
#' @method DIC default
#' @export DIC.default
#' @export 
DIC.default <- function(
    object,
    ...) {
  
  D  <- mean(deviance(object, value = "samples"))
  pD <- D - deviance(object, value = "mean")
    
  dic <- c("DIC" = D + pD, "deviance" = D, "pD" = pD)
  
  return(dic)
}
