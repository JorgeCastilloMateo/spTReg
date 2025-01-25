#' @title Deviance
#' @aliases deviance.iidm deviance.spTm
#' @description Returns the deviance of a fitted model object \code{"iidm"} or
#'   \code{"spTm"}.
#'   
#' @details The deviance is defined according to
#'   \deqn{D(\bm{\theta}) = -2 \log([\mathbf{y} \mid \bm{\theta}]),}
#'   where \eqn{\mathbf{y}} are the data, \eqn{\bm{\theta}} are the unknown
#'   parameters of the model, and \eqn{[\mathbf{y} \mid \bm{\theta}]} is the 
#'   likelihood function.
#'   
#'   The deviance is only defined up to an additive constant. Different 
#'   constants have conventionally been used for different
#'   purposes and so \code{deviance} may give different values.
#'   
#' @param object Object of class \code{"iidm"} or \code{"spTm"}.
#' @param value Model parameters used to compute the deviance (mean 
#'   \eqn{D(\bar{\bm{\theta}})} or samples \eqn{D(\bm{\theta}^{(b)})}). Can be 
#'   abbreviated.
#' @param ... currently no additional arguments.
#'  
#' @return The value of the deviance for each MCMC sample extracted from the 
#'   object \code{object}, or for the posterior mean of the MCMC samples 
#'   according to \code{value}.
#'   
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{DIC}}, \code{\link{iidm}}, \code{\link{spTm}}
#' @references 
#' Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A (2002). 
#' Bayesian measures of model complexity and fit (with discussion). 
#' \emph{Journal of the Royal Statistical Society, Series B}, \strong{64}(4), 583--639. 
#' \doi{10.1111/1467-9868.00353}
#' 
#' @rdname deviance
#' @method deviance iidm
#' @export deviance.iidm
#' @export 
deviance.iidm <- function(
    object,
    value = c("mean", "samples"),
    ...) {
  
  value <- match.arg(value)
  
  y <- model.response(model.frame(object))
  pred <- predict(object, type = "signal")
  
  if (object$method == "mean") {
    if (value == "mean") {
      pred  <- colMeans(pred)
      sigma <- mean(object$p.params.samples[,"sigma"])
      dev <- -2 * sum(stats::dnorm(y, mean = pred, sd = sigma, log = TRUE))
    } else {
      B <- nrow(object$p.params.samples)
      dev <- rep(NA, B)
      for (b in 1:B) {
        sigma <- object$p.params.samples[b, "sigma"]
        dev[b] <- -2 * sum(stats::dnorm(y, mean = pred[b, ], sd = sigma, log = TRUE))
      }
    }
  } else if (object$method == "quantile") {
    if (value == "mean") {
      pred  <- colMeans(pred)
      sigma <- mean(object$p.params.samples[,"sigma"])
      dev <- -2 * sum(dal(y, pred, sigma, object$quantile, log = TRUE))
    } else {
      B <- nrow(object$p.params.samples)
      dev <- rep(NA, B)
      for (b in 1:B) {
        sigma <- object$p.params.samples[b, "sigma"]
        dev[b] <- -2 * sum(dal(y, pred[b, ], sigma, object$quantile, log = TRUE))
      }
    }
  } else {
    stop("'method' in 'object' must be 'mean' or 'quantile'")
  }
  
  return(dev)
}

#' @rdname deviance
#' @method deviance spTm
#' @export deviance.spTm
#' @export 
deviance.spTm <- function(
    object,
    value = c("mean", "samples"),
    ...) {
  
  return(1)
}
