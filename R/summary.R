#' @title Summarizing model fits
#' @aliases summary.iidm
#' 
#' @description \code{summary} method for class \code{"iidm"}.
#' 
#' @description
#'   Function currently not working.
#'  
#' @param object an object of class \code{"iidm"}, usually, a result of a call 
#'   to \code{\link{iidm}}. 
#' @param ... further arguments passed to or from other methods.
#' @return The function \code{summary.iidm} computes and returns a list of 
#'   summary statistics of the fitted i.i.d. linear model given in 
#'   \code{object}, using the components (list elements) \code{"call"} and 
#'   \code{"terms"} from its argument, plus
#'   \item{r.squared}{if method="mean". A \code{coda} object of Bayesian \eqn{R^{2}} 
#'     (Gelman et al., 2018). The residual based \eqn{R^{2}} uses draws from 
#'     the residual distribution and is defined as 
#'     \deqn{R^{2,b} = \frac{Var_{\mu}^{b}}{Var_{\mu}^{b} + Var_{res}^{b}},}
#'     where
#'     \deqn{Var_{\mu}^{b} = V_{i=1}^{n} \hat{y}_{i}^{b}, \quad 
#'     Var_{res}^{b} = V_{i=1}^{n} \hat{e}_{i}^{b}.}
#'     The model based \eqn{R^{2}} uses draws from the modeled residual 
#'     variances and is defined as
#'     \deqn{Var_{res}^{b} = \sigma^{2,b}.}
#'   }
#'   \item{r.one}{if method="quantile". A \code{coda} object of Bayesian 
#'     \eqn{R^{1}(\tau)}. Instead of the residual sum of squares it uses the 
#'     residual weighted sum of absolute values associated with quantiles.}
#'   
#' @author Jorge Castillo-Mateo
#' 
#' @references 
#' Gelman A, Goodrich B, Gabry J, Vehtari A (2018). 
#' R-squared for Bayesian regression models.
#' \emph{The American Statistician}, \strong{73}(3), 307--309. 
#' \doi{10.1080/00031305.2018.1549100}.
#' 
#' @rdname summary
#' @method summary iidm
#' @export summary.iidm
#' @export 
summary.iidm <- function(object, ...) {
  
  z <- list()
  
  if (object$method == "mean") {
    v.mu <- apply(object$fitted.values, 1, stats::var)
    v.re <- apply(object$residuals, 1, stats::var)
    
    z$r.squared <- coda::mcmc(
      cbind(
        "residual" = v.mu / (v.mu + v.re),
        "model"    = v.mu / (v.mu + object$p.params.samples[, "sigma"]^2)
      ),
      start = attr(object$p.params.samples, "mcpar")[1],
      end = attr(object$p.params.samples, "mcpar")[2],
      thin = attr(object$p.params.samples, "mcpar")[3]
    )
  } else {
    v.mu <- apply(object$fitted.values, 1, .varq, quantile = object$quantile)
    v.re <- apply(object$residuals, 1, .varq, quantile = object$quantile)
    
    z$r.one <- coda::mcmc(
      cbind(
        "residual" = v.mu / (v.mu + v.re),
        "model"    = v.mu / (v.mu + object$p.params.samples[, "sigma"])
      ),
      start = attr(object$p.params.samples, "mcpar")[1],
      end = attr(object$p.params.samples, "mcpar")[2],
      thin = attr(object$p.params.samples, "mcpar")[3]
    )
  }
  
  class(z) <- "summary.iidm"
  
  return(z)
}

.varq <- function(x, quantile = 0.5) {
  
  q <- stats::quantile(x, prob = quantile)
  e <- x - q
  p <- x < q
  qhat <- (quantile - 1) * sum(e[p]) + quantile * sum(e[!p])
  
  return(qhat / (length(x) - 1))
}
