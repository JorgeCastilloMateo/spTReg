#' @title Spatial and temporal predictions for the space-time models
#' @aliases predict.iidm predict.spTm
#' 
#' @description This function is used to obtain spatial predictions and also to
#'   get the temporal forecasts using MCMC samples.
#'   
#' @details
#'   Function currently working for iidm, not working for spTm.
#'   
#'   If \code{type = "signal"}, the function computes the product 
#'   \eqn{\textbf{X}_{\text{new}} \bm{\beta}}.
#'   
#'   If \code{type = "response"}, the function computes the previous product 
#'   and adds a Gaussian or Laplacian error as appropriate. 
#' 
#' @param object Object of class \code{"iidm"} or \code{"spTm"}.
#' @param newdata An optional data frame in which to look for variables with 
#'   which to predict. If omitted, the fitted values are used.
#' 
#'   An optional data set providing the explanatory variables for
#'   spatial prediction or temporal forecasts. This data should have the same 
#'   space-time structure as the original \code{data}, see \code{\link{spTm}}.
#'   The index columns \code{t}, \code{l}, and \code{i} must follow the 
#'   numbering started at \code{data}, e.g., the first new location must be 
#'   $n+1$. If omitted, the fitted values are used.
#' @param newcoords The coordinates for the prediction or forecast sites. The 
#'   locations are in similar format to \code{coords}, see \code{\link{spTm}}.
#'   If omitted, the fitted values are used.
#' @param type Type of prediction (signal or response). Can be abbreviated.
#' @param ... currently no additional arguments.
#' 
#' @return A \code{coda} object holds the response variable posterior predictive 
#'   samples (or posterior quantile samples). The rows of this matrix 
#'   correspond to the posterior predictive samples and the columns are the 
#'   predicted values. 
#'   
#' @author Jorge Castillo-Mateo
#' 
#' @seealso \code{\link{iidm}}, \code{\link{spTm}}
#' 
#' @rdname predict
#' @method predict iidm
#' @export predict.iidm
#' @export 
predict.iidm <- function(
    object, 
    newdata, 
    type = c("signal", "response"), 
    ...) {
  
  tt <- terms(object)
  
  type <- match.arg(type)
  
  noData <- (missing(newdata) || is.null(newdata))
  if (noData) {
    class(object) <- "lm"
    mm <- X <- model.matrix(object)
    mmDone <- TRUE
    offset <- object$offset
  } else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    
    offset <- model.offset(m)
    if (!is.null(addO <- object$call$offset)) {
      addO <- eval(addO, newdata, environment(tt))
      offset <- if (length(offset)) 
        offset + addO
      else addO
    }
    
  }
  
  n <- nrow(X)
  p <- object$rank
  p1 <- seq_len(p)
  
  beta <- object$p.params.samples[,p1]
  
  predictor <- t(tcrossprod(X, beta))
  
  if (!is.null(offset)) 
    predictor <- predictor + offset
  
  if (type == "response") {
    sigma <- object$p.params.samples[, p + 1]
    B <- length(sigma)
    if (object$method == "mean") {
      predictor <- predictor + stats::rnorm(B * n, sd = sigma)
    } else if (object$method == "quantile") {
      predictor <- predictor + c(ralRcpp(rep.int(sigma, n), object$quantile))
    } else {
      stop("'method' in 'object' must be 'mean' or 'quantile'")
    }
  }
  
  predictor <- coda::mcmc(predictor, 
    start = object$mcmc$n.burnin + 1, 
    end = object$mcmc$n.burnin + object$mcmc$n.samples, 
    thin = object$mcmc$n.thin)
  
  return(predictor)
}

#' @rdname predict
#' @method predict spTm
#' @export predict.spTm
#' @export 
predict.spTm <- function(
    object, 
    newdata = NULL, 
    newcoords = NULL,
    type = c("signal", "response"), 
    ...) {
  
  return(1)
}
