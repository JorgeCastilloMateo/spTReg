#' @title MCMC sampling for the i.i.d. linear models
#' @description This function is used to draw MCMC samples using the Gibbs 
#'   sampler that fits the Bayesian linear regression models with i.i.d. 
#'   errors.
#'   
#' @details
#'   The function \code{iidm} can fit three types of model, depending on whether
#'   the mean or the quantile is of interest, or binary data is of interest.
#'   
#'   General model:
#'   \deqn{Y_{i} = \mathbf{x}_{i} \bm{\beta} + \epsilon_{i}, \qquad i=1,\ldots,n,}
#'   with \eqn{\epsilon_i \sim \text{i.i.d. } N(0, \sigma^{2})} for the mean, or 
#'   \eqn{\epsilon_i \sim \text{i.i.d. } AL(0, \sigma, \tau)} for the 
#'   \eqn{\tau}-th quantile.
#'   
#'   Priors:
#'   \deqn{\bm{\beta} \sim N_{k}(\bm{\mu}_{\bm{\beta}}, \bm{\Sigma}_{\bm{\beta}}^{-1})}
#'   \deqn{\sigma^{2} \sim IG(a_{\sigma}, b_{\sigma})}
#'   
#'   Models for binary data are not currently available.
#'   
#'   
#' @param formula an object of class \code{"\link[stats]{formula}"} (or one 
#'   that can be coerced to that class): a symbolic description of the model to
#'   be fitted. The details of model specification are given under ‘Details’.
#' @param data an optional data frame, list or environment (or object coercible 
#'   by \code{\link{as.data.frame}} to a data frame) containing the variables 
#'   in the model. If not found in data, the variables are taken from 
#'   \code{environment(formula)}, typically the environment from which 
#'   \code{iidm} is called.
#' @param subset an optional vector specifying a subset of observations to be 
#'   used in the fitting process. (See additional details about how this 
#'   argument interacts with data-dependent bases in the `Details' section of 
#'   the \code{\link[stats]{model.frame}} documentation.)
#' @param method the method to be used for fitting, \code{"mean"} or 
#'   \code{"quantile"} regression (default=\code{"mean"}).
#' @param quantile a numeric scalar containing the quantile level of interest 
#'   (default=0.5).
#' @param priors a list with each tag corresponding to a parameter name. Valid 
#'   tags are \code{"beta"} and \code{"sigma"}. The \code{"beta"} tag can be 
#'   either the vector \eqn{(a_{\bm{\beta}}, b_{\bm{\beta}})} such that 
#'   \eqn{\bm{\mu}_{\bm{\beta}} = a_{\bm{\beta}} \mathbf{1}_{k}} and 
#'   \eqn{\bm{\Sigma}_{\bm{\beta}}^{-1} = b_{\bm{\beta}} \mathbf{I}_{k}} or 
#'   another list with valid tags \code{"M"} and \code{"P"} corresponding to 
#'   the whole mean vector \eqn{\bm{\mu}_{\bm{\beta}}} and the whole precision 
#'   matrix \eqn{\bm{\Sigma}_{\bm{\beta}}^{-1}}, respectively. The 
#'   \code{"sigma"} tag can be the vector \eqn{(a_{\sigma}, b_{\sigma})}.
#' @param starting a list with each tag corresponding to a parameter name. Valid 
#'   tags are \code{"beta"} and \code{"sigma"}.
#' @param n.samples the number of MCMC iterations after \code{n.burnin}.
#' @param n.thin the number of MCMC iterations kept is 1 out of \code{n.thin} 
#'   iterations.
#' @param n.burnin the number of MCMC iterations discarded at the beginning.
#' @param verbose if \code{TRUE}, model specification and progress of the 
#'   sampler is printed to the screen. Otherwise, nothing is printed to the 
#'   screen.
#' @param n.report the interval to report MCMC progress.
#' @param model,x,y logicals. If \code{TRUE} the corresponding components of 
#'   the fit (the model frame, the model matrix, the response) are returned.
#' @param ... currently no additional arguments.
#' 
#' @return An object of class \code{iidm}, which is a list comprising:
#'   \item{p.params.samples}{a \code{coda} object of posterior samples for the 
#'     model parameters, that is \eqn{\bm{\beta}} and \eqn{\sigma}.}
#'   \item{residuals}{a \code{coda} object of the residuals, that is response 
#'     minus fitted values.}
#'   \item{fitted.values}{a \code{coda} object of the fitted mean or quantile 
#'     values.}
#'   \item{method}{the method used.}
#'   \item{quantile}{if method="quantile", the quantile level used.}
#'   \item{xlevels}{(only where relevant) a record of the levels of the factors
#'     used in fitting.}
#'   \item{call}{the matched call.}
#'   \item{terms}{the \code{\link[stats]{terms}} object used.}
#'   \item{y}{if requested, the response used.}
#'   \item{x}{if requested, the model matrix used.}
#'   \item{model}{if requested (the default), the model frame used.}
#' 
#' @author Jorge Castillo-Mateo
#' 
#' @examples
#' ## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
#' ## Page 9: Plant Weight Data.
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
#' weight <- c(ctl, trt)
#' set.seed(12345)
#' iidm.D9 <- iidm(weight ~ group, verbose = FALSE)
#' set.seed(12345)
#' iidm.D90 <- iidm(weight ~ group - 1, verbose = FALSE) # omitting intercept
#'
#' @export
iidm <- function(
  formula,
  data,
  subset,
  method = c("mean", "quantile"),
  quantile = 0.5,
  priors = list("beta" = c(0, 1 / 100^2), "sigma" = c(2, 1)), 
  starting = list("beta" = 0, "sigma" = 1), 
  n.samples = 1000, 
  n.thin = 1, 
  n.burnin = 1000, 
  verbose = TRUE,
  n.report = 100, 
  model = TRUE, 
  x = FALSE, 
  y = FALSE,
  ...
) {
  
  method <- match.arg(method)
  method.mean <- method == "mean"
  if (
    !method.mean && (
    (length(quantile) != 1) ||
    (!is.numeric(quantile)) ||
    (quantile <= 0) || (quantile >= 1))
    ) stop("'quantile' should be a number in (0, 1)")
  if (!all(c("beta", "sigma") %in% names(priors)))
    stop("'priors' should have the valid tags 'beta' and 'sigma'")
  if (!all(c("beta", "sigma") %in% names(starting)))
    stop("'starting' should have the valid tags 'beta' and 'sigma'")
  
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- stats::model.response(mf, "numeric")
  N <- length(y)
  if (stats::is.empty.model(mt)) {
    x <- NULL
    z <- list(p.params.samples = numeric(), 
              residuals = y, 
              fitted.values = 0 * y,
              method = method)
  } else {
    x <- stats::model.matrix(mt, mf)
    k <- ncol(x)
    if (is.list(priors$beta)) {
      if (!all(c("M", "P") %in% names(priors$beta)))
        stop("'priors$beta' should have the valid tags 'M' and 'P'")
      if (length(priors$beta$M) != k)
        stop("'priors$beta$M' should have 'length' equal to the number of regression coefficients")
      if (!all(dim(priors$beta$P) == k))
        stop("'priors$beta$P' should have both 'dim' equal to the number of regression coefficients")
    } else {
      if (length(priors$beta) != 2)
        stop("'priors$beta' should have 'length' equal to 2")
      priors$beta <- list("M" = rep(priors$beta[1], k),
                          "P" = priors$beta[2] * diag(k))
    }
    if (!(length(starting$beta) %in% c(1, k)))
      stop("'starting$beta' should have 'length' equal to 1 or to the number of regression coefficients")
    if ((length(starting$beta) == 1) && (k != 1))
      starting$beta <- rep(starting$beta, k)
    if (length(starting$sigma) != 1)
      stop("'starting$sigma' should have 'length' equal to 1")
    #keep <- matrix(nrow = n.samples / n.thin, ncol = k + 1 + 2 * N)
    keep <- matrix(nrow = n.samples / n.thin, ncol = k + 1)
    if (!verbose) 
      n.report <- n.samples + 1
    
    if (method.mean) {
      params <- iidMeanRcpp(
        y,
        x,
        priors$beta$M,
        priors$beta$P,
        priors$sigma[1],
        priors$sigma[2],
        starting$beta,
        1 / starting$sigma^2,
        N,
        k,
        keep,
        n.samples,
        n.thin,
        n.burnin,
        n.report
      )
      
      params[, k + 1] <- 1 / sqrt(params[, k + 1])
      
    } else {
      params <- iidQuantileRcpp(
        quantile,
        y,
        x,
        priors$beta$M,
        priors$beta$P,
        priors$sigma[1],
        priors$sigma[2],
        starting$beta,
        1 / starting$sigma,
        N,
        k,
        keep,
        n.samples,
        n.thin,
        n.burnin,
        n.report
      )
      
      params[, k + 1] <- 1 / params[, k + 1]
      
    }
    
    #if (is.null(colnames(x))) {
    #  colnames(params) <- c(paste0("V", 1:k), "sigma", 1:N, 1:N)
    #} else {
    #  colnames(params) <- c(colnames(x), "sigma", 1:N, 1:N)
    #}
    
    if (is.null(colnames(x))) {
      colnames(params) <- c(paste0("V", 1:k), "sigma")
    } else {
      colnames(params) <- c(colnames(x), "sigma")
    }
    
    params <- coda::mcmc(params, 
                         start = n.burnin + 1, 
                         end = n.burnin + n.samples, 
                         thin = n.thin)
    
    #z <- list(
    #  p.params.samples = params[, 1:(k + 1)],
    #  residuals = params[, k + 1 + N + 1:N],
    #  fitted.values = params[, k + 1 + 1:N],
    #  method = method
    #)
    
    z <- list(
      p.params.samples = params,
      method = method
    )
  }
  
  if (!method.mean)
    z$quantile <- quantile
  z$xlevels <- stats::.getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model) 
    z$model <- mf
  if (ret.x) 
    z$x <- x
  if (ret.y) 
    z$y <- y
  
  class(z) <- "iidm"
  
  return(z)
}