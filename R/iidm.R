#' @title MCMC sampling for the i.i.d. models
#' @description This function is used to draw MCMC samples using the Gibbs 
#'   sampler that fits the Bayesian linear regression models with i.i.d. 
#'   errors.
#'   
#' @details
#'   The function \code{iidm} can fit two types of model, depending on whether
#'   the mean or the quantile is of interest.
#'   
#'   General model:
#'   \deqn{Y_{i} = \mathbf{x}_{i} \bm{\beta} + \epsilon_{i}, \qquad i=1,\ldots,n,}
#'   with \eqn{\epsilon_i \sim \text{i.i.d. } N(0, \sigma^{2})} for the mean, or 
#'   \eqn{\epsilon_i \sim \text{i.i.d. } AL(0, \sigma, \tau)} for the 
#'   \eqn{\tau}-th quantile.
#'   
#'   Priors:
#'   \deqn{\bm{\beta} \sim N_{p}(\bm{\mu}_{\bm{\beta}}, \bm{\Sigma}_{\bm{\beta}})}
#'   \deqn{\sigma^{2} \text{ (mean)}, \sigma \text{ (quantile)} \sim IG(a_{\sigma}, b_{\sigma})}
#' 
#'   Parallel execution is only available if the package was compiled with 
#'   \code{OpenMP} support. This mode samples the exponential latent variables 
#'   for quantile regression in parallel, but the current implementation does 
#'   not improve computational times. Parallel and sequential execution are
#'   both reproducible with \code{\link{set.seed}}, and the numerical results 
#'   between the two modes should be the same.
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
#'   the \code{\link[stats]{model.frame}} documentation.
#' @param na.action a function indicating how to handle \code{NA}s; 
#'   \code{na.omit} deletes cases, \code{na.pass} keeps and imputes them within 
#'   the MCMC (similar to \code{NULL}). See \code{\link[base]{options}}.
#' @param method the method to be used for fitting, \code{"mean"} or 
#'   \code{"quantile"} regression (default=\code{"mean"}).
#' @param quantile a numeric scalar containing the quantile level of interest 
#'   (default=0.5).
#' @param priors a list with each tag corresponding to a parameter name. Valid 
#'   tags are \code{"beta"} and \code{"sigma"}. The \code{"beta"} tag can be 
#'   either the vector \eqn{(a_{\bm{\beta}}, b_{\bm{\beta}})} such that 
#'   \eqn{\bm{\mu}_{\bm{\beta}} = a_{\bm{\beta}} \mathbf{1}_{p}} and 
#'   \eqn{\bm{\Sigma}_{\bm{\beta}}^{-1} = b_{\bm{\beta}} \mathbf{I}_{p}} or 
#'   another list with valid tags \code{"M"} and \code{"P"} corresponding to 
#'   the whole mean vector \eqn{\bm{\mu}_{\bm{\beta}}} and the whole precision 
#'   matrix \eqn{\bm{\Sigma}_{\bm{\beta}}^{-1}}, respectively. The 
#'   \code{"sigma"} tag must be the vector \eqn{(a_{\sigma}, b_{\sigma})}. 
#'   (See the `Details' section to check the specific notation.)
#' @param starting a list with each tag corresponding to a parameter name. Valid 
#'   tags are \code{"beta"} and \code{"sigma"}. The \code{"beta"} tag can be 
#'   either a number \eqn{\beta^{(0)}} such that 
#'   \eqn{\bm{\beta}^{(0)} = \beta^{(0)} \textbf{1}_{p}} or directly the vector
#'   \eqn{\bm{\beta}^{(0)}}. The \code{"sigma"} tag must be a positive number 
#'   \eqn{\sigma^{(0)}}.
#' @param n.samples the number of MCMC iterations after \code{n.burnin}.
#' @param n.thin the number of MCMC iterations kept is 1 out of \code{n.thin} 
#'   iterations.
#' @param n.burnin the number of MCMC iterations discarded at the beginning.
#' @param verbose if \code{TRUE}, model specification and progress of the 
#'   sampler is printed to the screen. Otherwise, nothing is printed to the 
#'   screen.
#' @param n.report the interval to report MCMC progress.
#' @param parallel logical; if method="quantile", whether to run the 
#'   computation in parallel using \code{OpenMP} (default = \code{FALSE}).
#'   See `Details'.
#' @param n.threads integer; number of threads to use when 
#'   \code{parallel = TRUE}.  
#'   The default \code{0} uses all available physical cores.  
#'   This argument is ignored when \code{parallel = FALSE}.
#' @param model,x,y logicals; if \code{TRUE} the corresponding components of 
#'   the fit (the model frame, the model matrix, the response) are returned.
#' @param ... currently no additional arguments.
#' 
#' @return An object of class \code{iidm}, which is a list comprising:
#'   \item{p.params.samples}{a \code{coda} object of posterior samples for the 
#'     model parameters, that is \eqn{\bm{\beta}} and \eqn{\sigma}.}
#'   \item{Y.missing.samples}{a \code{coda} object of posterior samples for the 
#'     imputed missing response values.}
#'   \item{mcmc}{a list with information of the MCMC chain, that is 
#'     \code{n.samples}, \code{n.thin}, \code{n.burnin}, and CPU \code{time} 
#'     used.}
#'   \item{method}{the method used.}
#'   \item{quantile}{if method="quantile", the quantile level used.}
#'   \item{rank}{the numeric rank of the fitted linear model.}
#'   \item{xlevels}{(where relevant) a record of the levels of the factors
#'     used in fitting.}
#'   \item{call}{the matched call.}
#'   \item{terms}{the \code{\link[stats]{terms}} object used.}
#'   \item{y}{if requested, the response used.}
#'   \item{x}{if requested, the model matrix used.}
#'   \item{model}{if requested (the default), the model frame used.}
#'   \item{na.action}{(where relevant) information returned by 
#'     \code{\link[stats]{model.frame}} on the special handling of \code{NA}s.}
#' 
#' @author Jorge Castillo-Mateo
#' 
#' @seealso \code{\link{confint.iidm}}, \code{\link{predict.iidm}}, 
#'   \code{\link{print.iidm}}, \code{\link{summary.iidm}}
#' 
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
#' @examples
#' ## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
#' ## Page 9: Plant Weight Data.
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
#' weight <- c(ctl, trt)
#' set.seed(12345)
#' iidm.D9 <- iidm(weight ~ group)
#' set.seed(12345)
#' iidm.D90 <- iidm(weight ~ group - 1) # omitting intercept
#'
#' @export
iidm <- function(
  formula,
  data,
  subset,
  na.action,
  method = c("mean", "quantile"),
  quantile = 0.5,
  priors = list("beta" = c(0, 1e-08), "sigma" = c(0.1, 0.1)), 
  starting = list("beta" = 0.01, "sigma" = 1), 
  n.samples = 1000, 
  n.thin = 1, 
  n.burnin = 1000, 
  verbose = FALSE,
  n.report = 100,
  parallel = FALSE,
  n.threads = 0,
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

  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- stats::model.response(mf, "numeric")
  N <- length(y)
  if (stats::is.empty.model(mt)) {
    x <- NULL
    p <- 0L
    z <- list(p.params.samples = numeric(), 
              method = method)
  } else {
    x <- stats::model.matrix(mt, mf)
    p <- ncol(x)
    
    priors   <- .check_priors(priors, p)
    starting <- .check_starting(starting, p)
    
    if (!verbose) 
      n.report <- 0
    
    time <- Sys.time()
    if (method.mean) {
      res <- iidMeanRcpp(
        y,
        x,
        priors$beta$M,
        priors$beta$P,
        priors$sigma[1],
        priors$sigma[2],
        starting$beta,
        1 / starting$sigma^2,
        N,
        p,
        n.samples,
        n.thin,
        n.burnin,
        n.report
      )
      res$params[, p + 1] <- 1 / sqrt(res$params[, p + 1])
    } else {
      res <- iidQuantileRcpp(
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
        p,
        n.samples,
        n.thin,
        n.burnin,
        n.report,
        parallel,
        n.threads
      )
      res$params[, p + 1] <- 1 / res$params[, p + 1]
    }
    time <- Sys.time() - time
    
    if (is.null(colnames(x))) {
      colnames(res$params) <- c(paste0("V", 1:p), "sigma")
    } else {
      colnames(res$params) <- c(colnames(x), "sigma")
    }
    
    res$params <- coda::mcmc(
      res$params, 
      start = n.burnin + 1, 
      end = n.burnin + n.samples, 
      thin = n.thin)
    
    z <- list(p.params.samples = res$params)
    if (length(res) > 1) {
      colnames(res$missing) <- which(!is.finite(y))
      res$missing <- coda::mcmc(
        res$missing, 
        start = n.burnin + 1, 
        end = n.burnin + n.samples, 
        thin = n.thin)
      z$Y.missing.samples <- res$missing
    }    
    z$mcmc <- list(
      n.samples = n.samples, 
      n.thin = n.thin, 
      n.burnin = n.burnin, 
      time = time)
    z$method <- method
  }
  
  if (!method.mean)
    z$quantile <- quantile
  z$rank <- p
  z$na.action <- attr(mf, "na.action")
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
