#' @title MCMC sampling for the space-time models
#' @description This function is used to draw MCMC samples using the 
#'   Metropolis-within-Gibbs sampler that fits the Bayesian models.
#'   
#' @param y A numeric vector, or the name of a column of \code{data} to 
#'   which the model is to be fit.
#' @param data A data frame object with named columns giving the data to be 
#'   fit; any explanatory variable necessary for modeling any of the 
#'   parameters; and named columns \code{t}, \code{l}, and \code{i}, 
#'   representing the numeric index of the data in the long and short 
#'   temporal scales, and space, respectively. If not found in \code{data}, the
#'   variables are taken from \code{environment(formula)}, typically the 
#'   environment from which \code{spTm} is called.
#' @param location.fun,ar.fun,scale.fun An object of class 
#'   \code{\link{formula}} describing a model for each 
#'   parameter using columns from \code{data}. \code{data} must be supplied if 
#'   any of these arguments have an explanatory variable. Three options are 
#'   currently supported for \code{ar.fun}: the default \code{~ -1} with no 
#'   autoregressive term, \code{~ 1} with a common autoregressive term, and 
#'   \code{~ sp(1)} with a spatially-varying autoregressive term. Two options 
#'   are currently supported for \code{scale.fun}: the default \code{~ 1} with 
#'   a common scale term, and \code{~ sp(1)} with a spatially-varying scale 
#'   term.
#' @param coords an \eqn{n x 2} matrix of the observation coordinates in 
#'   \code{R^2} (e.g., easting and northing).
#' @param priors a list with each tag corresponding to a parameter name. Valid 
#'   tags are...?
#' @param starting a list with each tag corresponding to a parameter name. Valid 
#'   tags are...?
#' @param tuning a list with each tag corresponding to a parameter name. Valid 
#'   tags are...?
#' @param center.scale If \code{TRUE}, non-constant columns of \eqn{X} are 
#'   centered on zero and scaled to have variance one. If 
#'   \code{\link{predict.spTm}} is subsequently called this centering and 
#'   scaling is applied automatically.
#' @param n.samples The number of MCMC iterations after \code{n.burnin}.
#' @param n.thin The number of MCMC iterations kept is 1 out of \code{n.thin} 
#'   iterations.
#' @param n.burnin The number of MCMC iterations discarded at the beginning.
#' @param verbose If \code{TRUE}, model specification and progress of the 
#'   sampler is printed to the screen. Otherwise, nothing is printed to the 
#'   screen.
#' @param n.report The interval to report Metropolis sampler acceptance and 
#'   MCMC progress.
#' @param ... currently no additional arguments.
#' 
#' @return An object of class \code{spTm}, which is a list comprising:
#'   \item{coords}{the \eqn{n x 2} matrix specified by \code{coords}.}
#'   \item{p.params.samples}{A \code{coda} object of posterior samples for the 
#'     defined parameters.}
#'   \item{acceptance}{The Metropolis sampling acceptance percent.}
#'   
#'   The return object will include additional objects used for subsequent 
#'   Gaussian process kriging, data kriging, and model fit evaluation using 
#'   \code{\link{predict.spTm}}, \code{\link{validation.spTm}}, 
#'   respectively.
#' 
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{predict.spTm}}, \code{\link{predict.spTm}}
#' @references 
#' Castillo-Mateo J, Lafuente M, Asín J, Cebrián AC, Gelfand AE, Abaurrea J (2022). 
#' Spatial modeling of day-within-year temperature time series: an examination of daily maximum temperatures in Aragón, Spain. 
#' \emph{Journal of Agricultural, Biological and Environmental Statistics}, \strong{27}(3), 487--505. 
#' \doi{10.1007/s13253-022-00493-3}.
#' 
#' Castillo-Mateo J, Asín J, Cebrián AC, Gelfand AE, Abaurrea J. 
#' Spatial quantile autoregression for season within year daily maximum temperature data. 
#' \emph{Annals of Applied Statistics}, \strong{17}(3), 2305--2325. 
#' \doi{10.1214/22-AOAS1719}
#' 
#' @examples
#' 1 + 1
#'
#' @export
spTm <- function(
  y, data,
  model = c("mean", "quantile"),
  location.fun = ~ 1,
  ar.fun = ~ -1,
  scale.fun = ~ 1,
  coords, priors, starting, tuning, 
  center.scale = FALSE, 
  n.samples = 10000, n.thin = 1, n.burnin = 0, 
  verbose = TRUE, n.report = 100, ...) {
  
  
  
}