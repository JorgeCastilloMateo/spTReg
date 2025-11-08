#' @title Bayesian kriging for spatially varying coefficients
#' 
#' @description Bayesian kriging for a spatially varying coefficient in terms
#'   of the coefficients in the observation coordinates and the hyperparameters
#'   of the Gaussian process prior.
#'   
#' @details
#'   The Gaussian process with exponential covariance function is defined as: 
#'   
#'   \deqn{w(\mathbf{s}) \sim GP(\mu, C(\mathbf{s}, \mathbf{s}^{\prime}; \sigma^{2}, \phi)).}
#' 
#'   Parallel execution is only available if the package was compiled with 
#'   \code{OpenMP} support. Parallel and sequential execution are both 
#'   reproducible with \code{\link{set.seed}}, but the numerical results will 
#'   differ between the two modes because the underlying random-number streams 
#'   are independent.
#'  
#' @param w a \eqn{B \times n} matrix with columns 
#'   \eqn{w(\mathbf{s}_{1}),\ldots,w(\mathbf{s}_{n})}.
#' @param hp a \eqn{B \times 3} matrix with columns \eqn{\mu}, \eqn{\sigma},
#'   and \eqn{\phi}.
#' @param coords a \eqn{n \times 2} matrix of the observation coordinates 
#'   in \eqn{D \subset R^{2}}.
#' @param newcoords a \eqn{n_{0} \times 2} matrix of the new coordinates to 
#'   krige, this is a grid \eqn{G_{D} \subset D \subset R^{2}}.
#' @param parallel logical; whether to run the computation in parallel 
#'   using \code{OpenMP} (default = \code{FALSE}). See `Details'.
#' @param n.threads integer; number of threads to use when 
#'   \code{parallel = TRUE}.  
#'   The default \code{0} uses all available physical cores.  
#'   This argument is ignored when \code{parallel = FALSE}.
#' @param ... additional argument(s) for methods.
#' @return a \eqn{B \times n_{0}} matrix or \code{coda} object of posterior 
#'   samples for the kriged coefficient, that is \eqn{w(\mathbf{s}_{j})} for
#'   all \eqn{\mathbf{s}_{j} \in G_{D}}.
#'   
#' @author Jorge Castillo-Mateo
#' 
#' @export 
krigeBayes <- function(
    w,
    hp,
    coords,
    newcoords,
    parallel = FALSE,
    n.threads = 0,
    ...) {

  if (!is.matrix(w)) stop("'w' must be a matrix")
  if (!is.matrix(hp) || ncol(hp) != 3)
    stop("'hp' must have 3 columns (mu, sigma, decay)")
  if (!is.matrix(coords) || ncol(coords) != 2)
    stop("'coords' must be a matrix with exactly 2 columns.")
  if (!is.matrix(newcoords) || ncol(newcoords) != 2)
    stop("'newcoords' must be a matrix with exactly 2 columns.")
  
  krigeBayesRcpp(
    w, hp, coords, newcoords,
    parallel = parallel,
    nThreads = n.threads
  )
}