#' @title Bayesian kriging for spatially varying coefficients
#' 
#' @description Bayesian kriging for a spatially varying coefficient in terms
#'   of the coefficients in the observation coordinates and the hyperparameters
#'   of the Gaussian process prior.
#'   
#' @description
#'   Function currently working but under development. 
#' 
#' \deqn{w(\mathbf{s}) \sim GP(\mu, C(\mathbf{s}, \mathbf{s}^{\prime}; \sigma^{2}, \phi))}
#' 
#' 
#'   
#'  
#' @param w a \eqn{B \times n} matrix with columns 
#'   \eqn{w(\mathbf{s}_{1}),\ldots,w(\mathbf{s}_{n})}.
#' @param hp a \eqn{B \times 3} matrix with columns \eqn{\mu}, \eqn{\sigma},
#'   and \eqn{\phi}.
#' @param coords a \eqn{n \times 2} matrix of the observation coordinates 
#'   in \eqn{D \subset R^{2}}.
#' @param newcoords a \eqn{n_{0} \times 2} matrix of the new coordinates to 
#'   krige, this is a grid \eqn{G_{D} \subset D \subset R^{2}}.
#' @param ... additional argument(s) for methods.
#' @return a \code{coda} object of posterior samples for the kriged 
#'   coefficient, that is \eqn{w(\mathbf{s}_{j})} for all
#'   \eqn{\mathbf{s}_{j} \in G_{D}}.
#'   
#' @author Jorge Castillo-Mateo
#' 
#' @export 
krigeBayes <- function(w, hp, coords, newcoords, ...) {
  
  B  <- nrow(w)
  n  <- nrow(coords)
  n0 <- nrow(newcoords)
  w0 <- matrix(nrow = B, ncol = n0)
  
  d <- spTReg:::dist1(rbind(newcoords, coords))
  
  for (b in 1:B) {
    print(paste("Iteration", b, "out of", B))
    R22 <- exp(- hp[b, 3] * d[n0 + 1:n, n0 + 1:n])
    R22inv <- solve(R22)
    R11 <- exp(- hp[b, 3] * d[1:n0, 1:n0])
    R12 <- exp(- hp[b, 3] * d[1:n0, n0 + 1:n])
    R12R22inv <- R12 %*% R22inv
    R <- t(chol(R11 - R12R22inv %*% t(R12))) #lower triangular
    
    w0[b,] <- 
      hp[b, 1] + R12R22inv %*% (w[b,] - hp[b, 1]) + 
      R %*% rnorm(n0) * hp[b, 2]
  }
  
  return(w0)
}