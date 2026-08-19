#' @title Spatially-varying coefficient specification
#' @description Defines a spatially varying coefficient term to be used inside 
#'   a model formula for \code{\link{spTm}}.
#'
#' @details
#'   The function \code{sp()} creates an object that represents a spatially 
#'   varying coefficient of the form
#'
#'   \deqn{\alpha(s) = \mu(s) + \eta(s),}
#'
#'   where
#'
#'   \itemize{
#'     \item \eqn{\mu(s)} is a spatial mean defined by regressors in a formula
#'         such as \code{z ~ w1 + w2}.
#'     \item \eqn{\eta(s)} is a centered Gaussian process with a covariance
#'         kernel chosen by the user (exponential, Gaussian, Matérn, etc.).
#'   }
#'
#'   The left-hand side \code{z} multiplies the spatial coefficient inside the
#'   regression model; that is, the linear predictor includes a term
#'
#'   \deqn{z_{t\ell}(s_i)\,\alpha(s_i).}
#'
#'   The right-hand side regressors define the prior mean of \eqn{\alpha(s)}.
#'
#'   Example usage inside a formula:
#'
#'   \preformatted{
#'   y ~ x1 + x2 + sp(z ~ elev + slope, kernel="exponential", range=30)
#'   }
#'
#' @param formula a formula of the form \code{z ~ w1 + w2}, where the left-hand
#'   side is the covariate multiplied by the spatially varying coefficient and
#'   the right-hand side are spatial covariates defining the mean function
#'   \eqn{\mu(s)}.
#' @param kernel character string specifying the correlation kernel:
#'   currently available only \code{"exponential"}.
#' @param scale optional fixed standard deviation parameter. If 
#'   \code{NULL} (default), the scale is estimated in the MCMC.
#' @param decay optional fixed spatial decay (1 / range) parameter. If 
#'   \code{NULL} (default), the decay is estimated in the MCMC.
#'
#' @return An object of class \code{spCoef}.
#'
#' @seealso \code{\link{spTm}}, \code{\link{tp}}
#'
#' @examples
#' # Include a spatially varying coefficient of z with mean depending on an intercept and elevation:
#' sp(z ~ elev)
#'
#' # Inside a model formula:
#' # y ~ x1 + sp(z ~ elev)
#'
#' @export
sp <- function(formula,
               kernel = c("exponential"),
               scale = NULL,
               decay = NULL) {
  
  if (!inherits(formula, "formula"))
    stop("'formula' must be of the form z ~ w1 + w2")
  
  kernel <- match.arg(kernel)
  
  structure(
    list(
      formula = formula,
      kernel  = kernel,
      scale   = scale,
      decay   = decay
    ),
    class = "sp_term"
  )
}

# -------------------------------------------------------------------------
# Register the special term in formula processing
# -------------------------------------------------------------------------

# Tell R that "sp" is a special term (like s() in mgcv)
terms.sp_term <- function(x, ...) x

# Needed so that model.frame recognizes the special
model.frame.sp_term <- function(formula, data, ...) {
  # Expand the internal formula z ~ w
  mf <- model.frame(formula$formula, data = data, ...)
  mf
}

# Construct the design matrix for the mean μ(s)
model.matrix.sp_term <- function(object, data, contrasts.arg = NULL, ...) {
  mm <- model.matrix(object$formula, data = data,
                     contrasts.arg = contrasts.arg)
  attr(mm, "assign") <- NULL
  mm
}
