#' @title Space-varying coefficient in the formula
#' @description This function is used to define space-varying coefficients 
#'   within the formula for the space-time models.
#'   
#' @param x The explanatory variable for which space-varying coefficients are 
#'   defined.
#' 
#' @return An object of class \code{spCoef}.
#' 
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{spTm}}, \code{\link{tp}}
#' @export
sp <- function(x) {
  class(x) <- "spCoef"
  return(x)
}
