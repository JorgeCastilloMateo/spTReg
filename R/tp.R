#' @title Time-varying coefficient in the formula
#' @description This function is used to define time-varying coefficients 
#'   within the formula for the space-time models.
#'   
#' @param x The explanatory variable for which time-varying coefficients are 
#'   defined.
#' 
#' @return An object of class \code{tpCoef}.
#' 
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{spTm}}, \code{\link{sp}}
#' @export
tp <- function(x) {
  class(x) <- "tpCoef"
  return(x)
}
