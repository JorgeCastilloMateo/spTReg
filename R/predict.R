#' @title Spatial and temporal predictions for the space-time models
#' @description This function is used to obtain spatial predictions and also to
#'   get the temporal forecasts using MCMC samples.
#'   
#' @param model Object of class "spTm".
#' @param newdata An optional data set providing the explanatory variables for
#'   spatial prediction or temporal forecasts. This data should have the same 
#'   space-time structure as the original \code{data}, see \code{\link{spTm}}.
#'   The index columns \code{t}, \code{l}, and \code{i} must follow the 
#'   numbering started at \code{data}, e.g., the first new location must be 
#'   $n+1$. If omitted, the fitted values are used.
#' @param newcoords The coordinates for the prediction or forecast sites. The 
#'   locations are in similar format to \code{coords}, see \code{\link{spTm}}.
#'   If omitted, the fitted values are used.
#' 
#' @return A matrix that holds the response variable posterior predictive 
#'   samples (or posterior quantile samples). The rows of this matrix 
#'   correspond to the posterior predictive samples and the columns are the 
#'   predicted values. 
#'   
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{spTm}}, \code{\link{validation.spTm}}
#' @export
predict.spTm <- function(model, newdata, newcoords) {
  
  return(1)
}
