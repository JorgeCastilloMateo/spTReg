#' @title Validation metrics for the space-time models
#' @description This function is used to obtain validation metrics using MCMC
#'   samples.
#'   
#' @param model Object of class "spTm".
#'  
#' @return 
#'   \item{MSE}{Mean Squared Error.}
#'   \item{RMSE}{Root Mean Squared Error.}
#'   \item{MAE}{Mean Absolute Error.}
#'   \item{CRPS}{Continuous Rank Probability Score.}
#'   \item{CVG}{Coverage.}
#'   
#'   \item{WMAE}{Weighted Mean Absolute Error.}
#'   \item{p}{Probability of being under a quantile.}
#'   
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{spTm}}, \code{\link{predict.spTm}}
#' @export
validation.spTm <- function(model) {
  
  return(1)
}