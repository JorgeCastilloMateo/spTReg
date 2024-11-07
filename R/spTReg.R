#' @title \strong{spTReg}: Spatio-Temporal Mean, Quantile, and Binary auto-Regression Models
#' @aliases spTReg-package
#' @aliases spTReg
#' @description 
#'   \strong{spTReg} uses different hierarchical Bayesian spatio-temporal 
#'   modeling strategies with time- and spatially-varying coefficients, namely:
#'   
#'   (1) Mean regression using Gaussian processes and Gaussian errors,
#'   
#'   (2) Quantile regression using Gaussian processes and asymmetric Laplace errors,
#'   
#'   (3) Probit and logit binary regression using Gaussian processes.
#'   
#' @details 
#'   The back-end code of this package is built under \code{C++} language.
#'   
#'   Main functions used: 
#'   
#'   (1) \code{\link{iidm}} and \code{\link{predict.iidm}} fit and predict 
#'   models with independent and identically distributed errors.
#'   
#'   (2) \code{\link{spTm}} and \code{\link{predict.spTm}} fit and predict
#'   models for spatio-temporal dependent data.
#'   
#' @docType package
#' @author Jorge Castillo-Mateo <jorgecastillomateo@gmail.com>
#' @references 
#' #' Castillo-Mateo J, Asín J, Cebrián AC, Gelfand AE, Abaurrea J (2023). 
#' Spatial quantile autoregression for season within year daily maximum temperature data. 
#' \emph{Annals of Applied Statistics}, \strong{17}(3), 2305--2325. 
#' \doi{10.1214/22-AOAS1719}
#' 
#' Castillo-Mateo J, Gelfand AE, Gracia-Tabuenca Z, Asín J, Cebrián AC (\emph{in press}).
#' Spatio-temporal modeling for record-breaking temperature events in Spain. 
#' \emph{Journal of the American Statistical Association}.
#' 
#' Castillo-Mateo J, Lafuente M, Asín J, Cebrián AC, Gelfand AE, Abaurrea J (2022). 
#' Spatial modeling of day-within-year temperature time series: an examination of daily maximum temperatures in Aragón, Spain. 
#' \emph{Journal of Agricultural, Biological and Environmental Statistics}, \strong{27}(3), 487--505. 
#' \doi{10.1007/s13253-022-00493-3}.
#' 
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib spTReg
#' @name spTReg-package
NULL
