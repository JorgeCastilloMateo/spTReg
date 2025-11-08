#' @title Printing model fits
#' @aliases print.iidm
#'
#' @description \code{print} method for class \code{"iidm"}.
#'
#' @details
#' Displays a compact summary of an \code{"iidm"} model fit, including the 
#' original model call and the posterior means of the model parameters. 
#' For a more detailed output including diagnostics and posterior summaries,
#' see \code{\link{summary.iidm}}.
#'
#' @param x an object of class \code{"iidm"}, usually the result of a call to 
#'   \code{\link{iidm}}.
#' @param digits the number of significant digits to display (default is taken from
#'   \code{getOption("digits") - 3L}).
#' @param ... further arguments passed to or from other methods.
#'
#' @return The function \code{print.iidm} is used for its side effect: it prints
#'   a concise textual representation of the fitted model to the console.
#'   It returns the input object \code{x}, invisibly.
#'
#' @seealso \code{\link{iidm}}, \code{\link{summary.iidm}}
#'   
#' @author Jorge Castillo-Mateo
#' 
#' @rdname print
#' @method print iidm
#' @export print.iidm
#' @export 
print.iidm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  #cat("Method:", x$method)
  #if (!is.null(x$quantile))
  #  cat(" (τ = ", x$quantile, ")", sep = "")
  cat("\nCoefficients (posterior means):\n")
  print.default(format(colMeans(as.matrix(x$p.params.samples)), digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}
