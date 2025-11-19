#' @title MCMC sampling for the space-time models
#' @description This function is used to draw MCMC samples using the 
#'   Metropolis-within-Gibbs sampler that fits the Bayesian spatio-temporal
#'   models.
#'
#' @description
#'   Function currently working but under development.
#'   
#'   Model:
#'   
#'   \deqn{\mathbf{Y}_{t\ell} = \mathbf{O}_{t\ell} + \mathbf{P} (\mathbf{Y}_{t,\ell-1} - \mathbf{O}_{t,\ell-1}) + \bm{\epsilon}_{t\ell},}
#'   \deqn{\mathbf{O}_{t\ell} = \mathbf{X}_{t\ell} \bm{\beta} + \mathbf{U}_{t\ell} \bm{\gamma}_{t} + \sum_{m=1}^{r} \mathbf{V}_{t\ell,m} \bm{\alpha}_{m},}
#'   \deqn{\bm{\epsilon}_{t\ell} \sim \text{i.i.d. } N_{n}(\mathbf{0}_{n}, \bm{\Sigma}).}
#'
#'   Parallel execution is only available if the package was compiled with 
#'   \code{OpenMP} support. This mode samples the exponential latent variables 
#'   for quantile regression in parallel, but the current implementation does 
#'   not improve computational times. Parallel and sequential execution are
#'   both reproducible with \code{\link{set.seed}}, and the numerical results 
#'   between the two modes should be the same.
#'   
#' @param formula an object of class \code{"\link[stats]{formula}"} (or one 
#'   that can be coerced to that class): a symbolic description of the model to
#'   be fitted. The details of model specification are given under ‘Details’.
#' @param u,v (in addition to \code{formula}) \eqn{N \times q} and 
#'   \eqn{N \times r} matrices of regression variables accompanying a 
#'   time-varying coefficient and a spatially-varying coefficient, 
#'   respectively.
#' @param data a data frame object with named columns giving the data to be 
#'   fit; any explanatory variable necessary for modeling any of the 
#'   parameters; and named columns \code{t}, \code{l}, and \code{s}, 
#'   representing the numeric index of the data in the long and short 
#'   temporal scales, and space, respectively; sorted by space. If not found 
#'   in \code{data}, the variables are taken from \code{environment(formula)}, 
#'   typically the environment from which \code{spTm} is called.
#' @param location.fun,ar.fun,scale.fun (not working) an object of class 
#'   \code{\link{formula}} describing a model for each 
#'   parameter using columns from \code{data}. \code{data} must be supplied if 
#'   any of these arguments have an explanatory variable. Three options are 
#'   currently supported for \code{ar.fun}: the default \code{~ -1} with no 
#'   autoregressive term, \code{~ 1} with a common autoregressive term, and 
#'   \code{~ sp(1)} with a spatially-varying autoregressive term. Two options 
#'   are currently supported for \code{scale.fun}: the default \code{~ 1} with 
#'   a common scale term, and \code{~ sp(1)} with a spatially-varying scale 
#'   term.
#' @param coords a \eqn{n \times 2} matrix of the observation coordinates in 
#'   \code{R^2} (e.g., easting and northing).
#' @param priors a list with each tag corresponding to a parameter name. 
#' @param starting a list with each tag corresponding to a parameter name. 
#' @param center.scale (not currently available) if \code{TRUE}, non-constant 
#'   columns of \eqn{X} are 
#'   centered on zero and scaled to have variance one. If 
#'   \code{\link{predict.spTm}} is subsequently called this centering and 
#'   scaling is applied automatically.
#' @param n.samples the number of MCMC iterations after \code{n.burnin}.
#' @param n.thin the number of MCMC iterations kept is 1 out of \code{n.thin} 
#'   iterations.
#' @param n.burnin the number of MCMC iterations discarded at the beginning.
#' @param verbose if \code{TRUE}, model specification and progress of the 
#'   sampler is printed to the screen. Otherwise, nothing is printed to the 
#'   screen.
#' @param n.report the interval to report Metropolis sampler acceptance and 
#'   MCMC progress.
#' @param ... currently no additional arguments.
#' 
#' @return An object of class \code{spTm}, which is a list comprising (this 
#'   list of outputs must be updated):
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
#' 
#' @seealso \code{\link{confint.spTm}}, \code{\link{predict.spTm}}
#' 
#' @references 
#' Castillo-Mateo J, Lafuente M, Asín J, Cebrián AC, Gelfand AE, Abaurrea J (2022).
#' Spatial modeling of day-within-year temperature time series: an examination of daily maximum temperatures in Aragón, Spain. 
#' \emph{Journal of Agricultural, Biological and Environmental Statistics}, \strong{27}(3), 487--505. 
#' \doi{10.1007/s13253-022-00493-3}.
#' 
#' Castillo-Mateo J, Asín J, Cebrián AC, Gelfand AE, Abaurrea J (2023).
#' Spatial quantile autoregression for season within year daily maximum temperature data. 
#' \emph{Annals of Applied Statistics}, \strong{17}(3), 2305--2325. 
#' \doi{10.1214/22-AOAS1719}
#'
#' @export
spTm <- function(
  formula, 
  u,
  v,
  x_alpha, # list of design matrix for each process
  data, 
  subset,
  na.action,
  method = c("mean", "quantile"),
  quantile = 0.5,
  #ar.fun = ~ -1,
  #scale.fun = ~ 1,
  coords, 
  priors = list("beta" = c(0, 1 / 100), 
                "sigma" = c(0.1, 0.1),
                "rho" = c(0, 1 / 100),
                "hpGamma" = list("sigma" = c(0.1, 0.1), "rho" = c(0, 1 / 100)),
                "hpAlpha" = list("prec0" = 1 / 100, "sigma" = c(0.1, 0.1), "phi" = c(2, 100)),
                "hpSigma" = list("mu" = 0, "sigma" = c(0.1, 0.1), "phi" = c(2, 100)),
                "hpRho" = list("mu" = 0, "sigma" = c(0.1, 0.1), "phi" = c(2, 100))),
  starting = list("beta" = 0.01, 
                  "sigma" = 1, 
                  "rho" = 0,
                  "gamma" = 0, 
                  "alpha" = 0, 
                  "hpGamma" = c("sigma" = 1, "phi" = 3 / 100),
                  "hpAlpha" = c("sigma" = 1, "phi" = 3 / 100),
                  "hpSigma" = c("mu" = 0, "sigma" = 1, "phi" = 3 / 100),
                  "hpRho" = c("mu" = 0, "sigma" = 1, "phi" = 3 / 100)),
  n.samples = 1000, 
  n.thin = 1, 
  n.burnin = 1000, 
  verbose = FALSE,
  n.report = 100,
  parallel = FALSE,
  n.threads = 0,
  model = TRUE, 
  x = FALSE, 
  y = FALSE,
  ...) {
  
  method <- match.arg(method)
  method.mean <- method == "mean"
  if (
    !method.mean && (
      (length(quantile) != 1) ||
      (!is.numeric(quantile)) ||
      (quantile <= 0) || (quantile >= 1))
  ) stop("'quantile' should be a number in (0, 1)")
  if (!all(c("beta", "sigma", "phi", "mu") %in% names(priors)))
    stop("'priors' should have the valid tags 'beta', 'sigma', 'phi', and 'mu'")
  if (!all(c("beta", "sigma", "alpha", "hp") %in% names(starting)))
    stop("'starting' should have the valid tags 'beta', 'sigma', 'phi', and 'mu'")
  
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- stats::model.response(mf, "numeric")
  #attr(mt, "dataClasses") <- sapply(mf, class)
  N <- length(y)
  if (stats::is.empty.model(mt)) {
    x <- NULL
    z <- list(p.params.samples = numeric(), 
              residuals = y, 
              fitted.values = 0 * y,
              method = method)
  } else {
    x <- stats::model.matrix(mt, mf)
    #u <-
    #v <- 
    n <- nrow(coords)
    p <- ncol(x)
    if (missing(u) || is.null(u)) {
      q <- 0
    } else {
      q <- ncol(u) #sum(attr(mt, "dataClasses") == "tpCoef")
    }
    if (missing(v) || is.null(v)) {
      r <- 0
    } else {
      r <- ncol(v) #sum(attr(mt, "dataClasses") == "spCoef")
    }
    s <- rep(1:n, each = N / n)
    
    # check priors and starting of beta
    if (is.list(priors$beta)) {
      if (!all(c("M", "P") %in% names(priors$beta)))
        stop("'priors$beta' should have the valid tags 'M' and 'P'")
      if (length(priors$beta$M) != p)
        stop("'priors$beta$M' should have 'length' equal to the number of regression coefficients")
      if (!all(dim(priors$beta$P) == p))
        stop("'priors$beta$P' should have both 'dim' equal to the number of regression coefficients")
    } else {
      if (length(priors$beta) != 2)
        stop("'priors$beta' should have 'length' equal to 2")
      priors$beta <- list("M" = rep(priors$beta[1], p),
                          "P" = priors$beta[2] * diag(p))
    }
    if (!(length(starting$beta) %in% c(1, p)))
      stop("'starting$beta' should have 'length' equal to 1 or to the number of regression coefficients")
    if ((length(starting$beta) == 1) && (p != 1))
      starting$beta <- rep(starting$beta, p)
    if (length(starting$sigma) != 1)
      stop("'starting$sigma' should have 'length' equal to 1")
    # if () comprobar que starting$betat sea un numero o matriz de tamaño correcto
    # if () comprobar que starting$alpha sea un numero o matriz de tamaño correcto
    if (length(starting$alpha) == 1) 
      starting$alpha <- matrix(starting$alpha, nrow = n, ncol = r)
    if (length(starting$hp) == 3) 
      starting$hp <- matrix(starting$hp, nrow = 3, ncol = r)
    if (!verbose) 
      n.report <- n.samples + 1
    
    ### provisional ###
    p_alpha <- sapply(x_alpha, ncol, simplify = TRUE)
    M_beta_alpha <- list()
    P_beta_alpha <- list()
    beta_alpha <- list()
    for (m in 1:r) {
      M_beta_alpha[[m]] <- rep(0, p_alpha[m])
      P_beta_alpha[[m]] <- priors$mu[2] * diag(p_alpha[m])
      beta_alpha[[m]] <- rep(0, p_alpha[m])
    }
    ### ###
    
    dist <- dist1(coords)
    
    time <- Sys.time()
    if (method.mean) {
      res <- spMeanRcpp(
        y,
        x,
        u,
        v,
        x_alpha,
        dist,
        priors$beta$M,
        priors$beta$P,
        M_beta_alpha,
        P_beta_alpha,
        priors$phi[1],
        priors$phi[2],
        priors$sigma[1],
        priors$sigma[2],
        priors$mu[1],
        priors$mu[2],
        starting$beta,
        starting$gamma,
        starting$alpha,
        1 / starting$sigma^2,
        starting$hp,
        beta_alpha,
        N,
        n,
        p,
        q,
        r,
        p_alpha,
        s - 1,
        n.samples,
        n.thin,
        n.burnin,
        n.report
      )
      res$params[, p + 1] <- 1 / sqrt(res$params[, p + 1])
    } else {
      res <- spQuantileRcpp(
        quantile,
        y,
        x,
        v,
        x_alpha,
        dist,
        priors$beta$M,
        priors$beta$P,
        M_beta_alpha,
        P_beta_alpha,
        priors$phi[1],
        priors$phi[2],
        priors$sigma[1],
        priors$sigma[2],
        priors$mu[1],
        priors$mu[2],
        starting$beta,
        starting$alpha,
        1 / starting$sigma,
        starting$hp,
        beta_alpha,
        N,
        n,
        p,
        r,
        p_alpha,
        s - 1,
        n.samples,
        n.thin,
        n.burnin,
        n.report,
        parallel,
        n.threads
      )
      res$params[, p + 1] <- 1 / res$params[, p + 1]
    }
    time <- Sys.time() - time
    
    colnames(res$params) <- 1:ncol(res$params)
    if (is.null(colnames(x))) {
      colnames(res$params) <- c(paste0("V", 1:p), "sigma")
    } else {
      colnames(res$params) <- c(colnames(x), "sigma")
    }
    
    colnames(res$process) <- 1:ncol(res$process)
    for (m in 1:r) {
      if (is.null(colnames(x_alpha[[m]]))) {
        colnames(x_alpha[[m]]) <- paste0("V", 1:p_alpha[m])
      }
      idx <- (m - 1) * (n + mean(p_alpha[1:(m-1)]) + 3) + 1:(n + p_alpha[m] + 3)
      res$process[,idx[length(idx) - 1]] <- 1 / sqrt(res$process[,idx[length(idx) - 1]])
      colnames(res$process)[idx] <- 
        c(paste0("beta", m, "(s", 1:n, ")"), 
          paste0(colnames(x_alpha[[m]]), ",", m),
          paste0(c("mu", "sigma", "phi"), m))
    }
    
    res$params <- coda::mcmc(
      res$params, 
      start = n.burnin + 1, 
      end = n.burnin + n.samples, 
      thin = n.thin)
    
    res$process <- coda::mcmc(
      res$process, 
      start = n.burnin + 1, 
      end = n.burnin + n.samples, 
      thin = n.thin)
    
    z <- list(
      p.params.samples = res$params,
      p.process.samples = res$process
    )
    if (!is.null(res$missing)) {
      colnames(res$missing) <- which(!is.finite(y))
      res$missing <- coda::mcmc(
        res$missing, 
        start = n.burnin + 1, 
        end = n.burnin + n.samples, 
        thin = n.thin)
      z$Y.missing.samples <- res$missing
    }    
    z$mcmc <- list(
      n.samples = n.samples, 
      n.thin = n.thin, 
      n.burnin = n.burnin, 
      time = time)
    z$method <- method
  }
  
  if (!method.mean)
    z$quantile <- quantile
  z$rank <- p
  z$na.action <- attr(mf, "na.action")
  z$xlevels <- stats::.getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model) 
    z$model <- mf
  if (ret.x) 
    z$x <- x
  if (ret.y) 
    z$y <- y
  
  class(z) <- "spTm"
  
  return(z)
}
