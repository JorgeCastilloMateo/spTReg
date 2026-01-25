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
#'   \deqn{\mathbf{Y}_{t\ell} = \mathbf{O}_{t\ell} + \mathbf{P} (\mathbf{Y}_{t,\ell-1} - \mathbf{O}_{t,\ell-1}) + \mathbf{w}_{t\ell} +\bm{\epsilon}_{t\ell},}
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
  x_gamma, # list of design matrix for each tp process
  x_alpha, # list of design matrix for each sp process
  w.bool = FALSE, # daily dynamic spatial process
  data, 
  subset,
  na.action,
  method = c("mean", "quantile"),
  quantile = 0.5,
  coords, 
  dates,
  priors = list(
    "beta"    = c(0, 1e-08), 
    "sigma"   = c(0.1, 0.1),
    "tp" = list(
      "beta"  = c(0, 1e-08),
      "rho"   = c(0, 1e-08),
      "sigma" = c(0.1, 0.1)
    ),
    "sp" = list(
      "beta"  = c(0, 1e-08),
      "sigma" = c(0.1, 0.1),
      "phi"   = c(2, 100)
    ),
    "st" = list(
      "rho"   = c(0, 1e-08),
      "sigma" = c(0.1, 0.1),
      "phi"   = c(2, 100)
    )
  ),
  starting = list(
    "beta"    = 0.01, 
    "sigma"   = 1,
    "tp" = list(
      "tp"    = 0,
      "beta"  = 0,
      "rho"   = 0,
      "sigma" = 1
    ),
    "sp" = list(
      "sp"    = 0,
      "beta"  = 0,
      "sigma" = 1,
      "phi"   = 0.01
    ),
    "st" = list(
      "st"    = 0,
      "rho"   = 0,
      "sigma" = 1,
      "phi"   = 0.01
    )
  ),
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
      p_tp <- 0
      x_gamma <- list()
      u <- matrix(nrow = N, ncol = 0)
    } else {
      q <- ncol(u) #sum(attr(mt, "dataClasses") == "tpCoef")
      #stopifnot(is.list(x_gamma))
      #if (length(x_gamma) != q) 
      #  stop("'x_gamma' length != q")
      #if (!all(sapply(x_gamma, is.matrix, simplify = TRUE)))
      #  stop("'x_gamma' elements are not 'matrix'")
      #p_tp <- sapply(x_gamma, ncol, simplify = TRUE)
      p_tp <- 0
    }
    if (missing(v) || is.null(v)) {
      r <- 0
      p_sp <- 0
      x_alpha <- list()
      v <- matrix(nrow = N, ncol = 0)
    } else {
      r <- ncol(v) #sum(attr(mt, "dataClasses") == "spCoef")
      stopifnot(is.list(x_alpha))
      if (length(x_alpha) != r) 
        stop("'x_alpha' length != r")
      if (!all(sapply(x_alpha, is.matrix, simplify = TRUE)))
        stop("'x_alpha' elements are not 'matrix'")
      p_sp <- sapply(x_alpha, ncol, simplify = TRUE)
    }
    
    # date and site
    date_site <- data.frame(
      "date" = dates, 
      "site" = rep(1:n, each = N / n))
    date_site <- .group_dates(date_site)
    
    T <- length(unique(date_site$block))
    L <- tapply(date_site$num_within_block, date_site$block, 
                function(x) length(unique(x)))
    
    # check priors and starting
    priors   <- .check_priors(    priors, p, q, p_tp, r, p_sp, N, T, n, w.bool)
    starting <- .check_starting(starting, p, q, p_tp, r, p_sp, N, T, n, w.bool)
    
    # dist
    dist <- dist1(coords)
    
    if (!verbose) 
      n.report <- n.samples + 1
    
    # model fitting
    time <- Sys.time()
    if (method.mean) {
      res <- spMeanRcpp(
        y,
        x,
        u,
        v,
        x_gamma,
        x_alpha,
        dist,
        priors$beta$M,            ### priors
        priors$beta$P,
        priors$sigma,
        priors$tp$beta$M,
        priors$tp$beta$P,
        priors$tp$rho,              # 2 x q
        priors$tp$sigma,            # 2 x q
        priors$sp$beta$M,
        priors$sp$beta$P,
        priors$sp$sigma,            # 2 x r
        priors$sp$phi,              # 2 x r
        priors$st$rho,
        priors$st$sigma,
        priors$st$phi,
        starting$beta,            ### starting
        1 / starting$sigma^2,
        starting$tp$tp,             # n x q
        starting$tp$beta,           # list lengths p_tp
        rbind(starting$tp$rho,
          1 / starting$tp$sigma^2), # 2 x q
        starting$sp$sp,             # n x r
        starting$sp$beta,           # list lengths p_sp
        rbind(
          1 / starting$sp$sigma^2, 
          starting$sp$phi),         # 2 x r
        starting$st$st,             # N (vec)
        c(starting$st$rho, 
          1 / starting$st$sigma^2, 
          starting$st$phi),         # hp_w : rho_w prec_w phi_w
        N,
        n,
        T,
        L,
        p,
        q,
        r,
        p_tp,
        p_sp,
        w.bool,
        date_site$site - 1,             # s
        date_site$block - 1,            # t
        date_site$num_within_block - 1, # l
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
        u,
        v,
        x_gamma,
        x_alpha,
        dist,
        priors$beta$M,            ### priors
        priors$beta$P,
        priors$sigma,
        priors$tp$beta$M,
        priors$tp$beta$P,
        priors$tp$rho,              # 2 x q
        priors$tp$sigma,            # 2 x q
        priors$sp$beta$M,
        priors$sp$beta$P,
        priors$sp$sigma,            # 2 x r
        priors$sp$phi,              # 2 x r
        priors$st$rho,
        priors$st$sigma,
        priors$st$phi,
        starting$beta,            ### starting
        1 / starting$sigma^2,
        starting$tp$tp,             # n x q
        starting$tp$beta,           # list lengths p_tp
        rbind(starting$tp$rho,
          1 / starting$tp$sigma^2), # 2 x q
        starting$sp$sp,             # n x r
        starting$sp$beta,           # list lengths p_sp
        rbind(
          1 / starting$sp$sigma^2, 
          starting$sp$phi),         # 2 x r
        starting$st$st,             # N (vec)
        c(starting$st$rho, 
          1 / starting$st$sigma^2, 
          starting$st$phi),         # hp_w : rho_w prec_w phi_w
        N,
        n,
        T,
        L,
        p,
        q,
        r,
        p_tp,
        p_sp,
        w.bool,
        date_site$site - 1,             # s
        date_site$block - 1,            # t
        date_site$num_within_block - 1, # l
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
    
    res$params <- coda::mcmc(
      res$params, 
      start = n.burnin + 1, 
      end = n.burnin + n.samples, 
      thin = n.thin)
    z <- list(p.params.samples = res$params)
    
    if (!is.null(res$sp)) {
      colnames(res$sp) <- 1:ncol(res$sp)
      for (m in 1:r) {
        if (is.null(colnames(x_alpha[[m]]))) {
          colnames(x_alpha[[m]]) <- paste0("V", 1:p_sp[m])
        }
        idx <- (m - 1) * (n + mean(p_sp[1:(m-1)]) + 2) + 1:(n + p_sp[m] + 2)
        res$sp[,idx[length(idx) - 1]] <- 1 / sqrt(res$sp[,idx[length(idx) - 1]])
        colnames(res$sp)[idx] <- 
          c(paste0("beta", m, "(s", 1:n, ")"), 
            paste0(colnames(x_alpha[[m]]), ",", m),
            paste0(c("sigma", "phi"), m))
      }
      
      res$sp <- coda::mcmc(
        res$sp, 
        start = n.burnin + 1, 
        end = n.burnin + n.samples, 
        thin = n.thin)
      z$p.sp.samples <- res$sp
    }
    
    if (!is.null(res$st)) {
      colnames(res$st) <- 1:ncol(res$st)
      colnames(res$st)[1:N] <- date_site$name
      colnames(res$st)[N + 1:3] <- c("rho", "sigma", "phi")
      res$st[, N + 2] <- 1 / sqrt(res$st[, N + 2])
      
      res$st <- coda::mcmc(
        res$st, 
        start = n.burnin + 1, 
        end = n.burnin + n.samples, 
        thin = n.thin)
      z$p.st.samples <- res$st
    }
    
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
  z$rank <- c("p" = p)
  if (r > 1) 
    z$rank <- c(z$rank, "r" = r, "p_sp" = p_sp)
  if (w.bool) 
    z$rank <- c(z$rank, "w" = 1)
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

.check_priors <- function(
    priors, p, 
    q, p_tp, r, p_sp, N, T, n, w.bool) {
  
  stopifnot(is.list(priors))
  
  # beta
  if (!"beta" %in% names(priors)) {
    warning("'priors' (list) should contain 'beta'. Default = c(0, 1e-08)")
    priors$beta <- c(0, 1e-08)
  } else if (is.list(priors$beta)) {
    if (!all(c("M", "P") %in% names(priors$beta))) {
      warning("'priors$beta' (list) should contain 'M' and 'P'. Default = c(0, 1e-08)")
      priors$beta <- list(
        "M" = rep(0, p),
        "P" = 1e-08 * diag(p))
    } else {
      if (length(priors$beta$M) != p)
        stop("'priors$beta$M' (vec) length != p")
      if (!all(dim(priors$beta$P) == p))
        stop("'priors$beta$P' (mat) dim != p x p")
    }
  } else {
    if (length(priors$beta) != 2)
      stop("'priors$beta' (vec) length != 2")
    priors$beta <- list("M" = rep(priors$beta[1], p),
                        "P" = priors$beta[2] * diag(p))
  }
  
  # sigma
  if (!"sigma" %in% names(priors)) {
    warning("'priors' (list) should contain 'sigma'. Default = c(0.1, 0.1)")
    priors$sigma <- c(0.1, 0.1)
  } else if (length(priors$sigma) != 2)
    stop("'priors$sigma' (vec) length != 2")
  
  # tp
  if (!missing(q) && q > 0) {
    
    if (!"tp" %in% names(priors)) {
      warning(
        "'priors' (list) should contain 'tp'.
        \n Default = list('beta' = c(0, 1e-08), 'rho' = c(0, 1e-08), 'sigma' = c(0.1, 0.1))")
      priors$tp <- list(
        "beta" = c(0, 1e-08),
        "rho" = c(0, 1e-08),
        "sigma" = c(0.1, 0.1)
      )
    }
    
    if (sum(p_tp) > 0) {
      if (!"beta" %in% names(priors$tp)) {
        warning("'priors$tp' (list) should contain 'beta'. Default = c(0, 1e-08)")
        priors$tp$beta <- list(
          "M" = sapply(p_tp, rep, x = 0, simplify = FALSE),
          "P" = sapply(p_tp, diag, x = 1e-08, simplify = FALSE))
      } else if (is.list(priors$tp$beta)) {
        if (!all(c("M", "P") %in% names(priors$tp$beta))) {
          warning("'priors$tp$beta' (list) should contain 'M' and 'P'. Default = c(0, 1e-08)")
          priors$tp$beta <- list(
            "M" = sapply(p_tp, rep, x = 0, simplify = FALSE),
            "P" = sapply(p_tp, diag, x = 1e-08, simplify = FALSE))
        } else {
          if (!all(lengths(priors$tp$beta$M) == p_tp))
            stop("'priors$tp$beta$M' (list) lengts != p_tp")
          if (!all(lengths(priors$tp$beta$P) == p_tp^2))
            stop("'priors$tp$beta$P' (list) dims != p_tp x p_tp")
        }
      } else {
        if (length(priors$tp$beta) != 2)
          stop("'priors$tp$beta' (vec) length != 2")
        priors$tp$beta <- list(
          "M" = sapply(p_tp, rep, x = priors$tp$beta[1], simplify = FALSE),
          "P" = sapply(p_tp, diag, x = priors$tp$beta[2], simplify = FALSE))
      }
    } else {
      priors$tp$beta <- list("M" = list(), "P" = list())
    }
    
    if (!"rho" %in% names(priors$tp)) {
      warning("'priors$tp' (list) should contain 'rho'. Default = c(0, 1e-08)")
      priors$tp$rho <- matrix(c(0, 1e-08), nrow = 2, ncol = q)
    } else if (is.matrix(priors$tp$rho) && 
               !all(dim(priors$tp$rho) == c(2, q))) {
      stop("'starting$tp$rho' (mat) dim != c(2, q)")
    } else if (length(priors$tp$rho) == 2) {
      priors$tp$rho <- matrix(priors$tp$rho, nrow = 2, ncol = q)
    } else
      stop("'priors$tp$rho' (vec) length != 2")
    
    if (!"sigma" %in% names(priors$sp)) {
      warning("'priors$tp' (list) should contain 'sigma'. Default = c(0.1, 0.1)")
      priors$tp$sigma <- matrix(0.1, nrow = 2, ncol = q)
    } else if (is.matrix(priors$tp$sigma) && 
               !all(dim(priors$tp$sigma) == c(2, q))) {
      stop("'starting$tp$sigma' (mat) dim != c(2, q)")
    } else if (length(priors$tp$sigma) == 2) {
      priors$tp$sigma <- matrix(priors$tp$sigma, nrow = 2, ncol = q)
    } else
      stop("'priors$tp$sigma' (vec) length != 2")
    
  } else if (!missing(q)) {
    priors$tp <- list(
      "beta" = list("M" = list(), "P" = list()),
      "rho" = matrix(),
      "sigma" = matrix()
    )
  }
  
  # sp
  if (!missing(r) && r > 0) {
    
    if (!"sp" %in% names(priors)) {
      warning(
        "'priors' (list) should contain 'sp'.
        \n Default = list('beta' = c(0, 1e-08), 'sigma' = c(0.1, 0.1), 'phi' = c(2, 100))")
      priors$sp <- list(
        "beta" = c(0, 1e-08),
        "sigma" = c(0.1, 0.1),
        "phi" = c(2, 100))
    }
    
    if (sum(p_sp) > 0) {
      if (!"beta" %in% names(priors$sp)) {
        warning("'priors$sp' (list) should contain 'beta'. Default = c(0, 1e-08)")
        priors$sp$beta <- list(
          "M" = sapply(p_sp, rep, x = 0, simplify = FALSE),
          "P" = sapply(p_sp, diag, x = 1e-08, simplify = FALSE))
      } else if (is.list(priors$sp$beta)) {
        if (!all(c("M", "P") %in% names(priors$sp$beta))) {
          warning("'priors$sp$beta' (list) should contain 'M' and 'P'. Default = c(0, 1e-08)")
          priors$sp$beta <- list(
            "M" = sapply(p_sp, rep, x = 0, simplify = FALSE),
            "P" = sapply(p_sp, diag, x = 1e-08, simplify = FALSE))
        } else {
          if (!all(lengths(priors$sp$beta$M) == p_sp))
            stop("'priors$sp$beta$M' (list) lengts != p_sp")
          if (!all(lengths(priors$sp$beta$P) == p_sp^2))
            stop("'priors$sp$beta$P' (list) dims != p_sp x p_sp")
        }
      } else {
        if (length(priors$sp$beta) != 2)
          stop("'priors$sp$beta' (vec) length != 2")
        priors$sp$beta <- list(
          "M" = sapply(p_sp, rep, x = priors$sp$beta[1], simplify = FALSE),
          "P" = sapply(p_sp, diag, x = priors$sp$beta[2], simplify = FALSE))
      }
    } else {
      priors$sp$beta <- list("M" = list(), "P" = list())
    }
    
    if (!"sigma" %in% names(priors$sp)) {
      warning("'priors$sp' (list) should contain 'sigma'. Default = c(0.1, 0.1)")
      priors$sp$sigma <- matrix(0.1, nrow = 2, ncol = r)
    } else if (is.matrix(priors$sp$sigma) && 
               !all(dim(priors$sp$sigma) == c(2, r))) {
      stop("'starting$sp$sigma' (mat) dim != c(2, r)")
    } else if (length(priors$sp$sigma) == 2) {
      priors$sp$sigma <- matrix(priors$sp$sigma, nrow = 2, ncol = r)
    } else
      stop("'priors$sp$sigma' (vec) length != 2")
    
    if (!"phi" %in% names(priors$sp)) {
      warning("'priors$sp' (list) should contain 'phi'. Default = c(2, 100)")
      priors$sp$phi <- matrix(c(2, 100), nrow = 2, ncol = r)
    } else if (is.matrix(priors$sp$phi) && 
               !all(dim(priors$sp$phi) == c(2, r))) {
      stop("'starting$sp$phi' (mat) dim != c(2, r)")
    } else if (length(priors$sp$phi) == 2) {
      priors$sp$phi <- matrix(priors$sp$phi, nrow = 2, ncol = r)
    } else
      stop("'priors$sp$phi' (vec) length != 2")
    
  } else if (!missing(r)) {
    priors$sp <- list(
      "beta" = list("M" = list(), "P" = list()),
      "sigma" = matrix(),
      "phi" = matrix()
    )
  }
  
  # st
  if (!missing(w.bool) && w.bool) {
    
    if (!"rho" %in% names(priors$st)) {
      warning("'priors$st' (list) should contain 'rho'. Default = c(0, 1e-08)")
      priors$st$rho <- c(0, 1e-08)
    } else if (length(priors$st$rho) != 2)
      stop("'priors$st$rho' (vec) length != 2")
    
    if (!"sigma" %in% names(priors$st)) {
      warning("'priors$st' (list) should contain 'sigma'. Default = c(0.1, 0.1)")
      priors$st$sigma <- c(0.1, 0.1)
    } else if (length(priors$st$sigma) != 2)
      stop("'priors$st$sigma' (vec) length != 2")
    
    if (!"phi" %in% names(priors$st)) {
      warning("'priors$st' (list) should contain 'phi'. Default = c(2, 100)")
      priors$st$phi <- c(2, 100)
    } else if (length(priors$st$phi) != 2)
      stop("'priors$st$phi' (vec) length != 2")
    
  } else if (!missing(w.bool)) {
    priors$st <- list(
      "rho" = NA,
      "sigma" = NA,
      "phi" = NA
    )
  }
  
  priors
}

.check_starting <- function(
    starting, p, 
    q, p_tp, r, p_sp, N, T, n, w.bool) {
  
  stopifnot(is.list(starting))
  
  # beta
  if (!"beta" %in% names(starting)) {
    warning("'starting' should contain 'beta'. Default = 0.01")
    starting$beta <- rep(0.01, p)
  } else if (!length(starting$beta) %in% c(1, p)) {
    stop("'starting$beta' length != 1 or p")
  } else if (length(starting$beta) == 1)
    starting$beta <- rep(starting$beta, p)
  
  # sigma
  if (!"sigma" %in% names(starting)) {
    warning("'starting' should contain 'sigma'. Default = 1")
    starting$sigma <- 1
  } else if (length(starting$sigma) != 1)
    stop("'starting$sigma' length != 1")
  
  # tp
  if (!missing(q) && q > 0) {
    if (!"tp" %in% names(starting)) {
      warning(
        "'starting' should contain 'tp'.
        \n Default = list('tp' = 0, 'beta' = 0, 'rho' = 0, 'sigma' = 1)")
      starting$tp <- list(
        "tp" = 0,
        "beta" = 0,
        "rho" = 0,
        "sigma" = 1)
    }
    
    if (!"tp" %in% names(starting$tp)) {
      warning("'starting$tp' should contain 'tp'. Default = 0")
      starting$tp$tp <- matrix(0, nrow = T, ncol = q)
    } else if (!is.matrix(starting$tp$tp) && 
               !is.array(starting$tp$tp) && 
               length(starting$tp$tp) != 1) {
      stop("'starting$tp$tp' is numeric vector length != 1")
    } else if (is.matrix(starting$tp$tp) && 
               !all(dim(starting$tp$tp) == c(T, q))) {
      stop("'starting$tp$tp' is matrix dim != c(T, q)")
    } else if (length(starting$tp$tp) == 1) 
      starting$tp$tp <- matrix(starting$tp$tp, nrow = T, ncol = q)
    
    if (sum(p_tp) > 0) {
      if (!"beta" %in% names(starting$tp)) {
        warning("'starting$tp' should contain 'beta'. Default = 0")
        starting$tp$beta <- sapply(p_tp, rep, x = 0, simplify = FALSE)
      } else if (!is.matrix(starting$tp$beta) && 
                 !is.array(starting$tp$beta) && 
                 !is.list(starting$tp$beta) &&
                 length(starting$tp$beta) != 1) {
        stop("'starting$tp$beta' is numeric vector length != 1")
      } else if (is.list(starting$tp$beta) && 
                 !all(lengths(starting$tp$beta) == p_tp)) {
        stop("'starting$tp$beta' is list lengts != p_tp")
      } else if (!is.list(starting$tp$beta) && 
                 length(starting$tp$beta) == 1)
        starting$tp$beta <- sapply(p_tp, rep, x = starting$tp$beta, simplify = FALSE)
    } else {
      starting$tp$beta <- list()
    }
    
    if (!"rho" %in% names(starting$tp)) {
      warning("'starting$tp' should contain 'rho'. Default = 0")
      starting$tp$rho <- rep(0, q)
    } else if (!length(starting$tp$rho) %in% c(1, q)) {
      stop("'starting$tp$rho' length != 1 or q")
    } else if (length(starting$tp$rho) == 1)
      starting$tp$rho <- rep(starting$tp$rho, q)
    
    if (!"sigma" %in% names(starting$tp)) {
      warning("'starting$tp' should contain 'sigma'. Default = 1")
      starting$tp$sigma <- rep(1, q)
    } else if (!length(starting$tp$sigma) %in% c(1, q)) {
      stop("'starting$tp$sigma' length != 1 or q")
    } else if (length(starting$tp$sigma) == 1)
      starting$tp$sigma <- rep(starting$tp$sigma, q)
  } else if (!missing(q)) {
    starting$tp <- list(
      "tp" = matrix(nrow = 0, ncol = 0),
      "beta" = list(),
      "rho" = NA,
      "sigma" = NA
    )
  }
  
  # sp
  if (!missing(r) && r > 0) {
    if (!"sp" %in% names(starting)) {
      warning(
        "'starting' should contain 'sp'.
        \n Default = list('sp' = 0, 'beta' = 0, 'sigma' = 1, 'phi' = 0.01)")
      starting$sp <- list(
        "sp" = 0,
        "beta" = 0,
        "sigma" = 1,
        "phi" = 0.01)
    }
    
    if (!"sp" %in% names(starting$sp)) {
      warning("'starting$sp' should contain 'sp'. Default = 0")
      starting$sp$sp <- matrix(0, nrow = n, ncol = r)
    } else if (!is.matrix(starting$sp$sp) && 
               !is.array(starting$sp$sp) && 
               length(starting$sp$sp) != 1) {
      stop("'starting$sp$sp' is numeric vector length != 1")
    } else if (is.matrix(starting$sp$sp) && 
               !all(dim(starting$sp$sp) == c(n, r))) {
      stop("'starting$sp$sp' is matrix dim != c(n, r)")
    } else if (length(starting$sp$sp) == 1) 
      starting$sp$sp <- matrix(starting$sp$sp, nrow = n, ncol = r)
    
    if (sum(p_sp) > 0) {
      if (!"beta" %in% names(starting$sp)) {
        warning("'starting$sp' should contain 'beta'. Default = 0")
        starting$sp$beta <- sapply(p_sp, rep, x = 0, simplify = FALSE)
      } else if (!is.matrix(starting$sp$beta) && 
                 !is.array(starting$sp$beta) && 
                 !is.list(starting$sp$beta) &&
                 length(starting$sp$beta) != 1) {
        stop("'starting$sp$beta' is numeric vector length != 1")
      } else if (is.list(starting$sp$beta) && 
                 !all(lengths(starting$sp$beta) == p_sp)) {
        stop("'starting$sp$beta' is list lengts != p_sp")
      } else if (!is.list(starting$sp$beta) && 
                 length(starting$sp$beta) == 1)
        starting$sp$beta <- sapply(p_sp, rep, x = starting$sp$beta, simplify = FALSE)
    } else {
      starting$sp$beta <- list()
    }
    
    if (!"sigma" %in% names(starting$sp)) {
      warning("'starting$sp' should contain 'sigma'. Default = 1")
      starting$sp$sigma <- rep(1, r)
    } else if (!length(starting$sp$sigma) %in% c(1, r)) {
      stop("'starting$sp$sigma' length != 1 or r")
    } else if (length(starting$sp$sigma) == 1)
      starting$sp$sigma <- rep(starting$sp$sigma, r)
    
    if (!"phi" %in% names(starting$sp)) {
      warning("'starting$sp' should contain 'phi'. Default = 0.01")
      starting$sp$phi <- rep(0.01, r)
    } else if (!length(starting$sp$phi) %in% c(1, r)) {
      stop("'starting$sp$phi' length != 1 or r")
    } else if (length(starting$sp$phi) == 1)
      starting$sp$phi <- rep(starting$sp$phi, r)
  } else if (!missing(r)) {
    starting$sp <- list(
      "sp" = matrix(nrow = 0, ncol = 0),
      "beta" = list(),
      "sigma" = NA,
      "phi" = NA
    )
  }
  
  # st
  if (!missing(w.bool) && w.bool) {
    if (!"st" %in% names(starting)) {
      warning(
        "'starting' should contain 'st'.
        \n Default = list('st' = 0, 'rho' = 0, 'sigma' = 1, 'phi' = 0.01)")
      starting$st <- list(
        "st" = 0,
        "rho" = 0,
        "sigma" = 1,
        "phi" = 0.01)
    }
    if (!"st" %in% names(starting$st)) {
      warning("'starting$st' should contain 'st'. Default = 0")
      starting$st$st <- rep(0, N)
    } else if (!length(starting$st$st) %in% c(1, N)) {
      stop("'starting$st$st' length != 1 or N")
    } else if (length(starting$st$st) == 1) 
      starting$st$st <- rep(starting$st$st, N)
    if (!"rho" %in% names(starting$st)) {
      warning("'starting$st' should contain 'rho'. Default = 0")
      starting$st$rho <- 0
    } else if (length(starting$st$rho) != 1)
      stop("'starting$st$rho' length != 1")
    if (!"sigma" %in% names(starting$st)) {
      warning("'starting$st' should contain 'sigma'. Default = 1")
      starting$st$sigma <- 1
    } else if (length(starting$st$sigma) != 1)
      stop("'starting$st$sigma' length != 1")
    if (!"phi" %in% names(starting$st)) {
      warning("'starting$st' should contain 'phi'. Default = 0.01")
      starting$st$phi <- 0.01
    } else if (length(starting$st$phi) != 1)
      stop("'starting$st$phi' length != 1")
  } else if (!missing(w.bool)) {
    starting$st <- list(
      "st" = matrix(nrow = 0, ncol = 0),
      "rho" = NA,
      "sigma" = NA,
      "phi" = NA
    )
  }
  
  starting
}

.group_dates <- function(df) {
  
  df$date <- as.Date(df$date)
  
  df$block <- NA
  df$num_within_block <- NA
  
  site <- unique(df$site)
  
  structures <- list()
  
  for (s in site) {
    idx <- which(df$site == s)
    dates_site <- df$date[idx]
    
    ord <- order(dates_site)
    dates_site <- dates_site[ord]
    
    diff_days <- c(0, diff(dates_site))
    block <- 1 + cumsum(diff_days > 1)
    
    num_within <- ave(dates_site, block, FUN = seq_along)
    
    df$block[idx[ord]] <- block
    df$num_within_block[idx[ord]] <- num_within
    
    structures[[s]] <- paste(block, num_within, sep="-", collapse=",")
  }
  
  if (length(unique(structures)) > 1) {
    stop("'dates' differ between locations")
  }
  
  df$name <- paste0(
    "t", df$block,
    "l", df$num_within_block,
    "s", df$site
  )
  
  return(df)
}
