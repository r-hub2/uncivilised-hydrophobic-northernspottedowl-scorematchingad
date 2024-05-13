# estimates of the variance of the estimator for cppad_search
sme_estvar <- function(smdfun, estimate, Y, Yapproxcentres = NA * Y, approxorder = 10){
  stopifnot(inherits(smdfun, "ADFun"))
  # generate tapes with respect to the measurement (parameter is dynamic)
  p <- ncol(Y)
  Jsmdfun <- tapeJacobian(smdfun)
  Hsmdfun <- tapeJacobian(Jsmdfun)
  
  Jsmdfun_u <- tapeSwap(Jsmdfun)
  Hsmdfun_u <- tapeSwap(Hsmdfun)

  hess <- evaltape_wsum(Hsmdfun_u, xmat = Y, pmat = 0*estimate, xcentres = Yapproxcentres, approxorder = 10)
  sens <- -matrix(hess, byrow = TRUE, ncol = length(estimate))/nrow(Y)
  sensinv <- solve(sens)

  grads <- evaltape(Jsmdfun_u, xmat = Y, pmat = estimate, xcentres = Yapproxcentres, approxorder = 10)
  variability <- stats::cov(grads)

  Ginfinv <- sensinv %*% variability %*% sensinv /nrow(Y) #inverse of the Godambe information matrix, also called the sandwich information matrix
  return(Ginfinv) 
}

