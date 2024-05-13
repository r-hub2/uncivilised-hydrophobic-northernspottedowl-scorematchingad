# same parameters as cppad_closed
#' @title Iterative Score Matching Estimator Using Conjugate-Gradient Descent
#' @family generic score matching tools
#' @description 
#' Uses conjugate gradient descent to search for a vector of parameters such that gradient of the score matching discrepancy is within tolerance of zero.
#' Also estimates standard errors and covariance.
#' @inheritParams cppad_closed
#' @param theta The starting parameter set
#' @param control Control parameters passed to [`optimx::Rcgmin()`]
#' @details
#' The score matching discrepancy function and gradient of the score matching function are passed to [`optimx::Rcgmin()`]. 
#' The call to [`optimx::Rcgmin()`] uses the *sum* of observations (as opposed to the mean) to reduce floating point inaccuracies. This has implications for the meaning of the control parameters passed to `Rcgmin()` (e.g. `tol`). The results are converted into averages so the use of sums can be ignored when not setting control parameters, or studying the behaviour of Rcgmin. 
#'
#' Standard errors use the Godambe information matrix (aka sandwich method) and are only computed when the weights are constant.
#' The estimate of the sensitivity matrix \eqn{G} is
#' the negative of the average over the Hessian of `smdtape` evaluated at each observation in `Y`.
# \deqn{\hat{G(\theta)} = \hat{E} -H(smd(\theta;Y))),}
# where \eqn{smd} is the score matching discrepancy function represented by `smdtape`,
# \eqn{H} is the Hessian with respect to \eqn{\theta}, which is constant for quadratic-form functions,
# 
#' The estimate of the variability matrix \eqn{J} is then
#' the sample covariance (denominator of \eqn{n-1}) of the gradiant of `smdtape` evaluated at each of the observations in `Y` for the estimated \eqn{\theta}.
# \deqn{\hat{J}(\theta) = var(grad(w smd(\theta;Y))),}

#' The variance of the estimator is then estimated as
#' \eqn{G^{-1}JG^{-1}/n,}
# \deqn{\hat{G}(\theta)^{-1}\hat{J}(\theta)\hat{G}(\theta)^{-1}/n,}
#' where `n` is the number of observations.
#'
#' Taylor approximation is available because boundary weight functions and transformations of the measure in Hyvärinen divergence can remove singularities in the model log-likelihood, however evaluation at these singularities may still involve computing intermediate values that are unbounded.
#' If the singularity is ultimately removed, then Taylor approximation from a nearby location will give a very accurate evaluation at the removed singularity.
#'
#' # Warning
#' There appears to be floating point issues with evaluation of the gradient leading to slow or no convergence in many situations. Tweaking the convergence tolerance `tol` may be appropriate. If the closed form solution exists please use `cppad_closed()` instead.
#' @examples
#' smdtape <- buildsmdtape("sim", "sqrt", "sph", "ppi",
#'               ytape = rep(1/3, 3),
#'               usertheta = ppi_paramvec(p = 3),
#'               bdryw = "minsq", acut = 0.01,
#'               verbose = FALSE
#'               )$smdtape
#' Y <- rppi_egmodel(100)
#' cppad_search(smdtape, 0.9 * Y$theta, Y$sample)
#' sum((smvalues_wsum(smdtape, Y$sample, Y$theta)$grad/nrow(Y$sample))^2)
#' @export
cppad_search <- function(smdtape, theta, Y, Yapproxcentres = NA * Y, w = rep(1, nrow(Y)), approxorder = 10, control = list(tol = 1E-15, checkgrad = TRUE)){
  Jsmdfun <- tapeJacobian(smdtape)
  Hsmdfun <- tapeJacobian(Jsmdfun)
  
  smdfun_u <- tapeSwap(smdtape) #don't use a boundary point for taping here!
  Jsmdfun_u <- tapeSwap(Jsmdfun)

  smdobj <- function(atheta){
    evaltape_wsum(smdfun_u, xmat = Y, pmat = atheta, w = w, xcentres = Yapproxcentres, approxorder = approxorder)
  }
  smdgrad <- function(atheta){
    evaltape_wsum(Jsmdfun_u, xmat = Y, pmat = atheta, w = w, xcentres = Yapproxcentres, approxorder = approxorder)
  }

  # useful to debugging as Rcgmin hides the error codes
  #  and relatively cheap:
  # evaluating above functions at the start point
  stopifnot(is.finite(smdobj(theta)))
  stopifnot(all(is.finite(smdgrad(theta))))
  # now do the search 
  out <- optimx::Rcgmin(par = theta,
                        fn = smdobj,
                        gr = smdgrad,
                        control = control)
  if (out$convergence == 2){
    if (grepl("Initial point", out$message)){
      stop(paste(out$message, "Initial point was", paste(theta, collapse = " ")))
    } else {
      stop(paste(out$message, "Perhaps evaltape_wsum() generates a non-number?"))
    }
  }
  if (out$convergence != 0){warning("Optimisation did not converge.")}

  # return results as if averages, not sums were used
  out$value <- out$value / sum(w)

  out$SE <- "Not calculated."
  if (isTRUE(all(w[[1]] == w))){
    out$SE <- sqrt(diag(sme_estvar(smdtape, estimate = out$par, Y = Y, Yapproxcentres = Yapproxcentres, approxorder = approxorder)))
  }
  gradatest <- smdgrad(out$par) / sum(w)
  out$sqgradsize <- sum(gradatest^2)
  out$est <- out$par
  out$par <- NULL
  return(out)
}
#' @return A list of
#' + `est` The estimate
#' + `sqgradsize` The squared size of the gradient at the estimate
#' + `SE` The standard error of the estimates when `w` is constant. If `w` changes between observations then this slot contains the string "Not calculated."
#' + additional output from [`optimx::Rcgmin()`]

