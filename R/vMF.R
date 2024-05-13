#' @title Score Matching Estimator for the von-Mises Fisher Distribution
#' @family directional model estimators
#' @description
#' In general the normalising constant in von Mises Fisher distributions is hard to compute, so \insertCite{mardia2016sc;textual}{scorematchingad} suggested a hybrid method that uses maximum likelihood to estimate the mean direction and score matching for the concentration.
#' We can also estimate all parameters using score matching (`smfull` method), although this estimator is likely to be less efficient than the hybrid estimator.
#' On the circle the hybrid estimators were often nearly as efficient as maximum likelihood estimators \insertCite{mardia2016sc}{scorematchingad}.
#' For maximum likelihood estimators of the von Mises Fisher distribution, which all use approximations of the normalising constant, consider [`movMF::movMF()`].

#' @section von Mises Fisher Model: 
#' The von Mises Fisher density is proportional to
#' \deqn{\exp(\kappa \mu^T z),}
#' where \eqn{z} is on a unit sphere,
#' \eqn{\kappa} is termed the *concentration*,
#' and \eqn{\mu} is the *mean direction unit vector*.
#' The effect of the \eqn{\mu} and \eqn{\kappa} can be decoupled in a sense \insertCite{@p169, @mardia2000di}{scorematchingad}, allowing for estimating \eqn{\mu} and \eqn{\kappa} separately.

#' @details
#' The full score matching estimator (`method = "smfull"`) estimates \eqn{\kappa \mu}.
#' The hybrid estimator (`method = "Mardia"`) estimates \eqn{\kappa} and \eqn{\mu} separately.
#' Both use [`cppad_closed()`] for score matching estimation.
#' @param Y A matrix of multivariate observations in Cartesian coordinates. Each row is a multivariate measurement (i.e. each row corresponds to an individual).
#' @param paramvec `smfull` method only: Optional. A vector of same length as the dimension, representing the elements of the \eqn{\kappa \mu} vector. 
#' @param method Either "Mardia" or "hybrid" for the hybrid score matching estimator from \insertCite{mardia2016sc;textual}{scorematchingad}
#'  or "smfull" for the full score matching estimator.
#' @param w An optional vector of weights for each measurement in `Y`
#' @references
#' \insertAllCited{}
#' @examples
#' if (requireNamespace("movMF")){
#'   Y <- movMF::rmovMF(1000, 100 * c(1, 1) / sqrt(2))
#'   movMF::movMF(Y, 1) #maximum likelihood estimate
#' } else {
#'   Y <- matrix(rnorm(1000 * 2, sd = 0.01), ncol = 2)
#'   Y <- Y / sqrt(rowSums(Y^2))
#' }
#' vMF(Y, method = "smfull")
#' vMF(Y, method = "Mardia")
#' vMF(Y, method = "hybrid")
#' @return
#' A list of `est`, `SE` and `info`.
#'  * `est` contains the estimates in vector form, `paramvec`, and with user friendly names `k` and `m`.
#'  * `SE` contains estimates of the standard errors if computed. See [`cppad_closed()`].
#'  * `info` contains a variety of information about the model fitting procedure and results.
#' @export
vMF <- function(Y, paramvec = NULL, method = "Mardia", w = rep(1, nrow(Y))){
  fit <- NULL
  if (method == "smfull"){
    if (is.null(paramvec)){
      paramvec <- rep(NA, ncol(Y))
    }
    fit <- vMF_full(Y, paramvec, w=w)
  }
  if (method %in% c("Mardia", "mardia", "hybrid", "Hybrid")){
    if (!is.null(paramvec)){if (any(!is.na(paramvec))){stop("Mardia estimator cannot fix any elements of paramvec")}}
    fit <- vMF_Mardia(Y, w=w)
  }
  if (is.null(fit)){stop(sprintf("Method '%s' is not valid", method))}

  return(fit)
}


#for vMF_Mardia startk must be the value of the k parameter
vMF_Mardia <- function(sample, w = rep(1, nrow(sample))){
  mu <- vMF_m(sample, w = w)
  samplestd <- vMF_stdY(sample, m = mu, w = w)
  # check: mustd <- colMeans(samplestd); mustd <- mustd / sqrt(sum(mustd^2))
  kappaest <- vMF_kappa(Y = samplestd, w = w) 
  return(list(
    est = list(paramvec = kappaest$k * mu,
               k = kappaest$k,
               m = mu),
    SE = list(paramvec = "Not calculated.",
              k = kappaest$SE,
              m = "Not calculated."),
    info = kappaest
  ))
}


vMF_full <- function(sample, usertheta, w = NULL){
  p <- ncol(sample)
  stopifnot(length(usertheta) == p)

  tapes <- buildsmdtape("sph","identity", "sph", "vMF",
                        rep(1, p)/sqrt(p), 
                        usertheta = usertheta,
                        bdryw = "ones",
                        verbose = FALSE)
  out <- cppad_closed(tapes$smdtape, Y = sample, w=w)
  theta <- t_fu2t(out$est, usertheta)

  if (isa(out$SE, "numeric")){
    SE <- t_fu2t(out$SE, 0 * usertheta)
    SE <- list(paramvec = SE)
  } else {SE <- out$SE}
  return(list(
    est = list(paramvec = theta,
               k = sqrt(sum(theta^2)),
               m = theta / sqrt(sum(theta^2))),
    SE = SE,
    info = out
  ))
}

#vMF_paramvec <- function(m, k){
#  return(k*m)
#}

vMF_fromparamvec <- function(paramvec){
  k <- sqrt(sum(paramvec^2))
  m <- paramvec / k
  return(list(m = m, k = k))
}
