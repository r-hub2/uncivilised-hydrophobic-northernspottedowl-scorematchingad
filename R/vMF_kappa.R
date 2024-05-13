#' @noRd
#' @title Estimate the concentration for a von Mises Fisher distribution
#' @family directional model estimators
#' @description Using score matching, estimates the concentration \eqn{\kappa} from a sample with mean direction of `c(1, 0, 0, ..., 0)`. 
#' Results by \insertCite{mardia2016sc;textual}{scorematchingad} and some experiments of our own suggest that good implementations of the maximum likelihood estimator (e.g. [`movMF::movMF()`] ) will out perform `vMF_kappa()`.
#'
#' Often a sample with mean direction of `c(1, 0, 0, ...., 0)` is created by estimating the mean direction and rotating the data such that the mean direction equals `c(1, 0, 0, ...)`.
#' Performing this mean direction estimate, then rotation, then estimating concentration with score matching correponds to the hybrid estimator by \insertCite{mardia2016sc;textual}{scorematchingad}.
#' @inherit vMF sections
#' @references \insertAllCited()
#' @details
#' The function `vMF_kappa()` estimates *only* the concentration \eqn{\kappa}, and assumes that \eqn{\mu} is \eqn{(1, 0, 0, ..., 0)}.
#' @return
#' A list of `est`, `SE` and `info`.
#' `est` contains the estimate of the concentration in the slot `k` for easy use and `paramvec` for compatibility with other functions in this package.
#' `SE` contains estimates of the standard errors if computed by the estimating method. See [`cppad_closed()`].
#' `info` contains a variety of information about the model fitting procedure.
#' @param Y A data matrix, each row is an observation.
#' @param w Weights corresponding to each row of `Y`.
vMF_kappa <- function(Y, w = rep(1, nrow(Y))){
  # do estimate, where all but the first component of theta are fixed at zero
  # because kappa * e1 = (kappa, 0, 0, 0, ...)
    p <- ncol(Y)
    tapes <- buildsmdtape("sph","identity", "sph", "vMF",
                          ytape = rep(1, p)/sqrt(p),
                          usertheta = c(NA, rep(0, p-1)))
    sminfo <- cppad_closed(tapes$smdtape, Y, w = w)
    k <- sminfo$est
    SE <- sminfo$SE
  return(list(
    k = k,
    SE = SE,
    sminfo = sminfo
  ))
}

