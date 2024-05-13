#' @title Score Matching Estimator for Quadratic-Form Score Matching Discrepancies
#' @family generic score matching tools
#' @description 
#' For a tape (i.e. an `ADFun` object) of a quadratic-form score matching discrepancy function, calculates the vector of parameters such that the gradient of the score matching discrepancy is zero.
#' Also estimates standard errors and covariance.
#' Many score matching estimators have an objective function that has a quadratic form (see [`scorematchingtheory`]).
#' @param Yapproxcentres A matrix of Taylor approximation centres for rows of Y that require approximation. `NA` for rows that do not require approximation.
#' @param smdtape A `CppAD` tape of a score matching discrepancy function that has *quadratic form*. Test for quadratic form using [`testquadratic()`].
#' The `smdtape`'s independent variables are assumed to be the model parameters to fit
#' and the `smdtape`'s dynamic parameter is a (multivariate) measurement.
#' @param Y A matrix of multivariate observations. Each row is an observation.
#' @param w Weights for each observation.
#' @param approxorder The order of Taylor approximation to use.
#' @details
#' When the score matching discrepancy function is of quadratic form, then the gradient of the score matching discrepancy is zero at \eqn{H^{-1}b}{solve(H) %*% b},
#' where \eqn{H} is the average of the Hessian of the score matching discrepancy function evaluated at each measurement and 
#' \eqn{b} is the average of the gradient offset (see [`quadratictape_parts()`]) evaluated at each measurement. 
#' Both the Hessian and the gradient offset are constant with respect to the model parameters for quadratic-form score matching discrepancy functions.
#'
#' Standard errors use the Godambe information matrix (aka sandwich method) and are only computed when the weights are constant.
#' The estimate of the negative of the sensitivity matrix \eqn{-G} is
#' the average of the Hessian of `smdtape` evaluated at each observation in `Y`.
# \deqn{\hat{G(\theta)} = \hat{E} -H(smd(\theta;Y))),}
# where \eqn{smd} is the score matching discrepancy function represented by `smdtape`,
# \eqn{H} is the Hessian with respect to \eqn{\theta}, which is constant for quadratic-form functions,
# 
#' The estimate of the variability matrix \eqn{J} is 
#' the sample covariance (denominator of \eqn{n-1}) of the gradiant of `smdtape` evaluated at each of the observations in `Y` for the estimated \eqn{\theta}.
# \deqn{\hat{J}(\theta) = var(grad(w smd(\theta;Y))),}

#' The estimated variance of the estimator is then as
#' \eqn{G^{-1}JG^{-1}/n,}
# \deqn{\hat{G}(\theta)^{-1}\hat{J}(\theta)\hat{G}(\theta)^{-1}/n,}
#' where `n` is the number of observations.
#'
#' Taylor approximation is available because boundary weight functions and transformations of the measure in Hyvärinen divergence can remove singularities in the model log-likelihood, however evaluation at these singularities may still involve computing intermediate values that are unbounded.
#' If the singularity is ultimately removed, then Taylor approximation from a nearby location will give a very accurate evaluation at the removed singularity.
#' @examples
#' smdtape <- buildsmdtape("sim", "sqrt", "sph", "ppi",
#'               ytape = rep(1/3, 3),
#'               usertheta = ppi_paramvec(p = 3),
#'               bdryw = "minsq", acut = 0.01,
#'               verbose = FALSE
#'               )$smdtape
#' Y <- rppi_egmodel(100)
#' cppad_closed(smdtape, Y$sample)
#' @export
cppad_closed <- function(smdtape, Y, Yapproxcentres = NA * Y, 
                         w = rep(1, nrow(Y)),
                         approxorder = 10){
  stopifnot(nrow(Y) == length(w))
  stopifnot(nrow(Y) == nrow(Yapproxcentres))
  
  parts <- quadratictape_parts(smdtape, tmat = Y,
                               tcentres = Yapproxcentres,
                               approxorder = approxorder)

  # weight parts
  partsT <- lapply(parts, wcolSums, w = w)
  offset <- partsT$offset
  Hess <- partsT$Hessian
  Hess <- matrix(Hess, ncol = sqrt(ncol(parts$Hessian)))
  invHess <- tryCatch(solve(Hess),
               error = function(e) {
     if (grepl("system.*singular", e)){
          stop(paste("Hessian of the score-matching discrepancy function is not invertible.", e))
        } else {stop(e)}
     })
  root <- drop(-1 * invHess %*% offset)

  # compute SEs
  SE <- "Not calculated."
  covar <- "Not calculated."
  if (all(w[[1]] == w)){
    covar <- cppad_closed_estvar(Y, root, parts$offset, parts$Hess)
    SE <- sqrt(diag(covar))
  }
  return(list(
    est = root,
    Hessian = Hess,
    offset = offset,
    SE = SE,
    covar = covar
  ))
}

# function for getting standard errors for the above function
cppad_closed_estvar <- function(Y, theta, offsets, Hesss){
  sens <- - matrix(colMeans(Hesss), ncol = sqrt(ncol(Hesss)))
  sensinv <- solve(sens)

  # compute gradients
  grads <- lapply(1:nrow(offsets), function(i){
    drop(matrix(Hesss[i, ], ncol = ncol(offsets)) %*% theta + 
      t(offsets[i, , drop = FALSE]))
  })
  grads <- do.call(rbind, grads)
  variability <- stats::cov(grads)

  Ginfinv <- sensinv %*% variability %*% sensinv / nrow(offsets) #inverse of the Godambe information matrix, also called the sandwich information matrix
  return(Ginfinv) 
}


