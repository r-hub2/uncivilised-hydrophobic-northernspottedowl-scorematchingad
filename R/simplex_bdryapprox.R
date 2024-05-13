#' @noRd
#' @title Generate approximation centres for measurements on the boundary of the simplex
#' @description 
#' Move compositional measurements by `shiftsize` towards the centre of the simplex.
#' @details Measurements on the boundary of the simplex can lead to degenerate log-likelihood values, which are avoided in score matching by the boundary weight function. However, `CppAD` still computes the log-likelihood and derivatives, and (to KLH's knowledge) is not able to evaluate the value of \eqn{0/0} or similar.
#' Good approximations of the score matching discrepancy value can be obtained through Taylor approximations, which require a *centre*.
#' The function here generates this approximation centre.
#' @param Y The compositional measurements
#' @param shiftsize The distance to away from `u` to create the approximation centre.
#' @return A matrix the same size as `Y`, containing the measurements in `Y` translated `shiftsize` distance towards the simplex centre.
#' @examples
#' m <- rppi_egmodel(10)
#' simplex_boundaryshift(m$sample, shiftsize = 1E-5)
simplex_boundaryshift <- function(Y, shiftsize = 1E-4){
  p <- ncol(Y)
  middleofsimplex <- rep(1/p, p)
  shiftdir <- middleofsimplex - Y 
  shiftdir <- shiftdir / sqrt(rowSums(shiftdir^2)) #make a unit vector
  approxcentre <- Y + shiftsize * shiftdir
  approxcentre <- approxcentre / rowSums(approxcentre) # normalise to get back onto simplex
  return(approxcentre)
}

#' @noRd
#' @title Determine whether points are on the simplex boundary
#' @description Tests whether points are within `bdrythreshold` distance of the boundary, where distance is the size of minimum component. 
#' @param Y Measurement matrix.
#' @param bdrythreshold Measurements closer than `bdrythreshold` to the edge of the simplex will be considered boundary measurements.
#' @return A vector of `TRUE` or `FALSE` values. Rows of `Y` corresponding to `TRUE` are on the boundary of the simplex.
simplex_isboundary <- function(Y, bdrythreshold = 1E-15){
  isbdrypt <- apply(Y, MARGIN = 1, min) < bdrythreshold
  return(isbdrypt)
}


# OBSOLETE - JUST USED FOR TESTING
#' @noRd
#' @title Detect boundary observations and generate Taylor approximation centres
#' @description 
#' Splits the matrix of observations `Y` into boundary points and interior points using [`simplex_isboundary()`].
#' For the boundary points Taylor approximation centres are computed with [`simplex_boundaryshift()`].
#' A weight vector, if supplied, is split consistently with the above split.
#' @details Measurements on the boundary of the simplex can lead to degenerate log-likelihood values, which are avoided in score matching by the boundary weight function. However, `CppAD` still computes the log-likelihood and derivatives, and (to KLH's knowledge) is not able to evaluate the value of \eqn{0/0} or similar.
#' Very good approximations of the score matching discrepancy value can be obtained through Taylor approximations, which requires a centre.
#' The function converts the data matrix `Y` into the necessary sub-data sets.
#' @param Y A matrix of compositional measurements, each row is a multivariate measurement.
#' @param bdrythreshold See [`simplex_isboundary()`].
#' @param shiftsize The distance to away from a measurement to create the measurements approximation centre. Passed to [`simplex_boundaryshift()`]
#' @return A list of
#'  * `interior` Matrix of interior measurements.
#'  * `uboundary` Matrix of measurements on the boundary.
#'  * `boundaryapprox` Corresponding approximation centres for `uboundary`.
#'  * `winterior` Weights from `w` corresponding to `interior` (`NULL` if not weights supplied).
#'  * `wboundary` Weights from `w` corresponding to `uboundary` (`NULL` if not weights supplied).
#' @examples
#' m <- rppi_egmodel(10)
#' simplex_boundaryshift(m$sample, shiftsize = 1E-5)
simplex_boundarysplit <- function(Y, bdrythreshold = 1E-15, shiftsize = 1E-10, w = NULL){
  onbdry <- simplex_isboundary(Y, bdrythreshold = bdrythreshold)
  acentres <- simplex_boundaryshift(Y = Y[onbdry, , drop = FALSE], shiftsize = shiftsize)
  return(list(
    interior = Y[!onbdry, , drop = FALSE],
    uboundary = Y[onbdry, , drop = FALSE],
    boundaryapprox = acentres,
    winterior = w[!onbdry, drop = FALSE],  #when w NULL, then this is NULL too
    wboundary = w[onbdry, drop = FALSE]))
}

