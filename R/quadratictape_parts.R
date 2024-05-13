#' @title Evaluate the Hessian and Gradient Offset of a Taped Quadratic Function
#' @family tape evaluators
#' @description
#' When the score matching discrepancy function is quadratic then the gradient of the score matching discrepancy function can be written using the Hessian and an offset term. This can be useful for solving for the situation when the gradient is zero.
#' The Hessian and offset term are computed using `CppAD` tapes.
#' Taylor approximation can be used for locations at removed singularities (i.e. where intermediate values are unbounded).
#' `quadratictape_parts()` will error if `testquadratic(tape)` returns `FALSE`.
#' @details
#' A quadratic function can be written
#' \deqn{f(x; t) = \frac{1}{2} x^T W(t) x + b(t)^T x + c,}
#' where \eqn{t} is considered a vector that is constant with respect to the differentiation.
#' The Hessian of the function is with respect to \eqn{x} is
#' \deqn{H f(x; t) = \frac{1}{2}(W(t) + W(t)^T).}
#' The gradient of the function with respect to \eqn{x} can then be written
#' \deqn{\Delta f(x;t) = H f(x; t) x + b(t)^T x,}
#' where the Hessian and offset \eqn{b(t)} depend only on \eqn{t}.
#'
#' The functions here evaluate the Hessian and offset \eqn{b(t)} for many values of \eqn{t}.
#' Tapes of the Hessian and gradient offset are created using [`tapeHessian()`] and [`tapeGradOffset()`] respectively.
#' These tapes are then evaluated for every row of `tmat`.
#' When the corresponding `tcentres` row is not `NA`, then approximate (but very accurate) results are calculated using Taylor approximation around the location given by the row of `tcentres`.
#' 
#' For score matching \eqn{x} is the set of model parameters and the vector \eqn{t} is a (multivariate) measurement.
#' @param tape A tape of a quadratic function where the independent and dynamic parameters correspond to the \eqn{x} and \eqn{t} in the details section, respectively. For score matching `tape` should be a tape of the score matching discrepancy function \eqn{A(z) + B(z) + C(z)} in [`scorematchingtheory`] with \eqn{z} the *dynamic parameters* and the model parameters the *independent variable* (which is the usual for the return of [`buildsmdtape()`]).
#' @param tmat A matrix of vectors corresponding to values of \eqn{t} (see details). Each row corresponds to a vector. For score matching, these vectors are measurements.
#' @param tcentres A matrix of Taylor approximation centres for rows of `tmat` that require approximation. `NA` for rows that do not require approximation.
#' @param approxorder The order of the Taylor approximation to use.
#' @return A list of
#'  + `offset` Array of offsets \eqn{b(t)}, each row corresponding to a row in `tmat`
#'  + `Hessian` Array of vectorised \eqn{H f(x; t)} (see [`tapeHessian()`]), each row corresponding to a row in `tmat`. For each row, obtain the Hessian in matrix format by using `matrix(ncol = length(tape$xtape))`.
#' @examples
#' u <- rep(1/3, 3)
#' smdtape <- buildsmdtape("sim", "sqrt", "sph", "ppi",
#'               ytape = u,
#'               usertheta = ppi_paramvec(p = 3),
#'               bdryw = "minsq", acut = 0.01,
#'               verbose = FALSE
#'               )$smdtape
#' quadratictape_parts(smdtape, 
#'   tmat = rbind(u, c(1/4, 1/4, 1/2)))
#' @export
quadratictape_parts <- function(tape, tmat, tcentres = NA * tmat, approxorder = 10){
  stopifnot(inherits(tape, "ADFun"))
  stopifnot(nrow(tmat) == nrow(tcentres))
  stopifnot(testquadratic(tape))
  toapprox <- !is.na(tcentres[, 1])

  Hesstape <- tapeHessian(tape)
  OffsetTape <- tapeGradOffset(tape)
  Hesstape_switched <- tapeSwap(Hesstape) #Hesstape is wrt to x (which in smd world is actually the model parameter set), but we want it to be wrt to the dynamic parameter like OffsetTape is

  #exact and approximat evalution of Hess
  fakeparametermat <- matrix(0, 
                            ncol = length(tape$xtape),
                            nrow = nrow(tmat))   #fake because the results don't depend on it
  emptyparametermat <- matrix(vector(mode = "double"), ncol = 0, nrow = nrow(tmat)) 
  # for OffsetTape which has no dynamic parameters

  Hesss <- evaltape(Hesstape_switched, xmat = tmat, pmat = fakeparametermat, xcentres = tcentres, approxorder = approxorder)
  offsets <- evaltape(OffsetTape, xmat = tmat, pmat = emptyparametermat, xcentres = tcentres, approxorder = approxorder)

  return(list(
    offset = offsets,
    Hessian = Hesss
  ))
}

