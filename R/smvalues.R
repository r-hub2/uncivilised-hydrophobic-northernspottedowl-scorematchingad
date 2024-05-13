#' @title Compute Score Matching Discrepancy Value, Gradient, and Hessian
#' @family tape evaluators
#' @description Computes a range of relevant information for investigating score matching estimators.
#' @inheritParams evaltape
#' @param smdtape A taped score matching discrepancy. Most easily created by [`buildsmdtape()`].
#' @details Computes the score matching discrepancy function from [`scorematchingtheory`] or weighted sum of the score matching discrepancy function.
#' The gradient and Hessian are returned as arrays of row-vectors with each row corresponding to a row in `xmat` and `pmat`. 
#' Convert a Hessian row-vector to a matrix using `matrix(ncol = length(smdtape$xtape))`.
#' @examples
#' m <- rppi_egmodel(100)
#' smdtape <- buildsmdtape("sim", "sqrt", "sph", "ppi",
#'               ytape = rep(1/m$p, m$p),
#'               usertheta = ppi_paramvec(beta = m$beta),
#'               bdryw = "minsq", acut = 0.01)$smdtape
#' smvalues(smdtape, xmat = m$sample, pmat = m$theta[1:5])
#' smvalues_wsum(smdtape, m$sample, m$theta[1:5])$grad/nrow(m$sample)
#' @export
smvalues <- function(smdtape, xmat, pmat, xcentres = NA * xmat, approxorder = 10){
  stopifnot(inherits(smdtape, "ADFun"))
  # prepare tapes
  Jsmdfun <- tapeJacobian(smdtape)
  Hsmdfun <- tapeJacobian(Jsmdfun)
  
  smdfun_u <- tapeSwap(smdtape) #don't use a boundary for taping!
  Jsmdfun_u <- tapeSwap(Jsmdfun)
  Hsmdfun_u <- tapeSwap(Hsmdfun)

  smdvals <- evaltape(smdfun_u, xmat, pmat, xcentres = xcentres, approxorder = approxorder)

  gradvals <- evaltape(Jsmdfun_u, xmat, pmat, xcentres = xcentres, approxorder = approxorder)
  hessvals <- evaltape(Hsmdfun_u, xmat, pmat, xcentres = xcentres, approxorder = approxorder)

  return(list(
    obj = smdvals,
    grad = gradvals,
    hess = hessvals
  ))
}
#' @return
#' A list of 
#'  + `obj` the score matching discrepancy values
#'  + `grad` the gradient of the score matching discrepancy
#'  + `hess` the Hessian of the score matching discrepancy

#' @rdname smvalues
#' @param w Weights to apply to each row of `xmat` for computing the weighted sum. If `NULL` then each row is given a weight of `1`.
#' @export
smvalues_wsum <- function(tape, xmat, pmat, w=NULL, xcentres = NA * xmat, approxorder = 10){
  evals_l <- smvalues(tape, xmat = xmat, pmat = pmat,
                     xcentres = xcentres, approxorder = approxorder)
  
  # do weight checks afterwards so that eval results can be used to choose weights
  out <- lapply(evals_l, wcolSums, w = w)
  return(out)
}

