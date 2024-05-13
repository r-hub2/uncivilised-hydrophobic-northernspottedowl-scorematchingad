#' @title Evaluate a CppAD Tape Many Times
#' @family tape evaluators
#' @param tape An [`ADFun`] object (i.e. a tape of a function).
#' @param xmat A matrix of (multivariate) independent variables where each represents a single independent variable vector. Or a single independent variable vector that is used for all rows of `pmat`.
#' @param pmat A matrix of dynamic parameters where each row specifies a new set of values for the dynamic parameters of `tape`. Or a single vector of dynamic parameters to use for all rows of `xmat`.
#' @param xcentres A matrix of approximation for Taylor approximation centres for `xmat`. Use values of `NA` for rows that do not require Taylor approximation.
#' @param approxorder Order of Taylor approximation
#' @description Evaluates a tape exactly or approximately for an array of provided variable values and dynamic parameter values.
#' The function `evaltape_wsum()` computes the column-wise weighted sum of the result.
#' @details
#' Approximation is via Taylor approximation of the independent variable around the approximation centre provided in `xcentres`.
#' @return
#' A matrix, each row corresponding to the evaluation of the same row in `xmat`, `pmat` and `xcentres`.
#' @examples
#' u <- rep(1/3, 3)
#' tapes <- buildsmdtape("sim", "sqrt", "sph", "ppi",
#'               ytape = u,
#'               usertheta = ppi_paramvec(p = 3),
#'               bdryw = "minsq", acut = 0.01,
#'               verbose = FALSE
#'               )
#' evaltape(tapes$lltape, u, rppi_egmodel(1)$theta)
#' evaltape(tapes$smdtape, rppi_egmodel(1)$theta, u)
#' evaltape(tapes$lltape, rbind(c(0, 0, 1), c(0,0,1)), 
#'          rppi_egmodel(1)$theta, 
#'          xcentres = rbind(c(0.0005, 0.0005, 0.999), NA))
#' @export
evaltape <- function(tape, xmat, pmat, xcentres = NA * xmat, approxorder = 10){
  if (R6::is.R6(tape) && inherits(tape, "ADFun")){
     tape <- tape$ptr
  }
  stopifnot(nrow(xmat) == nrow(xcentres))
  if (is.vector(xmat)){xmat <- matrix(xmat, ncol = length(xmat))}
  if (is.vector(pmat)){pmat <- matrix(pmat, ncol = length(pmat))}
  if (isTRUE(nrow(xmat) == 1)){
    xmat <- matrix(xmat, byrow = TRUE, nrow = nrow(pmat), ncol = length(xmat))
    xcentres <- matrix(xcentres, byrow = TRUE, nrow = nrow(pmat), ncol = length(xcentres))
  }
  if (isTRUE(nrow(pmat) == 1)){
    pmat <- matrix(pmat, byrow = TRUE, nrow = nrow(xmat), ncol = length(pmat))
  }
  stopifnot(nrow(xmat) == nrow(pmat))
  toapprox <- !is.na(xcentres[, 1])

  evals_l <- list()
  # exact evaluations
  if (any(!toapprox)){
    evals_l[!toapprox] <- lapply(which(!toapprox), function(i){
      pForward0(tape, xmat[i, ], pmat[i, ])
    })
  }
  if (any(toapprox)){
    evals_l[toapprox] <- lapply(which(toapprox), function(i){
      pTaylorApprox(tape, xmat[i, ], xcentres[i, ], pmat[i, ], approxorder)
    })
  }

  evals <- do.call(rbind, evals_l)
  return(evals)
}

#' @rdname evaltape
#' @param w Weights to apply to each row of `xmat` for computing the weighted sum. If `NULL` then each row is given a weight of `1`.
evaltape_wsum <- function(tape, xmat, pmat, w=NULL, xcentres = NA * xmat, approxorder = 10){
  evals <- evaltape(tape, xmat = xmat, pmat = pmat,
                     xcentres = xcentres, approxorder = approxorder)
  
  # do weight checks afterwards so that eval results can be used to choose weights
  out <- wcolSums(evals, w)
  return(out)
}

