#' @noRd
#' @title Windham weights for a given parameter vector
#' @description Evaluates \eqn{f(z, c\circ\theta)}, which is the density \eqn{f} at observation \eqn{z}, given a parameter set \eqn{\theta} and vector of tuning constants \eqn{c}. The multiplication \eqn{\circ} is element-wise.
#' These are the weights used by [`Windham()`].
#' The density function is passed to `Windham_weights()` in log form.
#' @param cW The robustness tuning constants. A vector of the same length as `theta`. Easily created for the PPI model using [`ppi_cW()`] and [`ppi_cW_auto()`].
#' @param ldenfun A (possibly improper) log density function taking two arguments, `Y` and `theta`.
#' @param theta Parameters for the model. The element-wise multiplication `cW * theta` is passed to `ldenfun`.
#' @param Y A matrix of measurements. Each row a measurement. Passed to `ldenfun`.
#' @return
#' A vector of weights corresponding to the rows of `Y`.
#' The weights are normalised to sum to 1.
Windham_weights <- function(ldenfun, Y, theta, cW){
  if (is.null(ldenfun)){stop("ldenfun is NULL")}
  stopifnot(length(cW) == length(theta))
  stopifnot(is.numeric(cW))
  thetaforweights <- cW * theta #the elements of theta with FALSE inWW will be set to zero
  weights <- exp(ldenfun(Y = Y, theta = thetaforweights))
  weights=weights/sum(weights)
  if (any(!is.finite(weights))){#when lden values too huge for numerics, set them to the max non-Inf value
    stop(paste0("Non-finite weights for observations: ", 
               paste(which(!is.finite(weights)), collapse = ", "), " at theta of ",
               paste(theta, collapse = " "),
               ". Try a different weight vector 'cW', a different 'paramvec_start' or a different fixed point search method."))
  }
  return(weights)
}

