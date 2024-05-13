#' @title Robust Fitting of von Mises Fisher
#' @family directional model estimators
#' @description
#' Robust estimation for von Mises Fisher distribution using [`Windham()`].
#' @param cW Tuning constants for each parameter in the vMF parameter vector. If a single number then the constant is the same for each element of the parameter vector.
#' @param Y A matrix of observations in Cartesian coordinates.
#' @param ... Passed to [`Windham()`] and then passed onto [`vMF()`].
#' @family Windham robustness functions
#' @examples
#' if (requireNamespace("movMF")){
#'   Y <- movMF::rmovMF(1000, 100 * c(1, 1) / sqrt(2))
#' } else {
#'   Y <- matrix(rnorm(1000 * 2, sd = 0.01), ncol = 2)
#'   Y <- Y / sqrt(rowSums(Y^2))
#' }
#' vMF_robust(Y, cW = c(0.01, 0.01), method = "smfull")
#' vMF_robust(Y, cW = c(0.01, 0.01), method = "Mardia")
#' @export
vMF_robust <- function(Y, cW, ...){
  ellipsis::check_dots_used()
  extraargs <- list(...)

  # user friendly cW
  if (length(cW) == 1){
    if (!is.null(extraargs$paramvec)){
      isfixed = t_u2i(extraargs$paramvec)
      cW <- cW * !isfixed
    } else {
      cW <- rep(cW, ncol(Y))
    }
  }


  ldenfun <- function(Y, theta){ #here theta is km
    return(drop(Y %*% theta))
  }
  est <- Windham(Y = Y,
                     estimator = vMF,
                     ldenfun = ldenfun,
                     cW = cW,
                     ...)
  out <- list(
    est = c(list(paramvec = est$paramvec), vMF_fromparamvec(est$paramvec)),
    SE = "Not calculated.",
    info = est$optim
  )
  return(out)
}

vMF_kappa_robust <- function(Y, cW, ...){
  extraargs <- list(...)
  Y <- vMF_stdY(Y, w = extraargs$w)
  ellipsis::check_dots_used()
  ldenfun <- function(Y, theta){ #here theta is k and m is c(1, 0, ...)
    return(drop(theta * Y[, 1]))
  }
  est <- Windham(Y = Y,
                     estimator = vMF_kappa,
                     ldenfun = ldenfun,
                     cW = cW,
                     ...)
  out <- list(
    est = list(paramvec = est$paramvec, k = est$paramvec),
    SE = "Not calculated.",
    info = est$optim
  )
  return(out)
}
