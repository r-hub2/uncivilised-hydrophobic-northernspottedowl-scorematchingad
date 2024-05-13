#' @noRd
#' @title Mean direction and rotation for von Mises Fisher
#' @family directional model estimators
#' @description Computes the mean direction of a sample, which is also the maximum likelihood estimate for the mean direction parameter `m` in the von Mises Fisher distribution.
#' @param Y A matrix of observations in Cartesian coordinates. Each row is a single observation.
#' @param w Weights associated with each observation in `Y`.
vMF_m <- function(Y, w = NULL){
  if (is.null(w)){m <- colMeans(Y)}
  else {m <- apply(Y, MARGIN = 2, weighted.mean, w)}
  m <- m/sqrt(sum(m^2))
  return(m)
}

#' @noRd
#' @details `vMF_m()` Rotates a set of observations on a sphere to have mean direction of \eqn{(1, 0, 0, ..., 0)}.
#' @param m Mean direction. If omitted then it will be computed using `vMF_m()`.
#' @examples
#' if (requireNamespace("simdd"){
#'   Y <- simdd::rFisherBingham(nsim = 10, mu = c(1,2,3))
#' }
vMF_stdY <- function(Y, m = NULL, w = NULL){
  if(is.null(m)){
    m <- vMF_m(Y, w = w)
  }
  Rmat <- rotationmatrix(c(1, rep(0, ncol(Y)-1)), m)
  out <- Y %*% t(Rmat)
  return(out)
}

