#' @title Estimate the Fisher-Bingham Distribution
#' @family directional model estimators
#' @description Estimates parameters for the Fisher-Bingham distribution using score-matching.
#' @examples
#' p <- 3
#' A <- rsymmetricmatrix(p, -10, 10)
#' A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
#' m <- runif(p, -10, 10)
#' m <- m / sqrt(sum(m^2))
#' if (requireNamespace("simdd")){
#'   Y <- simdd::rFisherBingham(1000, 2 * m, A)
#'   FB(Y)
#' }
#' @inheritParams vMF
#' @param km Optional. A vector of the same length as the dimension, representing the parameter vector for the von Mises-Fisher component (i.e. the \eqn{\kappa \mu} see [`vMF()`]).
#' If supplied, the non-NA elements are fixed.
#' @param A Optional. The Bingham matrix. If supplied the non-NA elements of the Bingham matrix are fixed.
#' The final element of the diagonal of `A` must be NA as the software calculates this value to ensure the trace of the Bingham matrix is zero.
#' @details
#' The density of the Fisher-Bingham distribution is proportional to 
#' \deqn{\exp(z^TAz + \kappa\mu^Tz),}
#' where \eqn{A} is a matrix as in the Bingham distribution, and
#' \eqn{\kappa} and \eqn{\mu} are the concentration and mean direction, respectively, as in the von Mises-Fisher distribution.
#' 
#' # Warning: Slow Convergence with Sample Size
#' Score matching estimates of all elements of \eqn{A} and \eqn{\kappa\mu} converge slowly with sample size.
#' Even with a million simulated measurements,
#'  the gradient of the score matching discrepancy at the true parameters can have size (L2 Euclidean norm) more than 0.001, which is substantially non-zero.
#' @export
FB <- function(Y, km = NULL, A = NULL){
  p <- ncol(Y)
  if (is.null(A)){
    A <- matrix(NA, nrow = p, ncol = p)
  }
  stopifnot(all(dim(A) == c(p, p)))
  if (!is.na(A[p,p])){stop("The final diagonal element of matrix A cannot be fixed in this software. Please consider reordering your dimensions.")}
  if (is.null(km)){
    km <- rep(NA, p)
  }

  intheta <- FB_mats2theta(1, km, A)

  tapes <- buildsmdtape("sph","identity", "sph", "FB",
               rep(1, p)/sqrt(p), intheta,
               bdryw = "ones",
               verbose = FALSE)

  sminfo <- cppad_closed(tapes$smdtape, Y)
  theta <- intheta
  theta[is.na(intheta)] <- sminfo$est
  thetamat <- FB_theta2mats(theta)
  SE <- intheta * 0
  SE[is.na(intheta)] <- sminfo$SE
  SE <- FB_theta2mats(SE, isSE = TRUE)
  return(list(
   est = c(thetamat, list(paramvec = theta)),
   SE = c(SE, list(paramvec = sminfo$SE)),
   info = sminfo
  ))
}
#' @return
#' A list of `est`, `SE`, and `info`
#'  * `est` contains a slot for the estimate of \eqn{A} and \eqn{\kappa\mu}, and the vector form of the full Fisher-Bingham parameter set `paramvec`.
#'  * `SE` like `est`, but containing estimates of the standard error based on the sandwich method.
#'  * `info` a variety of information about the estimation process and results.

# non-normalised density function
#' @noRd
qdFB <- function(u, k, m, A){
  exp(u %*% A %*% u + k * m %*% u)
}

FB_mats2theta <- function(k, m, A){
  Binghamtheta <- Bingham_Amat2theta(A)
  theta <- c(Binghamtheta, k * m)
  return(theta)
}

FB_theta2mats <- function(theta, isSE = FALSE){
  #length of theta is:
  # l = p-1 + (p - 1) * p/2 + p
  # = p^2/2 + p * (1 - 1/2 + 1) - 1
  # = p^2/2 + p * (1.5) - 1
  # so
  # 0 = p^2 + 3*p - (2 + 2*l)
  # => p = (-3 pm sqrt(3^2 + 4(2+2l)))/2
  #  = (-3 pm sqrt(9 + 8 + 8l)))/2
  #  = (-3 pm sqrt(17 + 8l)))/2
  # = (-3 + sqrt(17 + 8l)))/2  # because sqrt(17 + 8l) >= 3
  p <- (-3 + sqrt(17 + 8 * length(theta)))/2
  A <- matrix(NA, ncol = p, nrow = p)
  diag(A) <- c(theta[1:(p-1)], NA)
  A[upper.tri(A)] <- theta[seq.int(p, length.out = (p - 1) * p/2)]
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  km <- theta[seq.int(p-1 + (p - 1) * p/2 + 1, length.out = p)]
  if (!isSE){
    A[p,p] <- -sum(diag(A)[1:(p-1)]) #trace is 0 constraint
    k <- sqrt(sum(km^2))
    m <- km/k
    return(list(
      k = k,
      m = m,
      km = km,
      A = A
    ))
  } else {
    A[p,p] <- NA #trace is 0 constraint
    return(list(
      km = km,
      A = A
    ))
  }
}
