#' @title Score Matching Estimators for the Bingham Distribution
#' @family directional model estimators
#' @inheritParams vMF
#' @param A For full score matching only: if supplied, then NA elements of `A` are estimated and the other elements are fixed. For identifiability the final element of `diag(A)` must be `NA`.
#' @description
#' Score matching estimators for the Bingham distribution's parameter matrix. Two methods are available: a full score matching method that estimates the parameter matrix directly and a hybrid method by \insertCite{mardia2016sc;textual}{scorematchingad} that uses score matching to estimate just the eigenvalues of the parameter matrix.
#' @details
#' The Bingham distribution has a density proportional to
#' \deqn{\exp(z^T A z),}
#' where \eqn{A} is a symmetric matrix and the trace (sum of the diagonals) of \eqn{A} is zero for identifiability \insertCite{@p181, @mardia2000di}{scorematchingad}.
#'
#' The full score matching method estimates all elements of \eqn{A} directly except the final element of the diagonal, which is calculated from the sum of the other diagonal elements to ensure that the trace of \eqn{A} is zero.
#'
#' The method by \insertCite{mardia2016sc;textual}{scorematchingad} first calculates the maximum-likelihood estimate of the eigenvectors \eqn{G} of \eqn{A}. 
#' The observations `Y` are then standardised to `Y`\eqn{G}. 
#' This standardisation corresponds to diagonalising \eqn{A}
#' where the eigenvalues of \eqn{A} become the diagonal elements of the new \eqn{A}.
#' The diagonal elements of the new \eqn{A} are then estimated using score matching, with the final diagonal element calculated from the sum of the other elements.
#' See \insertCite{mardia2016sc;textual}{scorematchingad} for details.
#' @references \insertAllCited{}
#' @examples
#' p <- 4
#' A <- rsymmetricmatrix(p)
#' A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
#' if (requireNamespace("simdd")){
#'   Y <- simdd::rBingham(100, A)
#'   Bingham(Y, method = "Mardia")
#' }
#' @return
#' A list of `est`, `SE` and `info`.
#'  * `est` contains the estimated matrix `A` and a vector form, `paramvec`, of `A` (ordered according to `c(diag(A)[1:(p-1)], A[upper.tri(A)])` ). For the Mardia method, the estimated eigenvalues of `A` (named `evals`) and eigenvectors of `A` (named `G`) are also returned.
#'  * `SE` contains estimates of the standard errors if computed. See [`cppad_closed()`].
#'  * `info` contains a variety of information about the model fitting procedure and results.
#' @export
Bingham <- function(Y, A = NULL, w = rep(1, nrow(Y)), method = "Mardia"){
  if (method == "smfull"){
    out <- Bingham_full(Y, A = A, w = w)}
  if (method %in% c("Mardia", "hybrid")){
    stopifnot(is.null(A))
    out <- Bingham_Mardia(Y, w = w)
    }
  return(out)
}

Bingham_full <- function(Y,  A = NULL, w = rep(1, nrow(Y))){
  p <- ncol(Y)
  if (is.null(A)){
    A <- matrix(NA, nrow = p, ncol = p)
  }
  stopifnot(all(dim(A) == c(p, p)))
  if (!is.na(A[p,p])){stop("The final diagonal element of matrix A cannot be fixed in this software. Please consider reordering your dimensions.")}

  intheta <- Bingham_Amat2theta(A)

  ytape <- rep(1, p) / sqrt(p)
  tapes <- buildsmdtape("sph","identity", "sph", "Bingham",
                           ytape, intheta,
                           bdryw = "ones")
  out <- cppad_closed(tapes$smdtape, Y = Y, w = w)
  theta <- intheta
  theta[is.na(intheta)] <- out$est
  A = Bingham_theta2Amat(theta)
  if (is.numeric(out$SE)){
    tSE <- intheta * 0
    tSE[is.na(intheta)] <- out$SE
    SE = Bingham_theta2Amat(tSE)
    SE[p, p] <- NA
  } else {
    tSE <- out$SE
    SE <- out$SE
  }

  return(list(
    est = list(A = A, paramvec = theta),
    SE = list(A = SE, paramvec = tSE),
    info = out
  ))
}

Bingham_Mardia <- function(Y, w = rep(1, nrow(Y))){
  Tmat <- t(Y) %*% (Y * w/sum(w))  #elements of Tmat are sample averages Yj*Yi being squares or cross products between the coordinates
  Tmat_es <- eigen(Tmat)
  Gammahat <- Tmat_es$vectors
  Ystd <- Y %*% Gammahat

  p <- ncol(Ystd)
  A <- matrix(0, nrow = p, ncol = p)
  diag(A) <- NA
  stopifnot(all(dim(A) == p))
  intheta <- Bingham_Amat2theta(A)
  tapes <- buildsmdtape("sph","identity", "sph", "Bingham",
                        rep(1, p) / sqrt(p), intheta,
                        bdryw = "ones")

  sm <- cppad_closed(tapes$smdtape, Y = Ystd, w = w)
  theta <- intheta
  theta[is.na(intheta)] <- sm$est
  Astd <- Bingham_theta2Amat(theta)
  Lambda <- diag(Astd)

  if (is.numeric(sm$SE)){
    SE <- intheta * 0
    SE[is.na(intheta)] <- sm$SE
    SE <- Bingham_theta2Amat(SE)
    SE[p, p] <- NA#final NA because final diagonal element is not estimated directly
    evalSE <- diag(SE)
  } else {
    evalSE <- sm$SE
  }
  A <- Gammahat %*% diag(Lambda) %*% t(Gammahat)
  return(list(
    est = list(A = A,
               evals = Lambda,
               G = Gammahat,
               paramvec = Bingham_Amat2theta(A)),
    SE = list(evals = evalSE),
    info = c(sm, list(Aforstddata = Astd))
  ))
}

Bingham_Amat2theta <- function(A){
  p <- ncol(A)
  stopifnot(isSymmetric(A))
  if(isTRUE(abs(sum(diag(A))) > 1E-8)){warning("Trace of A is not zero, final diagonal element of A will be altered.")}
  theta <- c(diag(A)[1:(p-1)], A[upper.tri(A)])
  return(theta)
}

Bingham_theta2Amat <- function(theta){
  #length of theta is p-1 + (p - 1) * p/2 = p2/2 - p/2 + p - 1 = p2/2 + p/2 - 1
  # => 0 = p2/2 + p/2 - 1 - length
  # => 0 = p2 + p - 2(1+length)
  # => p = (-1 +/- sqrt(1 + 4 * 2 * (1+length))) / 2
  p = (-1 + sqrt(1 + 4 * 2 * (1+length(theta)))) / 2  #the '-' alternative is always negative because sqrt(1 +epsilon) > 1
  A <- matrix(NA, ncol = p, nrow = p)
  diag(A) <- c(theta[1:(p-1)], NA)
  A[upper.tri(A)] <- theta[p:length(theta)]
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #trace is 0 constraint
  return(A)
}

