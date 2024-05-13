#' @rdname ppi_param_tools
#' @order 4
#' @details
#' The `Astar` parametrisation rewrites the PPI density as proportional to 
#' \deqn{\exp(z^TA^*z)\prod_{i=1}^p z_i^{\beta_i},}
#' where \eqn{A^*} (`Astar`) is a \eqn{p} by \eqn{p} matrix.
#' Because \eqn{z} lies in the simplex (in particular \eqn{\sum z_i = 1}), the density is the same regardless of the value of \eqn{1^T A^* 1}=`sum(Astar)`, where \eqn{1} is the vector of ones. Thus \eqn{A_L} and \eqn{b_L} specify \eqn{A^*} up to an additive factor. In the conversion `ppi_toAstar()`, \eqn{A^*} is returned such that \eqn{1^T A^* 1 = 0}.
#' `NULL` values or `NA` elements are not allowed for `ppi_toAstar()` and `ppi_fromAstar()`.
#' @examples
#' Astar <- rWishart(1, 6, diag(3))[,,1]
#' ppi_fromAstar(Astar)
#' ppi_toAstar(ppi_fromAstar(Astar)$AL, ppi_fromAstar(Astar)$bL)
#' @return `ppi_toAstar()`: The matrix \eqn{A^*}.
#' @export
ppi_toAstar <- function(AL, bL){
  # assumes DC = 0 initially
  p <- nrow(AL) + 1
  onesL <- rep(1, p-1)
  DB <- bL / 2
  DL <- AL + DB %*% t(onesL) + onesL %*% t(DB)
  Astar <- matrix(NA, p, p)
  Astar[1:(p-1), 1:(p-1)] <- DL
  Astar[1:(p-1), p] <- DB
  Astar[p,1:(p-1)] <- t(DB)
  Astar[p,p] <- 0
  # set everything to sum to 0 i.e 1Astar1 = 0
  Astar <- Astar - sum(Astar)/(p*p)
  return(Astar)
}

#' @rdname ppi_param_tools
#' @order 5
#' @param Astar The \eqn{A^*} matrix (a p by p symmetric matrix)
#' @return `ppi_fromAstar()`: A list of the matrix \eqn{A_L}, the vector \eqn{b_L} and a discarded constant.
#' @export
ppi_fromAstar <- function(Astar){
  stopifnot(isSymmetric(Astar))
  p = ncol(Astar)
  DL = Astar[1:(p-1), 1:(p-1)]
  DB = Astar[1:(p-1), p]
  DC = Astar[p, p]

  onesL <- rep(1, p-1)

  AL <- DL - DB %*% t(onesL) - onesL %*% t(DB) + DC * onesL %*% t(onesL)
  bL <- 2 * (DB - DC * onesL)
  return(list(
    AL = AL,
    bL = bL,
    const = DC))
}
