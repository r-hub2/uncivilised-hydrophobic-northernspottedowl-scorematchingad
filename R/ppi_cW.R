#' @title Quickly Generate a Vector of Windham Exponents for the PPI Model
#' @description These functions help to quickly generate a set of Windham exponents for use in [`ppi_robust()`] or [`Windham()`].
#' Rows and columns of \eqn{A_L} and \eqn{b_L} corresponding to components with strong concentrations of probability mass near zero have non-zero constant tuning exponent, and all other elements have a tuning constant of zero.
#' All elements of \eqn{\beta} have a tuning exponent of zero.
#' @param cW The value of the non-zero Windham tuning exponents.
#' @param ... Values of `TRUE` or `FALSE` in the same order of the components specifying that a component has probability mass concentrated near zero.
#' @return A vector of the same length as the parameter vector of the PPI model. Elements of \eqn{A_L}{A_L} will have a value of `cW` if both their row and column component has probability mass concentrated near zero. Similarly, elements of \eqn{b_L} will have a value of `cW` if their row corresponds to a component that has a probability mass concentrated near zero. All other elements are zero.
#' @details
#' The Windham robustifying method involes weighting observations by a function of the proposed model density \insertCite{windham1995ro}{scorematchingad}.
#' \insertCite{scealy2024ro;textual}{scorematchingad} found that only some of the tuning constants should be non-zero:
#' the tuning exponents corresponding to \eqn{\beta} should be zero to avoid infinite weights;and to improve efficiency any rows or columns of \eqn{A_L} corresponding to components without concentrations of probability mass (i.e. outliers can't exist) should have exponents of zero.
#' \insertCite{scealy2024ro;textual}{scorematchingad} set the remaining tuning exponents to a constant.
#' @references \insertAllCited{}
#' @export
ppi_cW <- function(cW, ...){
  ellipsis::check_dots_used()
  bools <- unlist(list(...))
  bools <- as.logical(bools)
  if (any(is.na(bools))){stop("An argument after 'cW' is neither TRUE nor FALSE")}
  #if (tail(bools, 1)){stop("The final component must be set to FALSE")}

  p <- length(bools) 
  ALs_ww <- matrix(0, p-1, p-1)
  ALs_ww[bools[-p], bools[-p]] <- 1
  inWW <- ppi_paramvec(p, AL = ALs_ww, bL = bools[-p], beta = FALSE)

  cW <- cW * inWW
  return(cW)
}

#' @param Y A matrix of observations
#' @rdname ppi_cW
#' @description The function `ppi_cW_auto()` automatically detects concentrations near zero by fitting a PPI distribution with \eqn{A_L=0} and \eqn{b_L=0} (i.e. a Dirichlet distribution) with the centred log-ratio transformation of the manifold.
#' @examples
#' Y <- rppi_egmodel(100)$sample
#' ppi_cW_auto(0.01, Y)
#' ppi_cW(0.01, TRUE, TRUE, FALSE)
#' @export
ppi_cW_auto <- function(cW, Y){
  anest <- ppi(Y = Y, paramvec = ppi_paramvec(p = ncol(Y), AL=0, bL=0),
    trans = "clr")
  betaest <- anest$est$beta
  cW <- do.call(ppi_cW, c(list(cW = cW), as.list(betaest < 0)))
  return(cW)
}


