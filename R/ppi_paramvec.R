#' @title PPI Parameter Tools
#' @order 1
#' @family PPI model tools
#' @description The defualt parameterisation of the PPI model is a symmetric covariance-like matrix \eqn{A_L}, a location-like vector \eqn{b_L} and a set of Dirichlet exponents \eqn{\beta}. For `p` components, \eqn{A_L} has `p-1` rows, \eqn{b_L} is a vector with `p-1` elements and \eqn{\beta} is a vector with `p` elements.
#' For score matching estimation this form of the parameters must be converted into a single parameter vector using `ppi_paramvec()`.
#' `ppi_paramvec()` also includes easy methods to set parameters to `NA` for estimation with [`ppi()`] (in [`ppi()`] the NA-valued elements are estimated and all other elements are fixed).
#' The reverse of `ppi_paramvec()` is `ppi_parammats()`.
#' An alternative parametrisation of the PPI model uses a single `p` by `p`  matrix \eqn{A^*} instead of \eqn{A_L} and \eqn{b_L}, and for identifiability \eqn{A^*} is such that \eqn{1^T A^* 1 = 0} where \eqn{1=(1,1,...,1)} and \eqn{0=(0,0,..., 0)} \insertCite{scealy2023sc}{scorematchingad}.
#' Convert between parametrisations using `ppi_toAstar()` and `ppi_fromAstar()`.
#' @inheritSection ppi PPI Model
#' @rdname ppi_param_tools
#' @name ppi_param_tools
NULL

#' @rdname ppi_param_tools
#' @order 2
#' @return `ppi_paramvec()`: a vector of length \eqn{p + (p-1)(2 + (p-1)/2)}.
#' @details
#' `ppi_paramvec()` returns a vector starting with the diagonal elements of \eqn{A_L}, then the off-diagonal elements extracted by [`upper.tri()`] (which extracts elements of \eqn{A_L} along each row, left to right, then top to bottom), then \eqn{b_L}, then \eqn{\beta}.
#' @param AL Either `NULL`, a p-1 x p-1 symmetric matrix, a number, or "diag".
#' If NULL then all \eqn{A_L} elements will be set to NA.
#' If a single number, then \eqn{A_L} will be fixed as a matrix of the given value.
#' If "diag" then the non-diagonal elements of \eqn{A_L} will be fixed to 0, and the diagonal will be `NA`.
#' @param bL Either `NULL`, a number, or a vector of length p-1.
#' If `NULL`, then all elements of \eqn{b_L} will be set to `NA`.
#' If a single number, then \eqn{b_L} will be fixed at the supplied value.
#' @param Astar  Either `NULL` or a p x p matrix.
#' If non-null, then overrides `AL` and `bL`.
#' If a matrix, all elements must be non-`NA` and `Astar` will be converted to `AL` and `bL` using [`ppi_fromAstar()`].
#' @param beta Either `NULL`, a number, or a vector of length p.
#' If NULL then all elements of \eqn{\beta} will be set to `NA`.
#' If a single number then the \eqn{\beta} elements will be fixed at the given number.
#' @param betaL Either `NULL`, a number, or a vector of length p-1.
#' If `NULL` then the 1...(p-1)th \eqn{\beta} elements will be set to `NA`.
#' If a single number then the 1...(p-1)th \eqn{\beta} elements will be fixed at the given number.
#' @param betap Either `NULL` or a number.
#' If `NULL` then the `p`th element of \eqn{\beta} will be set to `NA`, and [`ppi()`] will estimate it.
#' If a number, then the pth element of \eqn{\beta} will be fixed at the given value.
#' @param p The number of components. If `NULL` then `p` will be inferred from other inputs.
#' @examples
#' ppi_paramvec(AL = "diag", bL = 0, betap = -0.5, p = 3)
#' @export
ppi_paramvec <- function(p = NULL, AL = NULL, bL = NULL, Astar = NULL, beta = NULL, betaL = NULL, betap = NULL){
  # search for a p
  if (is.null(p)){
    if (!is.null(Astar)){p <- nrow(Astar)}
    else if (isTRUE(is.matrix(AL))){p <- nrow(as.matrix(AL)) + 1}
    else if (isTRUE((is.matrix(bL) || is.vector(bL)) && (length(bL) > 1))){p <- length(as.vector(bL)) + 1}
    else if (isTRUE( (is.matrix(betaL) || is.vector(betaL)) && (length(betaL) > 1))){p <- length(as.vector(betaL)) + 1}
    else if (isTRUE( (is.matrix(beta) || is.vector(beta)) && (length(beta) > 1))){p <- length(as.vector(beta))}
    else {stop("Could not guess 'p' from other arguments. Please specify 'p'.")}
  }
  stopifnot(is.vector(p, mode = "numeric"))
  stopifnot(length(p) == 1)

  # initialise parameter objects
  bLprep = rep(NA, p-1)
  betaLprep = rep(NA, p-1)
  betapprep = NA
  # AL, bL and A
  if (!is.null(Astar)){
    if (!(is.null(AL) & is.null(bL))){warning("AL, bL and Astar supplied. Astar argument will override AL and bL.")}
    translated <- ppi_fromAstar(Astar)
    ALprep <- translated$AL
    bLprep <- translated$bL
  } else {
    # AL
    if (is.null(AL)){
      ALprep = matrix(NA, nrow = p-1, ncol = p-1) #could also do nothing
    } else if (is.matrix(AL)){
    # If a matrix, then the NA elements will be estimated and the others will be fixed at the supplied value (i.e. not estimated).
      if(!isSymmetric.matrix(AL)){stop("AL must be symmetric.")}
      ALprep = AL
    } else if (is.numeric(AL)){#' If a single number, then AL will be fixed as a matrix of the given value.
      ALprep = matrix(AL, nrow = p-1, ncol = p-1)
    } else if (is.character(AL)){#' If "diag" then the non-diagonal elements of AL will be fixed to 0.
      stopifnot((AL == "diag") | (AL == "d") | (AL == "diagonal"))
      ALprep <- matrix(0, nrow = p-1, ncol = p-1)
      diag(ALprep) <- NA
    } else if (is.logical(AL)){
      ALprep = matrix(AL, nrow = p-1, ncol = p-1) #covers NA, TRUE, and FALSE
    } else {
      stop("AL is not of required type.")
    }
    #bL
    # If a number, then bL will be fixed at the supplied value.
    if (!is.null(bL)){
      if (is.matrix(bL)){
         if ((nrow(bL) == 1) || (ncol(bL) == 1)){bL <- as.vector(bL)}
         else {stop("bL is a matrix with multiple rows and columns")}
      }
      if (!is.vector(bL, mode = "any")){stop("bL must be a vector or value")}
      if (length(bL) == 1){
        bLprep = rep(bL, p-1)
      } else {#' If a vector, then the NA elements will be estimated and the others will be fixed at the supplied value.
        if (length(bL) != p-1){stop("bL must have length p-1")}
        bLprep = bL
      }
    }
  }

  # beta
  if (!is.null(betaL)){
    stopifnot(is.vector(betaL, "numeric") | is.vector(betaL, "logical"))
    stopifnot(is.null(beta))
    if (length(betaL) == 1){
      # If a number then the 1...(p-1) beta elements fixed at the given number.
      betaLprep = rep(betaL, p-1)
    } else if (length(betaL) == p-1){
      # If a vector, then the NA elements will be estimated and the others will be fixed at the supplied value.
      betaLprep = betaL
    } else {
      stop("betaL must have length p-1")
    }
  }
  if (!is.null(betap)){
    stopifnot(length(betap) == 1)
    stopifnot(is.numeric(betap) | is.logical(betap))
    stopifnot(is.null(beta))
    betapprep = betap
  }
  if (!is.null(beta)){
    stopifnot(is.null(betaL))
    stopifnot(is.null(betap))
    if (is.matrix(beta)){beta <- drop(beta)}
    stopifnot(is.vector(beta, "numeric") | is.vector(beta, "logical"))
    if (length(beta) == 1){
      betaLprep = rep(beta, p-1)
      betapprep = beta
    } else if (length(beta) == p){
      # If a vector, then the NA elements will be estimated and the others will be fixed at the supplied value.
      betaLprep = beta[1:(p-1)]
      betapprep = beta[p]
    } else {
      stop("beta must have length p")
    }
  }

  # combine above preparation into a vector, NA values to be estimated
  beta = c(betaLprep, betapprep)
  names(beta) <- paste0("beta", 1:length(beta))
  names(bLprep) <- paste0("bL", 1:length(bLprep))
  stopifnot(isSymmetric.matrix(ALprep))
  ALprep_vec <- fromsmatrix(ALprep)
  names(ALprep_vec) <- paste0("AL", names(ALprep_vec))
  theta <- c(ALprep_vec, bLprep, beta)
  return(theta)
}
