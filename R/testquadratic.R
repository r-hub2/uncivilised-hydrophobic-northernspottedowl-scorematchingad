#' @title Test Whether a CppAD Tape is a Quadratic Function
#' @family tape evaluators
#' @description
#' Uses the [`CppAD` parameter property](https://cppad.readthedocs.io/latest/fun_property.html#parameter) and derivatives (via [`tapeJacobian()`]) to test whether
#' the tape is quadratic.
#' @param tape An `ADFun` object.
#' @param xmat If non-`NULL` and `dynparamat` non-`NULL` then the third-order derivatives at independent variable values of the rows of `xmat` and dynamic parameters from the rows of `dynparammat` are tested.
#' @param dynparammat If non-`NULL` and `xmat` non-`NULL` then the third-order derivatives at independent variable values of the rows of `xmat` and dynamic parameters from the rows of `dynparammat` are tested.
#' @param verbose If TRUE information about the failed tests is printed.
#' @details Uses the `xtape` and `dyntape` values stored in `tape` to create new tapes.
#' A tape of the Hessian is obtained by applying [`tapeJacobian()`] twice, and then uses a [`CppAD` property](https://cppad.readthedocs.io/latest/fun_property.html#parameter) to test whether the Hessian is constant. A function of quadratic form should have constant Hessian.
#'
#' If `xmat` and `dynparammat` are non-`NULL` then `testquadratic()` also checks the Jacobian of the Hessian at `xmat` and `dynparammat` values. For quadratic form functions the Jacobian of the Hessian should be zero.
#' @return `TRUE` or `FALSE`
#' @examples
#' tapes <- buildsmdtape(
#'    "sim", "sqrt", "sph",
#'    ll = "ppi",
#'    ytape = c(0.2, 0.3, 0.5),
#'    usertheta = ppi_paramvec(p = 3), 
#'    bdryw = "minsq",
#'    acut = 0.1,
#'    verbose = FALSE)
#'
#'  testquadratic(tapes$smdtape)
#' @export
testquadratic <- function(tape, xmat = NULL, dynparammat = NULL, verbose = FALSE){
  stopifnot(inherits(tape, "ADFun"))
  tapeJ <- tapeJacobian(tape)
  tapeH <- tapeJacobian(tapeJ)

  #pParameter() test
  isparameter <- pParameter(tapeH$ptr)
  result_pParameter <- all(isparameter)
  if (verbose && !result_pParameter){
    message(sprintf("The Hessian was non-constant according to pParameter() for elements %s.",
                    paste(which(!isparameter), collapse = ", ")))
  }

  if (is.null(xmat) && is.null(dynparammat)){return(result_pParameter)}

  #try Jacobian
  stopifnot(isTRUE(nrow(xmat) == nrow(dynparammat)))
  thirdderivs <- lapply(1:nrow(xmat), function(i){
    pJacobian(tapeH$ptr, xmat[i, ], dynparammat[i, ])
  })
  isallzero <- unlist(lapply(thirdderivs, function(vec){all(vec == 0)}))
  result_thirdderiv <- all(isallzero)
  if (verbose && !result_thirdderiv){
    message(sprintf("The Jacobian of the Hessian was non-zero for row %s of xmat and dynparammat",
                    paste(which(!isallzero), collapse = ", ")))
  }

  finalresult <- result_thirdderiv && result_pParameter
  return(finalresult)
}


