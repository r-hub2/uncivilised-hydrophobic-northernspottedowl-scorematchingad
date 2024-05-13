# functions for building custom likelihoods: build, evaluate (to check it). (Then use tapell and check it more!)
#' @title Compile a custom log-likelihood function.
# the following is because RcppEigen is needed in imports for cppFunction() below, but package checking thinks RcppEigen isn't imported
#' @import RcppEigen
#' @param includes passed to [`Rcpp::cppFunction()`]. For internal use.
#' @param cacheDir passed to [`Rcpp::cppFunction()`].
#' @param rebuild passed to [`Rcpp::cppFunction()`].
#' @param showOutput passed to [`Rcpp::cppFunction()`].
#' @param verbose passed to [`Rcpp::cppFunction()`].
#' @param code `C++` code for a log-likehood function (with normalising constant omitted if desired). See details for more.
#' @description
#' Supply `C++` code to specify a custom log-likelihood, much like `TMB::compile()` is passed `C++` code that formulate models.
#' For score matching the normalising constant of the log-likelihood can be omitted.
#' Taping of the function currently only works with the `gcc` compiler, typically on Linux; taping works if `customll_test()` returns `TRUE`.
#' @details
#' The function uses [`RcppXPtrUtils::cppXPtr()`] and [`Rcpp::cppFunction()`]. 
#' It is good practice to check the returned object using [`evalll()`].
#' The first compilation in a session can be very slow.
#' # Code Argument
#' `code` must be `C++` that uses only `CppAD` and `Eigen`, which makes it very similar to the requirements of the input to `TMB::compile()` (which also uses `CppAD` and `Eigen`).
#' 
#' The start of `code` should always be "`a1type fname(const veca1 &x, const veca1 &theta){`" where `fname` is your chosen name of the log-likelihood function, `x` represents a point in the data space and `theta` is a vector of parameters for the log-likelihood. This specifies that the function will have two vector arguments (of type `veca1`) and will return a single numeric value (`a1type`).
#' 
#' The type `a1type` is a double with special ability for being taped by `CppAD`. The `veca1` type is a vector of `a1type` elements, with the vector wrapping supplied by the `Eigen` C++ package (that is an `Eigen` matrix with 1 column and dynamic number of rows).
#' 
#' The body of the function must use operations from Eigen and/or CppAD, prefixed by `Eigen::` and `CppAD::` respectively. 
#' There are no easy instructions for writing these as it is genuine `C++` code, which can be very opaque to those unfamiliar with `C++`.
#' See the [Eigen documentation](https://eigen.tuxfamily.org/dox/group__QuickRefPage.html) for quick reference to available operations from Eigen. Limited operations are available directly from `CppAD` without `Eigen`: [unary operations](https://cppad.readthedocs.io/latest/unary_standard_math.html) and [binary operations](https://cppad.readthedocs.io/latest/binary_math.html). 
#' For the purposes of score matching the operations should all be smooth to create a smooth log-likelihood and the normalising constant may be omitted.
#' @examples
#' customll_test()
#' 
#' myll <- customll("a1type dirichlet(const veca1 &u, const veca1 &beta) {
#'   size_t d  = u.size();
#'   a1type y(0.);  // initialize summation at 0
#'   for(size_t i = 0; i < d; i++)
#'   {   y   += beta[i] * log(u[i]);
#'   }
#'   return y;
#' }")
#' evalll(myll, rep(1/3, 3), rep(-0.5, 3))
#'
#' tapes <- buildsmdtape("sim", "identity", "sim", 
#'  myll, rep(1/3, 3), rep(NA, 3), 
#'  bdryw="minsq", acut = 0.01)
#' evaltape(tapes$lltape, rep(1/3, 3), rep(-0.5, 3))
#' 
#' @returns `customll()` returns an `adloglikelood` object (which is just an `externalptr` with attributes) for the compiled log-likelihood function. The returned object has an attribute `fname`.
#' @export
customll <- function(code, rebuild = FALSE, 
                     includes = character(),
                     cacheDir = getOption("rcpp.cache.dir", tempdir()), 
                     showOutput = verbose, verbose = getOption("verbose")){
  ptr <- RcppXPtrUtils::cppXPtr(code, depends = c("RcppEigen", "scorematchingad"), includes = includes, rebuild = rebuild, cacheDir = cacheDir, showOutput = showOutput, verbose = verbose)
  
  # in the spirit of RcppXPtrUtils::checkXPtr check output and arguments. Rewritten here for bespoke error messages.
  if (!all(grepl("^const veca1", attr(ptr, "args", exact = TRUE)))){
    stop("Arguments incorrect. There should be two arguments, both of type 'const veca1'.")
  }
  if (attr(ptr, "type", exact = TRUE) != "a1type"){
    stop(sprintf("The return type is %s but should be a1type", attr(ptr, "type", exact = TRUE)))
  }
  ptr <- new_adloglikelihood(ptr)
  return(ptr)
}

# adloglikelood class things
new_adloglikelihood <- function(ptr, fname = NULL){
  if (is.null(attr(ptr, "fname", exact = TRUE)) && !is.null(fname)){
    attr(ptr, "fname") <- fname
  } else if (!is.null(attr(ptr, "fname", exact = TRUE)) && !is.null(fname)){
    warning("overwriting objects fname attribute")
    attr(ptr, "fname") <- fname
  } else if (is.null(attr(ptr, "fname", exact = TRUE)) && is.null(fname)){
    stop("object needs an fname")
  }
  class(ptr) <- c("adloglikelihood", class(ptr))
  return(ptr)
}
validate_adloglikelihood <- function(ll){
  ll_a <- unclass(ll)
  if (typeof(ll_a) != "externalptr"){
    stop("adloglikelihood is not an 'externalptr'")
  }
  if (is.null(attr(ll_a, "fname", exact = TRUE))){
    stop("adloglikelihood objects must have a fname attribute")
  }
  return(ll)
}
#' @export
print.adloglikelihood <- function(x, ...){
  print(sprintf("log-likelihood function '%s' for automatic differentiation", attr(x, "fname", exact = TRUE)))
}

#' @name customll
#' @returns `customll_test()` returns `TRUE` if taping works in your `R` session. Otherwise it will return `FALSE` and generate warnings.
#' @export
customll_test <- function(){
  myll <- customll("a1type test(const veca1 &u, const veca1 &beta) {
  size_t d  = u.size();
  a1type y(0.);  // initialize summation at 0
  y = u[0] * beta[0];
  return y;
  }")
  tape <- tapell(myll, 0.1, NA, 
                   manifoldtransform("sim", "identity", "sim")$tran,
                 thetatape_creator = function(n){c(1)})
  t1 <- (abs(drop(evaltape(tape, 2, 5)) - 2*5) < sqrt(.Machine$double.eps))
  if (!t1){
    warning("taping fails")
    return(FALSE)
  }

  Jtape <- tapeJacobian(tape)
  t2 <- (abs(drop(evaltape(Jtape, 2, 10)) - 10) < sqrt(.Machine$double.eps))
  if (!t2){
    warning("taping of Jacobian fails")
    return(FALSE)
  }

  return(TRUE)
}
