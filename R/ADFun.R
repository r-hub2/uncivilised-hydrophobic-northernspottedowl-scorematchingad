# R6 ADFun class. Making it here because templating of the Cpp class means it is hard to expose via Rcpp
# Using R6 because it has modify in place semantics which match the ADFun objects
# following https://adv-r.hadley.nz/r6.html

#' @title A Class for Storing a CppAD Tape (ADFun) Object
#' @description
#' *This is a low level object useful for implementing score matching estimators.*
#' An `R6` class for storing a 'pointer' to a `CppAD` tape in `C++` (also called an `ADFun`) and associated information. 
#' Currently tools for modifying this information are not available in this package, however tools for creating new `ADFun` objects from an existing `ADFun` are available.
#' Typically an `ADFun` object will be created by [`buildsmdtape()`].
#' @field ptr A `Rcpp` external pointer to a `CppAD` `ADFun` object.
#' @param ptr A `Rcpp` external pointer to a `CppAD` `ADFun` object.
#' @field xtape The (numeric) vector of independent variable values used for taping.
#' @param xtape The (numeric) vector of independent variables used for taping.
#' @field dyntape The (numeric) vector of dynamic parameters used for taping.
#' @param dyntape The (numeric) vector of dynamic parameters used for taping.
#' @field usertheta A (numeric) vector of `NA` values and fixed values specifying the parameters of taped function that were considered dynamic parameters or fixed parameters respectively.
#' @param usertheta A (numeric) vector of `NA` values and fixed values specifying the inputs of the taped function that were considered independent variables or dynamic parameters respectively.
#' @field name An easy to read name for the taped function
#' @param name An easy to read name for the taped function
#' @inheritSection buildsmdtape Introduction to CppAD Tapes
#' @inheritSection buildsmdtape Warning: multiple CPU
#' @examples
#' tapes <- buildsmdtape(
#'   "sim", "sqrt", "sph",
#'   ll = "ppi",
#'   ytape =  rep(1/3, 3),
#'   usertheta = ppi_paramvec(p=3),
#'   bdryw = "minsq",
#'   acut = 0.01,
#'   verbose = FALSE)
#' tapes$smdtape$xtape
#' tapes$smdtape$dyntape
#' tapes$smdtape$name
#' tapes$smdtape$ptr
#' @export
ADFun <- R6::R6Class("ADFun",
  private = list( #private == not modifiable
    .ptr = NULL,
    .name = NULL,
    .xtape = vector("numeric"),
    .dyntape = vector("numeric"),
    .usertheta = vector("numeric")
  ),
  public = list(
    #' @description Create a new `ADFun` object from an external pointer to a `CppAD` `ADFun` object.
    initialize = function(ptr, name = NULL, xtape = vector("numeric"), dyntape = vector("numeric"), usertheta = rep(NA_real_, length(dyntape))){
      stopifnot(isa(ptr, "externalptr"))
      stopifnot(is.null(name) | isa(name, "character"))
      stopifnot(isa(xtape, "numeric"))
      stopifnot(isa(dyntape, "numeric"))
      stopifnot(isa(usertheta, "numeric"))
      private$.ptr <- ptr
      private$.name <- name
      private$.xtape <- xtape
      private$.dyntape <- dyntape
      private$.usertheta <- usertheta
    }
  ),
  active = list(
    ptr = function(value){if (missing(value)){private$.ptr} else {stop("`$ptr' is read only", call. = FALSE)}},
    name = function(value){if (missing(value)){private$.name} else {stop("`$name' is read only", call. = FALSE)}},
    xtape = function(value){if (missing(value)){private$.xtape} else {stop("`$xtape' is read only", call. = FALSE)}},
    dyntape = function(value){if (missing(value)){private$.dyntape} else {stop("`$dyntape' is read only", call. = FALSE)}},
    usertheta = function(value){if (missing(value)){private$.usertheta} else {stop("`$usertheta' is read only", call. = FALSE)}}
  ),
  cloneable = FALSE
)


