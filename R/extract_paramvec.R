#' @noRd
#' @title Try to extract estimated parameter vector from the result of an estimator
#' @param estobj The output from an estimator in this package.
#' @description If `estobj` is a list then looks first in the slot `est$paramvec` then the first element of the list. If the output is a numeric vector then it will use this. If none of these return a numeric vector then an error will be flagged.
extract_paramvec <- function(estobj){
  paramvecloc <- find_paramvec_location(estobj)
  extractor <- extract_paramvec_fun(paramvecloc)
  return(extractor(estobj))
}

# returns a FUNCTION
extract_paramvec_fun <- function(paramvec_location){
  fun <- switch(paramvec_location,
         "[['est']][['paramvec']]" = function(x){x$est$paramvec},
         "[[1]]" = function(x){x[[1]]},
         "[]" = function(x){x})
  return(fun)
}

find_paramvec_location <- function(estobj){
   if (is.list(estobj)){
      estparamvec <- try(estobj$est$paramvec, silent = TRUE)
      if (is.numeric(estparamvec) & is.vector(estparamvec)){return("[['est']][['paramvec']]")}
      estparamvec <- estobj[[1]]
      if (is.numeric(estparamvec) & is.vector(estparamvec)){return("[[1]]")}
    }
    estparamvec <- estobj
    if (is.numeric(estparamvec) & is.vector(estparamvec)){return("[]")}
    stop("Could not detect estimated parameter vector")
}


