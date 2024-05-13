#' @noRd
#' @param ... passed directly to estimator for testing purposes

Windham_assess_estimator <- function(estimator, Y, ..., w = NULL){
  if (is.null(estimator)){stop("estimator is NULL")}
  estimatorformals <- formals(estimator)
  if (length(setdiff(c("Y", "w"), names(estimatorformals))) != 0){
    stop("Estimator must have arguments 'Y' and 'w'.")
  }

  if (is.null(w)){
    w <- rep(1, nrow(Y)) / sqrt(nrow(Y))
  }

  estargs <- list(...)
  estargs$w <- w
  estargs$Y <- Y
  estobj <- do.call(estimator, estargs)
  estlocation <- find_paramvec_location(estobj)

  newparamvec <- extract_paramvec(estobj)

  # check that newparamvec has the correct length, if the length can be determined
  paramvec_length = 0
  paramvec = FALSE 
  if ("paramvec" %in% names(estimatorformals)){
    paramvec = TRUE
    if (!is.null(estargs$paramvec)){paramvec_length <- max(paramvec_length, length(estargs$paramvec))} #use max here to ignore times when the default is 'NULL'
    else if (isa(estimatorformals$paramvec, "call")){paramvec_length <- max(paramvec_length, length(eval(estimatorformals$paramvec)))}
  }
  paramvec_start = FALSE 
  if ("paramvec_start" %in% names(estimatorformals)){
    paramvec_start = TRUE
    if (!is.null(estargs$paramvec_start)){paramvec_length <- max(paramvec_length, length(estargs$paramvec_start))}
    else if (isa(estimatorformals$paramvec_start, "call")){paramvec_length <- max(paramvec_length, length(eval(estimatorformals$paramvec_start)))}
  } 
  if ((paramvec_length > 0) && (length(newparamvec) != paramvec_length)){
    stop(sprintf("Returned estimate has different length (%i) to input paramvec or paramvec_start (%i)", length(newparamvec), paramvec_length))
  }

  #check that fixing works, but only when paramvec is passed
  if (paramvec){
    if (!is.null(estargs$paramvec)){ 
      if (any(abs(newparamvec[t_u2i(estargs$paramvec)] - estargs$paramvec[t_u2i(estargs$paramvec)]) > sqrt(.Machine$double.eps))){
        stop("Fixed elements of the parameter vector are altered by estimator.")
    }}
  }

  return(list(
    paramvec = paramvec,
    paramvec_start = paramvec_start,
    estlocation = estlocation,
    paramvec_length_tested = (paramvec_length > 0),
    est = newparamvec
  ))
}



