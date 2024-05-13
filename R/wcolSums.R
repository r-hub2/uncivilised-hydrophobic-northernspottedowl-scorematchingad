#' @noRd 
#' @title An internal function for performing weighted colSums
#' @description This is used a number of times in the code base. Weights each row of `mat` by the corresponding `w` element, if `w` passed. Then does `colSums()`.
#' @param mat A 2D array.
#' @param w Either NULL or a vector of weights equal to the number of rows of w
wcolSums <- function(mat, w=NULL){
  if (is.null(w)){
    out <- colSums(mat)
  } else {
    stopifnot(length(w) == nrow(mat))
    wmat <- mat*w
    out <- colSums(wmat)
  }
  return(out)
}
