# tools for converting symmetric matrices to vectors and back

tosmatrix <- function(theta){
  # l = d + d * (d-1)/2
  # 2l = 2d + d^2 - d
  # 0 = d^2 + d - 2l
  # d = (-1 pm sqrt(1-4*(-2l)))/2
  #   = (-1 + sqrt(1+8l))/2
  d <- (-1+sqrt(8 * length(theta) + 1)) / 2
  mat <- matrix(NA, d, d)
  diag(mat) <- theta[1:d]
  mat[upper.tri(mat)] <- theta[(d + 1):length(theta)]
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  return(mat)
}

fromsmatrix <- function(mat){
  stopifnot(isSymmetric.matrix(mat))
  out <- c(diag(mat), mat[upper.tri(mat)])
  
  # element names
  namdiag <- paste0(1:nrow(mat),",",1:nrow(mat))
  indx <- which(upper.tri(mat), arr.ind = TRUE)
  nam <- paste0(indx[, "row"], ",", indx[, "col"])
  names(out) <- c(namdiag, nam)
  return(out)
}

# return the array of row/column indices, and the length
# for the upper triangle
indexcombinations <- function(d){
  AL <- matrix(NA, d, d)
  ind <- t(which(upper.tri(AL), arr.ind=TRUE))
  qind=length(ind[1,])
  return(list(ind = ind, qind = qind))
}

