#' @title Quickly Generate a Symmetric Matrix for Testing and Examples
#' @description A simple function for generating a symmetric matrix for use in examples.
#' The diagonal, and upper-triangular elements of the matrix are simulated independently from a uniform distribution. The lower-triangle of the output matrix is copied from the upper-triangle.
#' These matrices __do not__ represent the full range of possible symmetric matrices.
#' @param p The desired dimension of the matrix
#' @param min The minimum of the uniform distribution.
#' @param max The maximum of the uniform distribution
#' @return A `p` x `p` symmetric matrix.
#' @examples
#' rsymmetricmatrix(5)
#' @export
rsymmetricmatrix <- function(p, min = 0, max = 1){
  A <- matrix(NA, ncol = p, nrow = p)
  A[upper.tri(A)] <- stats::runif(sum(upper.tri(A)), max = max, min = min)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- stats::runif(p, min = min, max = max)
  return(A)
}
