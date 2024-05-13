#helper function for testing: given a sample, simulate integer weights, then return weights and a version where each measurement is replicated by the weight size
virtualweights <- function(Y, sizefactor = 1.5){
  ind <- sample(1:nrow(Y), ceiling(sizefactor*nrow(Y)), replace = TRUE)
  weights <- rep(0, nrow(Y))
  weights[as.numeric(names(table(ind)))] <- table(ind)
  newsample <- Y[ind, ]
  return(list(
    newY = newsample,
    w = weights
  ))
}
