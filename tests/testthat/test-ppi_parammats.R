test_that("to then from PPI param vector returns the same params", {
  set.seed(124)
  p = 5
  ALs <- rsymmetricmatrix(p-1, 0, 10)
  bL <- runif(p-1)
  beta <- runif(p, -1, 0)
  paramvec <- ppi_paramvec(AL = ALs,  bL=bL, beta=beta)
  params2 <- ppi_parammats(paramvec)
  expect_equal(params2$AL, ALs)
  expect_equal(params2$bL, bL, ignore_attr = "names")
  expect_equal(params2$beta, beta, ignore_attr = "names")
})

test_that("ppiltheta2p() matches other operations", {
  calculated_p <- vapply(3:10, function(p) ppiltheta2p(ppithetalength(p)), FUN.VALUE = 3.0)
  expect_equal(calculated_p, 3:10)

  p = 7
  ALs <- rsymmetricmatrix(p-1, 0, 10)
  bL <- runif(p-1)
  beta <- runif(p, -1, 0)
  expect_equal(ppiltheta2p(length(ppi_paramvec(AL=ALs, bL=bL, beta=beta))), p)
})

test_that("indexcombinations() for vectorising matrices matches ppi_paramvec and upper.tri()", {
  #skip("indexcombinations() not yet trying match to/ppi_parammats and upper.tri() - see reordering test below")
  p = 5
  Amat <- matrix(NA, nrow = p-1, ncol = p-1)
  diag(Amat) <- seq(1, p-1)
  Amat[upper.tri(Amat)] <- (1:sum(upper.tri(Amat))) + p -1
  Amat[lower.tri(Amat)] <- t(Amat)[lower.tri(Amat)]

  indinfo <- indexcombinations(p - 1)

  #vs upper.tri()
  expect_equal(indinfo$ind, t(which(upper.tri(Amat), arr.ind = TRUE)),
               ignore_attr = TRUE)

  #vs toPPI
  ppiparam <- ppi_paramvec(AL=Amat, bL=rep(0, p-1), beta=rep(0, p))
  indexcombparam <- c(diag(Amat),
  Amat[t(indinfo$ind)], #each column is the own dimension, each row single element to extract
  rep(0, p + p -1))
  expect_equal(ppiparam, indexcombparam, ignore_attr = "names")
})

test_that("from PPI param vector order", {
  set.seed(124)
  p = 5
  bL <- runif(p-1)
  beta <- runif(p, -1, 0)
  paramats <- ppi_parammats(c(seq(1:((p-1) + (p-1) * (p-2)/2)), bL, beta))

  #alternative
  sp <- p-1
  x=c(1:sp)
  ind=combn(x, 2, FUN = NULL, simplify = TRUE)
  qind=length(ind[1,])
  mat=matrix(NA,sp,sp)
  vals <- seq(1:((p-1) + (p-1) * (p-2)/2))
  for (j in 1:sp)
  {
    mat[j,j]=vals[j]
  }
  for (j in 1:qind)
  {
    mat[ind[1,j],ind[2,j]]=vals[sum(j,sp)]
    mat[ind[2,j],ind[1,j]]=vals[sum(j,sp)]
  }

  expect_false(all(mat == paramats$AL))

})
