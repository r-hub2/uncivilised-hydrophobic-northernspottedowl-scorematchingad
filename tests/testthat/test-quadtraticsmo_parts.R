test_that("Hess + Offset match gradient for a PPI Example", {
  mod <- rppi_egmodel(100)
  Y <- mod$sample

  tapes <- buildsmdtape(
     start = "sim",
     tran = "alr",
     end = "Euc",
     ll = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3, betap=0.5), 
     verbose = FALSE)
  smdtape <- tapes$smdtape

  values <- quadratictape_parts(smdtape, Y)

  # expect results to match for gradient
  theta <- head(mod$theta, 7) #theta could be anything though
  gradorig <- t(apply(Y, MARGIN = 1, function(x){pJacobian(smdtape$ptr, theta, x)}))

  gradpoly <- lapply(1:nrow(values$offset), function(i){
    drop(matrix(values$Hessian[i, ], ncol = length(theta)) %*% theta + 
      t(values$offset[i, , drop = FALSE]))
  })
  gradpoly <- do.call(rbind, gradpoly)
  expect_equal(gradorig, gradpoly)

  # if W is symmetric then should equal the smd up to a constant wrt theta
  # and my expectation is that W is symmetric for ppi
  smdorig <- evaltape(smdtape, theta, Y)

  smdpoly <- lapply(1:nrow(values$offset), function(i){
    drop(0.5 * theta %*% matrix(values$Hessian[i, ], ncol = length(theta)) %*% theta + 
      values$offset[i, , drop = FALSE] %*% theta)
  })
  smdpoly <- unlist(smdpoly)
  constant <- smdorig-smdpoly
  
  #test constant by trying another theta
  theta2 <- theta+1
  smdorig2 <- evaltape(smdtape, theta2, Y)
  smdpoly2 <- lapply(1:nrow(values$offset), function(i){
    drop(0.5 * theta2 %*% matrix(values$Hessian[i, ], ncol = length(theta2)) %*% theta2 +
      values$offset[i, , drop = FALSE] %*% theta2)
  })
  smdpoly2 <- unlist(smdpoly2)
  expect_equal(smdorig2-smdpoly2, constant)
})

test_that("quadratictape_parts with approx centres is close to quadratic_parts for simplex interior points", {
  mod <- rppi_egmodel(100)
  Y <- mod$sample
  Ycen <- simplex_boundaryshift(Y)

  tapes <- buildsmdtape(
     start = "sim",
     tran = "alr",
     end = "Euc",
     ll = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3), 
     verbose = FALSE)
  smdtape <- tapes$smdtape
  
  valuesexact <- quadratictape_parts(smdtape, Y)
  valuesapprox <- quadratictape_parts(smdtape, Y, tcentres = Ycen, approxorder = 1)
  expect_equal(valuesexact, valuesapprox, tolerance = 1E-5)
  # but still an approximation
  expect_gt(max(abs(valuesexact$offset - valuesapprox$offset)), 1E-10)

  #expect higher order approximation to be closer
  valuesapprox2 <- quadratictape_parts(smdtape, Y, tcentres = Ycen, approxorder = 10)
  expect_lt(sum((valuesexact$offset - valuesapprox2$offset)^2), 
            sum((valuesexact$offset - valuesapprox$offset)^2))
})


