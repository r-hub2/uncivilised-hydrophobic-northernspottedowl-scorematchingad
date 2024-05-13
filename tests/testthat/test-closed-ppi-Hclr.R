test_that("Fitting ppi via clr transform gets close estimates via sqrt", {
  set.seed(1234)
  model <- rppi_egmodel(100, maxden = 4)

  out <- ppi(Y = model$sample,
             paramvec = ppi_paramvec(p = 3, betap = -0.5),
             trans = "clr", bdrythreshold = 1E-5)
  expect_absdiff_lte_v(out$est$paramvec, model$theta, out$SE$paramvec * 3)
  out2 <- ppi(Y = model$sample,
             paramvec = ppi_paramvec(p = 3, betap = -0.5),
             trans = "sqrt", acut = 0.01, bdryw = "minsq")
  expect_absdiff_lte_v(out$est$paramvec, out2$est$paramvec, 
       out$SE$paramvec * c(1, 0.1, 0.5, 1, 0.1, 1, 1, 1))
})


test_that("Fitting ppi all parameters via clr transform can get close to true values for eg ppi, and is insensitive to bdrythreshold", {
  set.seed(12345)
  model <- rppi_egmodel(10000, maxden = 4)

  expect_warning({out <- ppi(Y = model$sample,
             trans = "clr", bdrythreshold = 1E-10, constrainbeta = TRUE)}, "Beta estimates.*-1") #default
  expect_warning({out2 <- ppi(Y = model$sample,
             trans = "clr", bdrythreshold = 1E-5, constrainbeta = TRUE)}, "Beta estimates.*-1")

  expect_equal(out$est$paramvec, out2$est$paramvec)
  expect_equal(out$SE$paramvec, out2$SE$paramvec)
  expect_equal(out$info$covar, out2$info$covar)

  expect_absdiff_lte_v(out$est$paramvec, model$theta, out$SE$paramvec * 3)
  # and that the SE are small
  expect_lte_v(abs(out$SE$paramvec), c(rep(20, 5), 0.05, 0.05, 20))
})


test_that("Fitting ppi smd values via clr transform is insensitive to bdrythreshold", {
  set.seed(12345)
  m <- rppi_egmodel(1000, maxden = 4)
  #add some zeroes
  pushtozero <- function(x){
    whichmin <- which.min(x)
    x[whichmin] <- 0
    x <- x / sum(x) #normalise
    return(x)
  }
  newsample <- rbind(
    t(apply(m$sample[1:50, ], MARGIN = 1, pushtozero)),
    m$sample[51:1000, ])
  out <- ppi_smvalues(Y = newsample, evalparam = m$theta,
             trans = "clr", bdrythreshold = 1E-10) #default
  out2 <- ppi_smvalues(Y = newsample, evalparam = m$theta,
             trans = "clr", bdrythreshold = 1E-5)

  expect_equal(out, out2)
})

# method of moments from Scealy and Wood 2022 supplementary
dirichmom <- function(X) {
  # Method of Moments estimates of the Dirichlet Distribution
  temp <- dim(X); n <- temp[1]; m <- temp[2]
  # X <- cbind(X,matrix(1-apply(X,1,sum)))
  mom <- apply(X,2,mean)*(mean(X[,1])-mean(X[,1]^2))/(mean(X[,1]^2) - ((mean(X[,1]))^2))
  return(mom - 1)
}

test_that("Fitting high concentration Dirichlet with 0 and 1E-5 bdrythreshold gets same answer and is close to moment estimator and true value", {
  set.seed(13511)
  beta0 <- c(-0.8, -0.1, 0) #the largest element is going to look a bit like beta>0 
  Y <- MCMCpack::rdirichlet(10000, beta0+1)
  out <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0),           
             trans = "clr", bdrythreshold = 0)
  out2 <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0),           
             trans = "clr", bdrythreshold = 1E-5)

  expect_equal(out$est$paramvec, out2$est$paramvec)
  expect_equal(out$SE$paramvec, out2$SE$paramvec)
  expect_equal(out$info$covar, out2$info$covar)
  

  expect_absdiff_lte_v(dirichmom(Y), out$est$beta, rep(1E-2, 3))
  expect_absdiff_lte_v(dirichmom(Y), out$est$beta, out$SE$beta)
  expect_absdiff_lte_v(out$est$beta, beta0, 3 * out$SE$beta)
})


test_that("High concentration Dirichlet with zeros estimated broadly similar to moment", {
  set.seed(1311)
  beta0 <- c(-0.8, -0.1, 0) #the largest element is going to look a bit like beta>0 
  Y <- MCMCpack::rdirichlet(10000, beta0+1)
  #add some zeroes
  pushtozero <- function(x){
    whichmin <- which.min(x)
    x[whichmin] <- 0
    x <- x / sum(x) #normalise
    return(x)
  }
  newY <- rbind(
    t(apply(Y[1:8000, ], MARGIN = 1, pushtozero)),
    Y[8001:10000, ])
  mean(apply(newY, 1, min) == 0) #80% have a zero

  out <- ppi(Y = newY,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0),           
             trans = "clr", bdrythreshold = 1E-10)
  expect_absdiff_lte_v(dirichmom(newY), out$est$beta, rep(3E-1, 3))
  expect_absdiff_lte_v(dirichmom(newY), out$est$beta, out$SE$beta * 20)
})


test_that("Hess + Offset match gradient for Hclr in interior", {
  set.seed(1245)
  mod <- rppi_egmodel(100)
  Y <- mod$sample
  Y <- rbind(Y, rep(1/3, 3))

  tapes <- buildsmdtape("sim","clr", "Hn111", "ppi", ytape = c(0.2, 0.3, 0.5), 
                        usertheta = ppi_paramvec(p = 3))

  # find boundary points and remove them - expect results pJacobian to give good results when no points need approximation
  isbdry <- simplex_isboundary(Y, 1E-5)
  Y <- Y[!isbdry, ]
  values <- quadratictape_parts(tapes$smdtape, Y)

  # expect results to match for gradient
  gradorig <- t(apply(Y, MARGIN = 1, function(x){pJacobian(tapes$smdtape$ptr, mod$theta, x)}))

  gradpoly <- lapply(1:nrow(values$offset), function(i){
    drop(matrix(values$Hessian[i, ], ncol = length(mod$theta)) %*% mod$theta + 
      t(values$offset[i, , drop = FALSE]))
  })
  gradpoly <- do.call(rbind, gradpoly)
  expect_equal(gradorig, gradpoly)
})


test_that("W is symmetric for ppi with clr, fitting all parameters", {
  # if W is symmetric then should equal the smd up to a constant wrt theta
  # and my expectation is that W is symmetric for ppi with Hclr
  set.seed(1245)
  mod <- rppi_egmodel(100)
  Y <- mod$sample
  Y <- rbind(Y, rep(1/3, 3))
  isbdry <- simplex_isboundary(Y, 1E-5)
  Y <- Y[!isbdry, ]

  usertheta = ppi_paramvec(p = 3, beta = mod$beta)
  ftheta <- t_ut2f(usertheta, mod$theta)
  tapes <- buildsmdtape("sim","clr", "Hn111", "ppi", ytape = c(0.2, 0.3, 0.5), 
                        usertheta = usertheta)

  values <- quadratictape_parts(tapes$smdtape, Y)

  smdorig <- evaltape(tapes$smdtape, pmat = Y, xmat = ftheta)

  smdpoly <- lapply(1:nrow(values$offset), function(i){
    drop(0.5 * ftheta %*% matrix(values$Hessian[i, ], ncol = length(ftheta)) %*% ftheta + 
      values$offset[i, , drop = FALSE] %*% ftheta)
  })
  smdpoly <- unlist(smdpoly)
  constant <- smdorig-smdpoly
  
  #test constant by trying another theta
  ftheta2 <- ftheta+1
  smdorig2 <- evaltape(tapes$smdtape, pmat = Y, xmat = ftheta2)
  smdpoly2 <- lapply(1:nrow(values$offset), function(i){
    drop(0.5 * ftheta2 %*% matrix(values$Hessian[i, ], ncol = length(ftheta2)) %*% ftheta2 +
      values$offset[i, , drop = FALSE] %*% ftheta2)
  })
  smdpoly2 <- unlist(smdpoly2)
  expect_equal(smdorig2-smdpoly2, constant)
})


