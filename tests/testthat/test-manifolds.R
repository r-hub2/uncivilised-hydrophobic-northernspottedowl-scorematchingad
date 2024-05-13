test_that("Manifold objects can be created, and member functions run", {
  mod <- Rcpp::Module("manifolds", PACKAGE="scorematchingad")
  Euc <- new(mod$man_ad, "Euc")
  expect_s4_class(Euc, "Rcpp_man_ad")
  z <- c(0.1, 0.5)
  expect_equal(dim(Euc$Pmatfun(z)), c(2,2))
  expect_type(Euc$Pmatfun(z), "double")
  expect_equal(dim(Euc$dPmatfun(z, 1)), c(2,2))
})

test_that("Transform objects can be created, and member functions run", {
  mod <- Rcpp::Module("manifolds", PACKAGE="scorematchingad")
  obj <- mod$transform_ad
  alr <- new(obj, "alr")
  expect_s4_class(alr, "Rcpp_transform_ad")
  u <- c(0.1, 0.3, 0.6)
  z <- alr$toM(u)
  expect_equal(alr$fromM(z), u)
  expect_true(is.finite(alr$logdetJfromM(z)))
})

mod <- Rcpp::Module("manifolds", PACKAGE="scorematchingad")

test_that("Sphere manifold object matches analytic results", {
  sph <- new(mod$man_ad, "sph")
  z <- c(0.2, 0.8, 0.3)
  z <- z/sqrt(sum(z^2))
  
  Pmatz <- sph$Pmatfun(z)
  expect_equal(t(z) %*% Pmatz %*% runif(3), 0, ignore_attr = TRUE)
  
  #check of first dPmatfun for first component
  z1 <- z[1]
  numgrad <- numericDeriv(quote(sph$Pmatfun(c(z1, z[2:3]))), theta = "z1")
  expect_equal(drop(attr(numgrad, "gradient")), sph$dPmatfun(z, 0))
})

test_that("sqrt transform matches other calculations", {
  sqrt <- new(mod$transform_ad, "sqrt")
  u <- c(0.1, 0.3, 0.6)
  z <- sqrt$toM(u)
  expect_equal(z, sqrt(u))
  expect_equal(sum(z^2), 1)
  expect_equal(sqrt$fromM(z), u)
  
  #check determinant using integration over a unitcube
  integrand <- function(zmat){#each column is a measurement
    xmat <- apply(zmat, MARGIN = 2, sqrt$fromM)
    inunitcube <- (colSums(xmat < 0) == 0) * (colSums(xmat > 1) == 0)
    Jdets <- exp(apply(zmat, MARGIN = 2, sqrt$logdetJfromM))
    return(matrix(Jdets * inunitcube, nrow = 1))
  }
  if (!requireNamespace("cubature")){skip("Need cubature package")}
  volumeviaM <- cubature::hcubature(
    f = integrand,
    lowerLimit = c(0, 0),
    upperLimit = c(1, 1),
    vectorInterface = TRUE,
    fDim = 1)
  expect_equal(volumeviaM$integral, 1, tolerance = 1E-5)
})

test_that("Simplex manifold object matches analytic results", {
  sim <- new(mod$man_ad, "sim")
  u <- c(0.1, 0.3, 0.6)
  expect_equal(sim$Pmatfun(u), diag(rep(1, 3)) - rep(1, 3) %*% t(rep(1, 3))/3)

  # check projection
  set.seed(1342)
  x <- runif(3)
  u2 <- sim$Pmatfun(u) %*% x
  expect_equal(t(u2) %*% rep(1, 3), 0, ignore_attr = TRUE)

  expect_equal(sim$dPmatfun(u, 2), matrix(0, nrow = 3, ncol = 3))
})


test_that("Identity transform matches analytic results", {
  iden <- new(mod$transform_ad, "identity")
  u <- c(0.1, 0.3, 0.6)
  expect_equal(u, iden$toM(u))
  expect_equal(u, iden$fromM(u))
  expect_equal(iden$logdetJfromM(u), log(1))
})

test_that("alr transform matches analytic results",{
  alr <- new(mod$transform_ad, "alr")
  u <- c(0.1, 0.3, 0.6)
  z <- alr$toM(u)
  expect_length(z, 2)
  expect_equal(alr$fromM(z), u)
  
  # check determinant
  integrand <- function(zmat){#each column is a measurement
    Jdets <- exp(apply(zmat, MARGIN = 2, alr$logdetJfromM))
    return(matrix(Jdets, nrow = 1))
  } 
  if (!requireNamespace("cubature")){skip("Need cubature package")}
  volume <- cubature::hcubature(
    f = integrand,
    lowerLimit = c(-1, -1) * 1E2,
    upperLimit = c(1, 1) * 1E2,
    vectorInterface = TRUE,
    fDim = 1)
  expect_equal(volume$integral, 0.5, tolerance = 1E-3) #0.5 seems to be the area of simplex under local coordinates (got 0.5 from integrating the simplex against local coordinates (e1, e2))
})


test_that("Euc manifold object matches analytic results",{
  Euc <- new(mod$man_ad, "Euc")
  z <- c(0.5, 3)
  expect_equal(Euc$Pmatfun(z), diag(rep(1, 2)))
  expect_equal(Euc$dPmatfun(z, 2), matrix(0, nrow = 2, ncol = 2))
})

test_that("Hn111 manifold passes analytic tests",{
  Hn111 <- new(mod$man_ad, "Hn111")
  # check projection matrix
  z <- c(0.1, 0.5, -0.6)
  set.seed(1342)
  x <- runif(3)
  x2 <- Hn111$Pmatfun(z) %*% x
  expect_equal(t(x2) %*% rep(1, 3), 0, ignore_attr = TRUE)

  expect_equal(Hn111$dPmatfun(z, 2), matrix(0, nrow = 3, ncol = 3))
})

test_that("clr transform passes analytic tests",{
  clr <- function(Y){
  logY <- log(Y)
    lgeommean <- rowSums(logY)/3
    return(logY - lgeommean)
  }
  clrinv <- function(Z){
    expZ <- exp(Z)
    return(expZ/rowSums(expZ))
  }
  
  ldetJfromM <- function(z){
    u <- as.vector(clrinv(matrix(z, nrow = 1)))
    sum(log(u)) + log(length(u))
  }

  clr_cpp <- new(mod$transform_ad, "clr")
  u <- c(0.1, 0.3, 0.6)
  z <- clr_cpp$toM(u)
  expect_equal(z, clr(matrix(u, nrow = 1)), ignore_attr = TRUE)
  expect_equal(rep(1, 3) %*% z, 0, ignore_attr = TRUE)
  expect_equal(clr_cpp$fromM(z), u)

  # test determinant
  expect_equal(ldetJfromM(matrix(z, nrow = 1)), clr_cpp$logdetJfromM(z))

  # by integration check determinant
  integrand <- function(znpmat){#each column is a measurement
    zmat = rbind(znpmat, 0-colSums(znpmat))
    Jdets <- exp(apply(zmat, MARGIN = 2, clr_cpp$logdetJfromM))
    return(matrix(Jdets, nrow = 1))
  } 
  if (!requireNamespace("cubature")){skip("Need cubature package")}
  volume <- cubature::hcubature(
    f = integrand,
    lowerLimit = c(-1, -1) * 1E1,
    upperLimit = c(1, 1) * 1E1,
    vectorInterface = TRUE,
    fDim = 1)
  expect_equal(volume$integral, 0.5, tolerance = 1E-3) #0.5 seems to be the area of simplex under local coordinates (got 0.5 from integrating the simplex against local coordinates (e1, e2))
})


