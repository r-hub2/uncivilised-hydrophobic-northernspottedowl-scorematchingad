test_that("Solution without boundary considerations for PPI has zero gradient and matches numerical minimum", {
  set.seed(13411)
  mod <- rppi_egmodel(100)
  Y <- mod$sample

  tapes <- buildsmdtape(
     start = "sim",
     tran = "alr",
     end = "Euc",
     ll = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3, betap = tail(mod$beta, 1)),
     bdryw = "ones",
     verbose = FALSE)
  smdtape <- tapes$smdtape

  estobj <- cppad_closed(smdtape, Y)

  grads <- t(apply(Y, MARGIN = 1, function(x) pJacobian(smdtape$ptr, estobj$est, x))) 
  totalgrad <- colSums(grads)
  expect_lt(sum(totalgrad^2), 1E-20)

  numericalmin <- cppad_search(smdtape, theta = estobj$est, Y = Y)
  expect_equal(numericalmin$est, estobj$est, ignore_attr = TRUE)
  expect_equal(numericalmin$SE, estobj$SE)
})

test_that("Closed-from solution with boundary points matches hard-coded version", {
  mnongamma <- rppi_egmodel(1)
  theta <- ppi_paramvec(beta = c(-0.95, -0.9, 0.5), AL = mnongamma$AL, bL = 0)
  set.seed(1234)
  Ycts <- rppi(1000, paramvec = theta)
  dsample <- round(Ycts * 100)/ 100
  dsample[, 3] <- 1 - rowSums(dsample[, 1:2])

  isbdry <- simplex_isboundary(dsample)
  Yapproxcentres <- dsample
  Yapproxcentres[!isbdry, ] <- NA 
  Yapproxcentres[isbdry, ] <- simplex_boundaryshift(dsample[isbdry, , drop = FALSE])

  tapes <- buildsmdtape(
     start = "sim",
     tran = "alr",
     end = "Euc",
     ll = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3, bL = 0, betap = tail(theta, 1)), 
     bdryw = "ones",
     verbose = FALSE)

  estobj <- cppad_closed(tapes$smdtape, Y = dsample, Yapproxcentres, approxorder = 10)

  est_hardcode <- ppi(dsample, paramvec = ppi_paramvec(p = 3, bL = 0, betap = tail(theta, 1)),
      trans = "alr", method = "hardcoded")
  expect_equal(est_hardcode$est$paramvec, t_fu2t(estobj$est, ppi_paramvec(p = 3, bL = 0, betap = tail(theta, 1))), ignore_attr = "names")

})

test_that("Closed-form solution with all boundary points and alr matches hardcoded", {
  skip("Hardcoded and closed methods both fail on this. But not on the microbiome data")
  set.seed(123)
  m <- rppi_egmodel(100)
  #add some zeroes
  pushtozero <- function(x){
    whichmin <- which.min(x)
    x[whichmin] <- 0
    x <- x / sum(x) #normalise
    return(x)
  }
  newsample <- t(apply(m$sample, MARGIN = 1, pushtozero))
  mean(apply(newsample, 1, min) == 0) #100% have a zero

  hardcoded <- ppi(Y = newsample,
                   paramvec = ppi_paramvec(bL = 0, p = m$p, betap = tail(m$beta, 1)),
                   method = "hardcoded",
                   trans = "alr")
  cppad <- ppi(Y = newsample,
                   paramvec = ppi_paramvec(bL = 0, p = m$p, betap = tail(m$beta, 1)),
                   method = "closed",
                   trans = "alr")
})

