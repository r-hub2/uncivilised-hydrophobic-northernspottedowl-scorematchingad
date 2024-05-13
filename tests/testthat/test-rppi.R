# test that rppi old matches old rppi new

#sample size
n=5000

#parameters for the PPI model
m <- rppi_egmodel(2)

test_that("current PPI simulation method gives samples with similar empirical density estimates as the original simulation method", {
  skip_on_cran() #accuracy of the method is tested by all the estimators
  # simulate using old method
  time_historic <- system.time(samp2 <- rppi_singly(n,m$p,m$beta,m$AL,m$bL,4))
  if (!requireNamespace("ks")){skip("Need ks package")}
  H <- ks::Hpi(samp2$samp3[, -m$p])
  kde_historic <- ks::kde(samp2$samp3[, -m$p], H)
  #simulate sample from PPI model
  time_current <- system.time(samp3 <- rppi(n,beta = m$beta,AL=m$AL,bL=m$bL, maxden = 4))
  H <- ks::Hpi(samp3[, -m$p])
  kde_current <- ks::kde(samp3[, -m$p], H)

  testresult <- ks::kde.test(samp2$samp3[, -m$p], samp3[, -m$p])
  expect_gte(testresult$pvalue, 0.01)
  time_ratio <- time_current/time_historic
  expect_lt(time_ratio[["user.self"]], 0.05)
})

test_that("rppi() is fixed by set.seed()", {
  m <- rppi_egmodel(2)

  set.seed(3212)
  Y1 <- rppi(100, beta = m$beta, AL= m$AL, bL= m$bL, maxden=4)

  set.seed(3212)
  Y2 <- rppi(100, beta = m$beta, AL= m$AL, bL= m$bL, maxden=4)

  expect_equal(Y1, Y2)
})

test_that("rppi() passed a paramvec works", {

  set.seed(3212)
  Y1 <- rppi(100, beta= m$beta, AL= m$AL, bL = m$bL, maxden=4)
  
  set.seed(3212)
  Y2 <- rppi(100, paramvec = m$theta, maxden=4)
  
  expect_equal(Y1, Y2)
})

test_that("rppi() errors appropriately", {

  expect_error(rppi(100, beta = m$beta, AL = m$AL, bL = m$bL, paramvec = m$theta, maxden = 4))
  expect_error(rppi(100, beta = c(m$beta, m$beta), AL = m$AL, bL = m$bL, maxden = 4))
  expect_error(rppi(100, beta = m$beta, bL = m$bL, maxden = 4))
})

test_that("dppi() produces -Inf results outside simplex", {
  m <- rppi_egmodel(1)
  prop <- matrix(c(-1, 1, 1, 0.1, 0.1, 0.1, 0.8, 0.8, 0.7), ncol = 3, byrow = TRUE)
  suppressWarnings(logdens <- dppi(prop, beta = m$beta, AL = m$AL, bL = m$bL) )
  expect_equal(logdens, rep(-Inf, 3))
})

test_that("rppi() passed a zero value for bL works", {
  set.seed(1)
  expect_equal(rppi(1, AL = rsymmetricmatrix(5-1), beta = runif(5), bL = 0),
  c(0.1580301, 0.2347493, 0.2087618, 0.3172094, 0.0812494),
  ignore_attr = TRUE,
  tolerance = 1E-5)
})
