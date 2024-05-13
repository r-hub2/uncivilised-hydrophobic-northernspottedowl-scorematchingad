test_that("ppi_alr_gengamma matches CppAD closed method for constant weight, p = 3", {
  set.seed(1234)
  m <- rppi_egmodel(100, maxden = 4)

  est_cppad <- ppi(m$sample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta[3]), trans = "alr", method = "closed", bdryw = "ones")
  est_hardcoded <- ppi(m$sample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta[3]), trans = "alr", method = "hardcoded")

  expect_equal(est_hardcoded$est$paramvec, est_cppad$est$paramvec, ignore_attr = "names")
})

test_that("ppi_alr_gengamma matches CppAD closed method for constant weight and data with zeros, p = 3", {
  set.seed(1234)
  m <- rppi_egmodel(100, maxden = 4)
  dsample <- round(m$sample * 100)/100
  dsample[, 3] <- 1 - rowSums(dsample[, 1:2])

  est_hardcoded <- ppi(dsample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta[3]), trans = "alr", method = "hardcoded")
  est_cppad <- ppi(dsample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta[3]), trans = "alr", method = "closed", bdryw = "ones",
                         bdrythreshold = 1E-100) #1E-200 was too small for some reason - it produced NaN values

  expect_equal(est_hardcoded$est$paramvec, est_cppad$est$paramvec, ignore_attr = "names")
})

test_that("ppi_alr_gengamma matches CppAD closed method for constant weight, p = 5", {
  skip_on_cran() #slow but useful
  set.seed(1273)
  p = 5
  ALs <- rsymmetricmatrix(p-1, -4, 4)
  bL <- rep(0, p-1)
  beta <- c(-0.7, -0.8, -0.3, 0, 0)
  # set.seed(1345) #this seed leads to samples with that give reasonable estimates
  set.seed(1111) #this seed leads to some ginormous elements for the second diagonal element of ALs
  prop <- rppi(1000, beta=beta, AL=ALs, bL=bL, maxden=5) #rppi_singly took 1005 seconds, rppi() took 13seconds

  est_cppad <- ppi(prop, ppi_paramvec(bL = bL, betap = beta[p]), trans = "alr", method = "closed", bdryw = "ones",
                         bdrythreshold = 1E-20)
  expect_absdiff_lte_v(est_cppad$est$AL, ALs, 3 * est_cppad$SE$AL)
  expect_absdiff_lte_v(est_cppad$est$beta, beta, 3 * est_cppad$SE$beta)
  #expect that the SE are small relative to size of the coefficients
  expect_lt(median(abs(est_cppad$SE$paramvec/est_cppad$est$paramvec), na.rm = TRUE), 0.3)

  est_hardcoded <- ppi_alr_gengamma(prop, betap = beta[p], w = rep(1, nrow(prop)))

  # check that estimates via cppad are close to hardcoded
  expect_equal(est_hardcoded$est$paramvec, est_cppad$est$paramvec, ignore_attr = TRUE)
})


test_that("ppi_alr_gengamma matches for simulated weights", {
  set.seed(1234)
  m <- rppi_egmodel(100, maxden = 4)
  #simulate weights
  ind <- sample(1:100, 150, replace = TRUE)
  weights <- rep(0, 100)
  weights[as.numeric(names(table(ind)))] <- table(ind)
  newsample <- m$sample[ind, ]

  est_sim <- ppi_alr_gengamma(newsample, betap = m$beta[3], w = rep(1, nrow(newsample)))
  est_hardcoded <- ppi_alr_gengamma(m$sample, betap = m$beta[3], w = weights)
  expect_equal(est_hardcoded$est$paramvec, est_sim$est$paramvec)
})

test_that("ppi_alr_gengamma() and cppad match for a randomly selected weight vector", {
  set.seed(1234)
  m <- rppi_egmodel(100, maxden = 4)
  set.seed(1212)
  w <- runif(100)

  est_cppad <- ppi(m$sample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta[3]), trans = "alr", method = "closed", bdryw = "ones",
                         w = w,
                         )
  est_hardcoded <- ppi(m$sample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta[3]), trans = "alr", method = "hardcoded", w = w)

  expect_equal(est_hardcoded$est$paramvec, est_cppad$est$paramvec, ignore_attr = "names")
})


