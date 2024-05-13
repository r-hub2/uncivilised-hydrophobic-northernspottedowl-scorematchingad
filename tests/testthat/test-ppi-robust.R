test_that("ppi including betaL with cW gives correct params on simulated data, with two outliers. p=3", {
  set.seed(1273)
  m <- rppi_egmodel(1000, maxden = 4)
  outlier1 <- c(0.9, 0.9, 0.01)
  outlier1 <- outlier1/sum(outlier1)
  outlier2 <- c(0.9, 0.1, 0.01)
  outlier2 <- outlier2/sum(outlier2)
  m$sample <- rbind(m$sample, outlier1, outlier2)

  #non-robust estimates
  est_norobust <- ppi(m$sample, ppi_paramvec(p=3, betap = tail(m$beta, 1)), acut=0.1, method = "hardcoded", trans = "sqrt", bdryw = "minsq")
  est_norobust2 <- ppi_robust(Y = m$sample, paramvec = ppi_paramvec(p=3, betap = tail(m$beta, 1)), acut=0.1, method = "hardcoded", trans = "sqrt", bdryw = "minsq", cW = ppi_cW_auto(0, m$sample))

  expect_equal(est_norobust2$est$paramvec, est_norobust$est$paramvec)

  #robust
  est_robust1 <- ppi_robust(Y = m$sample, paramvec = ppi_paramvec(p=3, betap = tail(m$beta, 1)), acut=0.1, method = "hardcoded", trans = "sqrt", bdryw = "minsq", cW = ppi_cW_auto(1E-1, m$sample))

  rmse <- function(v1, v2){sqrt(mean((v1 - v2)^2))}
  expect_gt(rmse(m$theta, est_norobust$est$paramvec),
            rmse(m$theta, est_robust1$est$paramvec))
})

test_that("Robustness runs for hardcoded and cppad methods", {
  set.seed(1273)
  m <- rppi_egmodel(50, maxden = 4)  # at 20 got singularities
  outlier1 <- c(0.9, 0.9, 0.01)
  outlier1 <- outlier1/sum(outlier1)
  outlier2 <- c(0.9, 0.1, 0.01)
  outlier2 <- outlier2/sum(outlier2)
  m$sample <- rbind(m$sample, outlier1, outlier2)

  #Ralr
  out <- ppi_robust(Y = m$sample, paramvec = ppi_paramvec(p=3, bL = 0, betap = -0.5), method = "hardcoded", trans = "alr", cW = ppi_cW_auto(1E-1, m$sample))
  expect_gt(out$info$fpevals, 1)

  #dir minsq : with AL=0 and bL=0 the default weights are 1, but customisation of cW should alter this
  out <- ppi_robust(Y = m$sample, paramvec = ppi_paramvec(p=3, AL = 0, bL = 0), method = "hardcoded",
             trans = "sqrt", acut = 0.1, bdryw = "minsq", cW = 1E-1 * ppi_paramvec(p=3, AL=0, bL=0, beta=1))
  expect_gt(out$info$fpevals, 1) #errors
  #dir prodsq
  out <- ppi_robust(Y = m$sample, paramvec = ppi_paramvec(p=3, AL = 0, bL = 0), method = "hardcoded",
             trans = "sqrt", acut = 0.1, bdryw = "prodsq", cW = 1E-1 * ppi_paramvec(p=3, AL=0, bL=0, beta=1))
  expect_gt(out$info$fpevals, 1) #errors

  # estimator1 bL=0
  out <- ppi_robust(Y = m$sample, paramvec = ppi_paramvec(bL = 0, beta = m$beta), acut=0.1, method = "hardcoded", trans = "sqrt", bdryw = "minsq", cW = ppi_cW_auto(1E-1, m$sample))
  expect_gt(out$info$fpevals, 1)
  # estimator1 bL=0 prodsq
  out <- ppi_robust(Y = m$sample, paramvec = ppi_paramvec(bL = 0, beta = m$beta),
             acut=0.1, method = "hardcoded", trans = "sqrt", bdryw = "prodsq",
             cW = ppi_cW_auto(1E-1, m$sample))
  expect_gt(out$info$fpevals, 1)

  # estimator1 bL!=0
  out <- ppi_robust(Y = m$sample, paramvec = ppi_paramvec(beta = m$beta), acut=0.1, method = "hardcoded", trans = "sqrt", bdryw = "minsq", cW = ppi_cW_auto(1E-1, m$sample))
  expect_gt(out$info$fpevals, 1)
  # estimator1 bL!=0 prodsq
  out <- ppi_robust(Y = m$sample, paramvec = ppi_paramvec(beta = m$beta),
             acut=0.1, method = "hardcoded", trans = "sqrt", bdryw = "prodsq", cW = ppi_cW_auto(1E-2, m$sample))
  expect_gt(out$info$fpevals, 1)
  # estimatorall1 betap fixed
  out <- ppi_robust(Y = m$sample, paramvec = ppi_paramvec(p=3, betap = tail(m$beta, 1)),
             acut=0.1, method = "hardcoded", trans = "sqrt", bdryw = "minsq",
             cW = ppi_cW_auto(1E-2, m$sample))
  expect_gt(out$info$fpevals, 1)

  # estimatorall1 betap fitted
  suppressWarnings({out <- ppi_robust(Y = m$sample,
             acut=0.1, method = "hardcoded", trans = "sqrt", bdryw = "minsq",
             cW = ppi_cW_auto(1E-2, m$sample), constrainbeta = TRUE)})
  expect_gt(out$info$fpevals, 1)

  # cppad - the default takes a long time
  out <- ppi_robust(Y = m$sample,
             acut=0.1, method = "closed", trans = "sqrt", bdryw = "minsq",
             control = list(tol = 1E-10, maxit = 100),
             fpcontrol = list(Method = "Simple", MaxIter = 100),
             cW = ppi_cW_auto(1E-2, m$sample))
  expect_gt(out$info$fpevals, 1)
})
