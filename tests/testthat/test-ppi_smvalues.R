test_that("PPI ALR hardcoded estimate has low smd and smgrad values and constant Hessian and offset", {
  mnongamma <- rppi_egmodel(1)
  theta <- ppi_paramvec(beta = c(-0.95, -0.9, 0.5), AL = mnongamma$AL, bL = 0)
  set.seed(1234)
  Ycts <- rppi(1000, paramvec = theta)
  dsample <- round(Ycts * 100)/ 100
  dsample[, 3] <- 1 - rowSums(dsample[, 1:2])
  colMeans(dsample == 0)
  mean(apply(dsample, 1, min) == 0)  #0.96

  usertheta = ppi_paramvec(p = ncol(dsample),
                           bL = 0,
                           betap = tail(theta, 1))
  est_hardcoded <- ppi(dsample, paramvec = usertheta, trans = "alr", method = "hardcoded")

  hardcodedvals <- ppi_smvalues(dsample, paramvec = usertheta, evalparam = est_hardcoded$est$paramvec, trans = "alr")

  modelvals <- ppi_smvalues(dsample, paramvec = usertheta, evalparam = theta, trans = "alr")

  # for any given sample the estimate would be better than the true value, and can test this
  expect_lt(hardcodedvals$obj, modelvals$obj)
  expect_lt_v(abs(hardcodedvals$grad), abs(modelvals$grad))
  expect_lt_v(abs(hardcodedvals$grad), rep(1E-10, length(hardcodedvals$grad)))

  # expect the Hessian and offset to be the same - because they are independent of theta, only depend on the data
  expect_equal(hardcodedvals$hess, modelvals$hess)
  expect_equal(hardcodedvals$offset, modelvals$offset)


  perobservation <- ppi_smvalues(dsample, paramvec = usertheta, evalparam = theta, trans = "alr", average = FALSE)
  expect_equal(lapply(perobservation, colMeans), modelvals)
})



