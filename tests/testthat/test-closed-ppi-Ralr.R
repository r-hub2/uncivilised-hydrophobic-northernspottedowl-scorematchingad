test_that("Fitting ppi via alr transform with fixed beta gets close to true values", {
  skip_on_cran()
  set.seed(1234)
  model <- rppi_egmodel(1000, maxden = 4)

  acut = 0.1 #not needed for Ralr per se, but code still expects it
  out <- ppi(Y = model$sample,
             paramvec = ppi_paramvec(beta = model$beta),
             trans = "alr")

  expect_absdiff_lte_v(out$est$paramvec, model$theta, out$SE$paramvec * 3)
})

test_that("Fitting ppi via alr inc all beta gets close to true values", {
  set.seed(111)
  model <- rppi_egmodel(1000, maxden = 4)

  out <- ppi(Y = model$sample,
             trans = "alr")

  expect_absdiff_lte_v(out$est$paramvec, model$theta, out$SE$paramvec * 3)
  # eestimating the -0.8 betas is poor
})
