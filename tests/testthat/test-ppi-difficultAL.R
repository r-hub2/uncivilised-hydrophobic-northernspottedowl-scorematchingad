skip_on_cran()
test_that("full ppi estimates are mostly within 3 SE for difficult AL with large maxden, p = 5", {
  set.seed(1273)
  p = 5
  ALs <- exp(rsymmetricmatrix(p-1, -4, 4))
  bL <- rep(0, p-1)
  beta <- c(-0.7, -0.8, -0.3, 0, 0)
  set.seed(1111)
  suppressWarnings(prop <- rppi(500, beta = beta, AL=ALs, bL=bL, maxden=35)) #rppi_singly took 1005 seconds, rppi() took 13seconds
  #prop %>% as_tibble() %>% tidyr::pivot_longer(everything()) %>% ggplot() + facet_wrap(vars(name)) + geom_freqpoly(aes(x=value))

  est_cppad <- ppi(prop, ppi_paramvec(bL = bL), trans = "sqrt", bdryw = "minsq",
                         method = "closed",
                         acut = 0.01,
                         bdrythreshold = 1E-20,
                         control = list(tol = 1E-10, maxit = 1000))
  expect_gt(mean(abs(est_cppad$est$paramvec - ppi_paramvec(AL = ALs, bL = bL, beta = beta)) <= 3 * est_cppad$SE$paramvec), 0.85)

  #don't expect that the beta are within a fraction of the true values
  expect_error(expect_equal(est_cppad$est$beta, beta, tolerance = 1E-1))
})

test_that("full ppi estimates are within 3 SE of beta for difficult AL with large maxden, p = 3", {
  set.seed(12735)
  p = 3
  ALs <- exp(rsymmetricmatrix(p-1, -4, 4))
  bL <- rep(0, p-1)
  beta <- c(-0.7, -0.3, 0)
  set.seed(11112)
  suppressWarnings(prop <- rppi(1000, beta=beta, AL=ALs, bL=bL, maxden=20))
  expect_equal(colMeans(prop), c(0.1, 0.99, 0.1), tolerance = 0.5)

  est_cppad <- ppi(prop, ppi_paramvec(bL = bL), trans = "sqrt", bdryw = "minsq",
                         method = "closed",
                         acut = 0.01,
                         bdrythreshold = 1E-20,
                         control = list(tol = 1E-10))
  expect_absdiff_lte_v(est_cppad$est$AL, ALs, 3 * est_cppad$SE$AL)
  expect_absdiff_lte_v(est_cppad$est$beta, beta, 3 * est_cppad$SE$beta)
})
