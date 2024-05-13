test_that("Fitting ppi via ppi_alrsqrt_robust version gets closer with outliers", {
  set.seed(1111)
  model <- rppi_egmodel(100, maxden = 4)

  model$sample[1, ] <- c(0.5, 0.4, 0.1)
  model$sample[2, ] <- c(0.5, 0.4, 0.1)

  out <- ppi_alrsqrt(Y = model$sample,
             paramvec = ppi_paramvec(p = ncol(model$sample), betap = -0.5),
             acut = 0.01)

  out_robust <- ppi_alrsqrt_robust(Y = model$sample,
             cW = ppi_cW_auto(0.001, model$sample), #c(0.1, 0.01, 0.01, 0.1, 0.01, 0, 0, 0), # 
             paramvec = ppi_paramvec(p = ncol(model$sample), betap = -0.5),
             acut = 0.01,
             paramvec_start = model$theta
             )

  expect_absdiff_lte_v(out_robust$est$paramvec, model$theta, abs(out$est$paramvec - model$theta))
})


