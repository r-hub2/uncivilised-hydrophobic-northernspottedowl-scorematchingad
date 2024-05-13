test_that("WindamRobust() via ppi_robust() gives correct params on simulated data, with two outliers. p=3", {
  skip_on_cran()#"only extra checks are the variable cW ones"
  set.seed(1273)
  m <- rppi_egmodel(1000, maxden = 4)
  outlier1 <- c(0.9, 0.9, 0.01)
  outlier1 <- outlier1/sum(outlier1)
  outlier2 <- c(0.9, 0.1, 0.01)
  outlier2 <- outlier2/sum(outlier2)
  m$sample <- rbind(m$sample, outlier1, outlier2)

  #check non-robust estimates, excluding the outliers
  est_simple <- ppi(m$sample[1:1000, ], acut=0.1, method = "closed", trans = "sqrt", bdryw = "minsq")
  #get non-robust estimates with the outliers
  est_simple_outlier <- ppi(m$sample, acut=0.1, method = "closed", trans = "sqrt", bdryw = "minsq")

  #calculate robust estimates
  suppressWarnings({est <- ppi_robust(Y = m$sample, 
           cW = 0.1 * ppi_paramvec(m$p, AL = TRUE, bL = FALSE, beta = FALSE),
           acut=0.1, method = "closed", trans = "sqrt", bdryw = "minsq")})

  # variable c, expect estimates to be different
  cW <- ppi_paramvec(m$p, AL = matrix(c(0.1, 1E-3, 1E-3, 0.1), nrow = 2, ncol = 2),
                                 bL = 0, beta = 0)
  suppressWarnings({est_varcW <-  ppi_robust(Y = m$sample, 
           cW = cW,
           acut=0.1, method = "closed", trans = "sqrt", bdryw = "minsq")})


  expect_equal(est$est$paramvec, est_simple$est$paramvec, tolerance = 0.1, ignore_attr = TRUE)
  #below checks that the non-robust estimate with outliers is much different to the robust estimate
  expect_error(expect_equal(est$est$paramvec, est_simple_outlier$est$paramvec, tolerance = 0.1, ignore_attr = TRUE))

  # expect that the different cW values would lead to different estimates
  expect_gt(mean(abs(est$est$paramvec - est_varcW$est$paramvec)), 10)
})

test_that("robust ppi() with Ralr transform gives correct params on simulated, no outlier, data. p=3", {

  set.seed(1273)
  p = 3
  ALs <- exp(rsymmetricmatrix(p-1, -4, 4))
  bL <- matrix(0, nrow = p-1)
  beta <- c(-0.8, -0.3, 0)
  set.seed(1345) #this seed generates samples that are representative-enough for estimatorlog_ratio() to give close estimates
  prop <- rppi(1000, beta=beta, AL=ALs, bL=bL, maxden=4)

  #check non-robust estimates
  est_unload <- ppi_alr_gengamma(prop, betap = beta[p], w = rep(1, nrow(prop)))
  # ppi_parammats(est_unload$ppi)$ALs #fairly terrible at the AL
  # ppi_parammats(est_unload$ppi)$beta #pretty good at beta

  #calculate robust estimates
  cW=0.001
  est1 = ppi_robust(Y = prop, paramvec = ppi_paramvec(p=3, bL = 0, betap = 0),
             method = "closed", trans = "alr",
             cW = ppi_cW_auto(cW, prop))
  expect_equal(est1$est$AL, ALs, tolerance = 1)
  expect_equal(est1$est$beta, beta, tolerance = 1E-1, ignore_attr = "names")
  rmse <- function(v1, v2){sqrt(mean((v1 - v2)^2))}
  # expect the non-robust estimate to be equal or poorer in accuracy:
  expect_gt(rmse(beta, est_unload$est$beta) + 1E-6, rmse(beta, est1$est$beta))
})

test_that("robust ppi gives correct params on simulated, no outlier, data. p = 5", {
  set.seed(1273)
  p = 5
  ALs <- rsymmetricmatrix(p-1, -4, 4)
  bL <- matrix(0, nrow = p-1)
  beta <- c(-0.7, -0.8, -0.3, 0, 0)
  set.seed(13456) #this seed generates samples that are representative-enough for estimatorlog_ratio() to give close estimates
  prop <- rppi(10000, beta=beta, AL=ALs, bL=bL, maxden=4)
  # prop %>% as_tibble() %>% tidyr::pivot_longer(everything()) %>% ggplot() + facet_wrap(vars(name)) + geom_freqpoly(aes(x=value))

  #calculate robust estimates
  cW=0.1
  est1=ppi_robust(Y = prop, paramvec = ppi_paramvec(bL = 0, betap = tail(beta, 1), p=5), cW = ppi_cW(cW, 1, 1, 1, 0, 0), trans = "alr", method = "closed")
  expect_equal(est1$est$AL, ALs, tolerance = 1E0)
  expect_equal(est1$est$beta, beta, tolerance = 1E-1, ignore_attr = "names")
})

test_that("Windham_assess_estimator works", {
  set.seed(3121)
  Y <- matrix(runif(10*5), nrow = 10, ncol = 5)
  w <- NULL

  assessment <- Windham_assess_estimator(function(Y, w = rep(1, nrow(Y))){
    return(colMeans(Y * w))
  }, Y = Y, w = NULL)
  expect_equal(assessment[1:4], list(paramvec = FALSE,
       paramvec_start = FALSE,
       estlocation = "[]",
       paramvec_length_tested = FALSE))

  assessment <- Windham_assess_estimator(
    function(Y, paramvec, w = rep(1, nrow(Y))){
      m <- colMeans(Y*w)
      m[t_u2i(paramvec)] <- paramvec[t_u2i(paramvec)]
      return(m)
    },
    Y = Y, w = NULL, paramvec = rep(NA, ncol(Y)))
  expect_equal(assessment[1:4], list(paramvec = TRUE,
       paramvec_start = FALSE,
       estlocation = "[]",
       paramvec_length_tested = TRUE))
  
  goodfun <- function(Y, paramvec, paramvec_start = NULL, w = rep(1, nrow(Y))){
    m <- colMeans(Y*w)
    m[t_u2i(paramvec)] <- paramvec[t_u2i(paramvec)]
    return(m)
  }
  assessment <- Windham_assess_estimator(goodfun, Y = Y, w = NULL, paramvec = rep(NA, ncol(Y)))
  expect_equal(assessment[1:4], list(paramvec = TRUE,
       paramvec_start = TRUE,
       estlocation = "[]",
       paramvec_length_tested = TRUE))

  assessment <- Windham_assess_estimator(goodfun, Y = Y, w = NULL, paramvec = NULL)
  expect_equal(assessment[1:4], list(paramvec = TRUE,
       paramvec_start = TRUE,
       estlocation = "[]",
       paramvec_length_tested = FALSE))


  # test the estimation location detection
  assessment <- Windham_assess_estimator(function(Y, paramvec = rep(NA, ncol(Y)), paramvec_start = NULL, w = rep(1, nrow(Y))){
    m <- goodfun(Y, paramvec, paramvec_start = paramvec_start, w = w)
    return(m)
  } , Y = Y, w = NULL)
  expect_equal(assessment$estlocation, "[]")
  
  assessment <- Windham_assess_estimator(function(Y, paramvec, paramvec_start = NULL, w = rep(1, nrow(Y))){
    m <- goodfun(Y, paramvec, paramvec_start = paramvec_start, w = w)
    return(list(est = list(paramvec = m)))
  }, Y = Y, w = NULL, paramvec = NULL)
  expect_equal(assessment$estlocation, "[['est']][['paramvec']]")
  
  assessment <- Windham_assess_estimator(function(Y, paramvec, paramvec_start = NULL, w = rep(1, nrow(Y))){
    m <- goodfun(Y, paramvec, paramvec_start = paramvec_start, w = w)
    return(list(est = m))
  }, Y = Y, w = NULL, paramvec = NULL)
  expect_equal(assessment$estlocation, "[[1]]")
})

test_that("ppi_robust() matches specialist Windham_alrgengamma() for simple data with p = 5", {
  set.seed(1222)
  p = 5
  ALs <- rsymmetricmatrix(p-1, -4, 4)
  bL <- matrix(0, nrow = p-1)
  beta <- c(-0.7, -0.8, -0.3, 0, 0)
  set.seed(13456) #this seed generates samples that are representative-enough for estimatorlog_ratio() to give close estimates
  prop <- rppi(100, beta=beta, AL=ALs, bL=bL, maxden=4)
  # prop %>% as_tibble() %>% tidyr::pivot_longer(everything()) %>% ggplot() + facet_wrap(vars(name)) + geom_freqpoly(aes(x=value))

  #calculate robust estimates with ppi_robust()
  cW=0.1
  est1=ppi_robust(Y = prop, paramvec = ppi_paramvec(bL = 0, betap = tail(beta, 1), p=5), cW = ppi_cW(cW, 1, 1, 1, 0, 0), trans = "alr", method = "closed")

  # use the specialist function
  est2 = ppi_robust_alrgengamma(Y = prop,
                    paramvec = ppi_paramvec(bL = 0, betap = tail(beta, 1), p=5),
                    cW = ppi_cW(cW, 1, 1, 1, 0, 0),
                    method = "hardcoded")

  expect_equal(est1$est, est2$est, tolerance = 1E-5, ignore_attr = "names")
})
