#### Preparing microbiome data #############################################

list2env(ppi_microbiomedata_cleaned_TCAP(), globalenv())

#a_c for h function:
acut=0.01

#### Test including b_L ####
test_that("estimator1 and SE is historically correct with b_L included (Scealy and Wood, 2022, Table 4)", {


  #calculate scoring estimate:
  estimator= estimator1(propreal,acut,1, beta0, computeSE = TRUE)
  estimate1=estimator$est$paramvec
  estimate1 <- estimate1[1:(length(estimate1) - p)]
  dim(estimate1) <- c(length(estimate1), 1)
  #check historically
  expect_snapshot_value(signif(estimate1, 8), style = "json2") #8 is the default number of digits for jsonlite::serializeJSON

  #check it matches cppad ppi()
  est_cppad <- ppi(Y = propreal, acut = acut,
                   method = "closed",
                   trans = "sqrt", bdryw = "minsq",
                   bdrythreshold = 1E-15, shiftsize = 1E-10,
                   paramvec = ppi_paramvec(beta = beta0))
  expect_equal(est_cppad$est$paramvec, estimator$est$paramvec, ignore_attr = TRUE)

  #estimate of W matrix
  W_est=estimator$info$W
  expect_snapshot_value(round(max(W_est), 8), style = "json2") #have to use round here because the json conversion doesn't necessarily show it in scientific notation
  expect_snapshot_value(signif(mean(W_est), 8), style = "json2")
  expect_snapshot_value(signif(which.max(W_est), 8), style = "json")

  #standard errors
  std1= estimator$SE$paramvec
  std1 <- std1[1:(length(std1) - p)]
  expect_snapshot_value(signif(std1, 8), style = "json2")

  #estimated parameters
  thetamats <- ppi_parammats(estimator$est$paramvec)
  AL <- thetamats$AL
  bL <- thetamats$bL
  dim(bL) <- c(length(bL), 1)
  #values in Table 4 in the article:
  expect_snapshot_value(signif(AL[upper.tri(AL, diag = TRUE)], 8), style = "json2")
  expect_snapshot_value(signif(bL, 8), style = "json2")
  expect_snapshot_value(signif(estimate1/std1, 8), style = "json2")
})


test_that("alr and cppad closed estimator for this data set are consistent", {
  #check alr estimators too
  est_hardcoded <- ppi(Y = propreal, method = "hardcoded",
                 trans = "alr", 
                 paramvec = ppi_paramvec(p = ncol(propreal), bL = 0, betap = tail(beta0, 1)))

  est_closed <- ppi(Y = propreal, method = "closed",
                 trans = "alr", 
                 paramvec = ppi_paramvec(p = ncol(propreal), bL = 0, betap = tail(beta0, 1)))
  expect_equal(est_hardcoded$est$paramvec, est_closed$est$paramvec, ignore_attr = "names")
})

#### Test omitting b_L ####
test_that("estimator1 and SE is historically correct with b_L ommitted (article table 3)", {

  #calculate scoring estimate:
  estimator=estimator1(propreal,acut,0, beta0, computeSE = TRUE)
  estimate1=estimator$est$paramvec
  #rearrange to historical ordering
  estimate1 <- estimate1[1:(length(estimate1) - p - (p-1))]
  dim(estimate1) <- c(length(estimate1), 1)
  expect_snapshot_value(signif(estimate1, 8), style = "json2")

  #estimate of W matrix
  W_est=estimator$info$W
  expect_snapshot_value(round(max(W_est), 8), style = "json2") #have to use round here because the json conversion doesn't necessarily show it in scientific notation
  expect_snapshot_value(signif(mean(W_est), 8), style = "json2")
  expect_snapshot_value(signif(which.max(W_est), 8), style = "json")

  #standard errors
  std1= estimator$SE$paramvec
  std1 <- std1[1:(length(estimate1))] #rearrange to combn ordering for historical comparison

  #estimated parameters
  thetamats <- ppi_parammats(estimator$est$paramvec)
  AL <- thetamats$AL

  #values in Table 3 in the article:
  expect_snapshot_value(signif(AL[upper.tri(AL, diag = TRUE)], 8), style = "json2")
  expect_snapshot_value(round(estimate1/std1, 8), style = "json2")

})

