skip_on_cran() #too slow

##### Prepare First Data Set #####
list2env(ppi_microbiomedata_TCAP(), globalenv())

test_that("hardcoded alr estimator matches cppad calculations for Cyanobacteria/Chloroplast, Actinobacteria, Proteobacteria and pooled", {

# quick check of closed vs hardcoded
est_hardcoded=ppi(Y = propreal,
         method = "hardcoded", trans = "alr",
         paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0))
  hardcodedvals <- ppi_smvalues(propreal, 
                      paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                      evalparam = est_hardcoded$est$paramvec,
                      trans = "alr")
  expect_lt_v(hardcodedvals$grad, rep(1E-15, length(hardcodedvals$grad)))

  est_cppad=ppi(Y = propreal,
         method = "closed", trans = "alr",
         paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
         bdrythreshold = 1E-20, shiftsize = 1E-20,
         control = list(maxit = 1E5, tol = 1E-20 * n))
  expect_equal(est_hardcoded$est$paramvec, est_cppad$est$paramvec, ignore_attr = TRUE)
})

test_that("robust ppi via alr estimator matches historical results on dataset with Cyanobacteria/Chloroplast, Actinobacteria, Proteobacteria and pooled", {
# prepare for robustness
#initial values for robust estimators
ALs=matrix(0,p-1,p-1)
bL=matrix(0,p-1,1)
ALs[1,1]= -127480.0929
ALs[1,2]= 14068.39057
ALs[1,3]= 1782.261826
ALs[1,4]=  -240.076568
ALs[2,2]= -8191.17253
ALs[2,3]=  -8.002680
ALs[2,4]= 374.693979
ALs[3,3]= -46.638659
ALs[3,4]= 9.027633
ALs[4,4]= -39.208915
ALs[lower.tri(ALs)] <- t(ALs)[lower.tri(ALs)]
beta0=matrix(0,p,1)
beta0[1]=-0.80
beta0[2]=-0.85
beta0[3]=0
beta0[4]=-0.2
beta0[5]=0
bL_est=bL
ALs_est=ALs
beta0_est=beta0
sp=p-1

#try simulating to see required maxden
sim <- rppi(1000, beta= beta0, AL=ALs, bL=bL, maxden = 0)


#calculate robust estimates
cW=0.7
cWvec <- ppi_cW(cW, TRUE, TRUE, FALSE, FALSE, FALSE)
est1=ppi_robust(Y = propreal,
                cW = cWvec,
                method = "hardcoded", trans = "alr",
                paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                paramvec_start = ppi_paramvec(AL = ALs_est, bL = bL_est, beta = beta0_est))
#estimate of A_L:
expect_snapshot_value(signif(est1$est$AL,6), style = "json2")
#estimate of beta:
expect_snapshot_value(signif(est1$est$beta,6), style = "json2")


# check that restarting fp ends up in a quick finish
est1b=ppi_robust(Y = propreal,
                cW = cWvec,
                method = "hardcoded", trans = "alr",
                paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                paramvec_start = est1$est$paramvec,
                fpcontrol = list(Method = "Simple", ConvergenceMetricThreshold = 1E-9)) #slightly larger to beat est1
expect_equal(est1b$est$paramvec, est1$est$paramvec)

# check that ppi_robust with 'closed' method also works
est1c=ppi_robust(Y = propreal,
                cW = cWvec,
                method = "closed", trans = "alr",
                paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                paramvec_start =  ppi_paramvec(AL = ALs_est, bL = bL_est, beta = beta0_est)) 
  expect_equal(est1c$est$paramvec, est1$est$paramvec)

# check specialist ppi_robust_alrgengamma() matches results too
est1d=ppi_robust_alrgengamma(Y = propreal,
                cW = cWvec,
                method = "closed", 
                paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                paramvec_start =  ppi_paramvec(AL = ALs_est, bL = bL_est, beta = beta0_est)) 
  expect_equal(est1d$est$paramvec, est1$est$paramvec, tolerance = 1E-5)
})

#### Test Second Data Set ####

test_that("robust ppi via alr estimator matches historical results on dataset with Spirochates, Verrucomicrobia, Cyanobacteria/Chloroplast, TM7 and pooled", {

  list2env(ppi_microbiomedata_SVCTP(), globalenv())

  ##Estimation

  #initial values for robust estimators
  ALs=matrix(-10000,p-1,p-1)
  bL=matrix(0,p-1,1)
  beta0=matrix(0,p,1)
  beta0[1]=-0.80
  beta0[2]=-0.80
  beta0[3]=-0.80
  beta0[4]=-0.80
  bL_est=bL
  ALs_est=ALs
  beta0_est=beta0

  #calculate robust estimates
  cW=1.25
  est1=ppi_robust(Y = propreal,
                   cW = ppi_cW(cW, TRUE, TRUE, TRUE, TRUE, FALSE),
                   method = "hardcoded", trans = "alr",
                   paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                   paramvec_start = ppi_paramvec(AL = ALs_est, bL = bL_est, beta = beta0_est),
                   fpcontrol = list(Method = "Simple", MaxIter = 50))
  #estimate of A_L:
  expect_snapshot_value(signif(ppi_parammats(est1$est$paramvec)$AL,6), style = "json2")
  #estimate of beta:
  expect_snapshot_value(signif(ppi_parammats(est1$est$paramvec)$beta,6), style = "json2")

  # check that closed method works too
  est1b=ppi_robust(Y = propreal,
                   cW = ppi_cW(cW, TRUE, TRUE, TRUE, TRUE, FALSE),
                   method = "closed", trans = "alr",
                   paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                   paramvec_start = ppi_paramvec(AL = ALs_est, bL = bL_est, beta = beta0_est),
                   fpcontrol = list(Method = "Simple", MaxIter = 50))
  expect_equal(est1b$est$paramvec, est1$est$paramvec)

  # check specialist ppi_robust_alrgengamma() matches results too
  est1c=ppi_robust_alrgengamma(Y = propreal,
                cW = ppi_cW(cW, TRUE, TRUE, TRUE, TRUE, FALSE),
                method = "closed", 
                paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                paramvec_start =  ppi_paramvec(AL = ALs_est, bL = bL_est, beta = beta0_est)) 
  expect_equal(est1c$est$paramvec, est1$est$paramvec, tolerance = 1E-3)
})


test_that("A simulate then estimate of TCAP errors appropriately", {
  list2env(ppi_microbiomedata_TCAP(), globalenv())
  #initial values for robust estimators
  ALs=matrix(0,p-1,p-1)
  bL=matrix(0,p-1,1)
  ALs[1,1]= -127480.0929
  ALs[1,2]= 14068.39057
  ALs[1,3]= 1782.261826
  ALs[1,4]=  -240.076568
  ALs[2,2]= -8191.17253
  ALs[2,3]=  -8.002680
  ALs[2,4]= 374.693979
  ALs[3,3]= -46.638659
  ALs[3,4]= 9.027633
  ALs[4,4]= -39.208915
  ALs[lower.tri(ALs)] <- t(ALs)[lower.tri(ALs)]
  beta0 <- c(-0.80,-0.85,0,-0.2,0)

  #calculate robust estimates
  cW=0.7
  cWvec <- ppi_cW(cW, TRUE, TRUE, FALSE, FALSE, FALSE)
  est1=ppi_robust_alrgengamma(Y = propreal,
                cW = cWvec,
                method = "hardcoded",
                paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                paramvec_start = ppi_paramvec(AL = ALs, bL = bL, beta = beta0))

  #simulate from estimate
  set.seed(129)
  samp <- rppi(n = nrow(propreal), paramvec = est1$est$paramvec, maxden = 0)
  # convert latent variables to multinomial values
  mnprop <- t(apply(samp, MARGIN = 1, FUN = function(probs){rmultinom(1, 2000, prob = probs)}, simplify = TRUE)/2000)

  options("show.error.messages" = FALSE)

  expect_warning({est2 <- ppi_robust_alrgengamma(Y = mnprop,
                  cW = cWvec,
                  method = "closed",
                  paramvec = ppi_paramvec(p=p, bL = 0, betap = 0),
                  paramvec_start =  est1$est$paramvec)},
                  "extremely small")
  options("show.error.messages" = TRUE)
})

