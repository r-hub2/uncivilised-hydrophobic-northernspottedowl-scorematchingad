####model in Section 2.3 with beta=(-0.8,-0.8,?) here a_c is set to 20%####
# edited from 'modelA.R' in the original code

#### Setup ####
list2env(rppi_egmodel(1), globalenv())
theta <- c(diag(AL), AL[upper.tri(AL)], bL)



#### Tests ####
test_that("Score1ac estimator estimates beta0[0] and other consistently with cppad version", {
  #sample size
  n=100
  beta0=matrix(-0.8,p,1)
  beta0[p]=-0.5

  #simulate sample from PPI model
  set.seed(321)
  samp3=rppi(n,beta=beta0,AL=AL,bL=bL,maxden=4)
  theta <- ppi_paramvec(p,
              AL = AL, bL = drop(bL), beta = drop(beta0))

  ####Score1ac estimator##

  #a_c for h function:
  acut=1.6e-02

  #calculate scoring estimate for full model
  estimator=estimatorall1(samp3,acut, betap = NULL)
  estimate1all=estimator$estimator1

  # Get SE from CppAD methods
  intheta <- ppi_paramvec(p)
  tapes <- buildsmdtape("sim", "sqrt", "sph", "ppi",
                        samp3[1, ], intheta,
                        bdryw = "minsq",
                        acut = acut)
  SE <- sqrt(diag(sme_estvar(tapes$smdtape, drop(estimate1all), samp3)))

  expect_absdiff_lte_v(estimate1all, theta, 3 * SE)

  # compare to ppi via cppad
  expect_lt(sum(smvalues_wsum(tapes$smdtape, samp3, drop(estimate1all))$grad^2), 1E-14)
  est2 <- ppi(samp3, trans = "sqrt", bdryw = "minsq", acut = acut, bdrythreshold = 1E-20, control = list(tol = 1E-12), method = "closed")
  expect_equal(est2$est$paramvec, drop(estimate1all), tolerance = 1E-2) #within 1% of each other roughly
  expect_absdiff_lte_v(drop(estimate1all), est2$est$paramvec, 2*est2$SE$paramvec)
})

test_that("Score1ac estimator can estimate beta0[1:(p-1)] for beta0[p] larger than -0.5", {
  #sample size
  n=100
  beta0=matrix(-0.8,p,1)
  beta0[p] = 5
  theta <- ppi_paramvec(p, AL = AL, bL = drop(bL), beta = drop(beta0))

  #simulate sample from PPI model
  set.seed(124)
  samp3=rppi(n,beta=beta0,AL=AL,bL=bL,maxden=4)

  ####Score1ac estimator##

  #a_c for h function:
  acut=1.6e-02

  #calculate scoring estimate for full model, beta[p] correctly fixed
  estimator=estimatorall1(samp3,acut, betap = beta0[p])
  estimate1all=estimator$estimator1

  # SE from cppad
  intheta <- ppi_paramvec(p, betap = beta0[p])
  tapes <- buildsmdtape("sim", "sqrt", "sph", "ppi",
                        samp3[1, ], intheta,
                        bdryw = "minsq",
                        acut = acut)
  SE <- sqrt(diag(sme_estvar(tapes$smdtape, drop(estimate1all), samp3)))

  #3*SE bounds
  expect_absdiff_lte_v(estimate1all, theta[is.na(intheta)], 3 * SE)
})

test_that("Score1ac estimator can estimate beta0[1:(p-1)] for beta0[p] large but misspecified", {
  #sample size
  n=100
  beta0=matrix(-0.8,p,1)
  beta0[p]= 5
  theta <- ppi_paramvec(p, AL = AL, bL = drop(bL), beta = drop(beta0))

  #simulate sample from PPI model
  samp3=rppi(n,beta=beta0,AL=AL,bL=bL,maxden=4)

  ####Score1ac estimator##

  #a_c for h function:
  acut=1.6e-02

  #calculate scoring estimate for full model, beta[p] correctly fixed
  estimator=estimatorall1(samp3,acut, betap = -0.5)
  estimate1all=estimator$estimator1

  # SE from cppad
  intheta <- ppi_paramvec(p, betap = -0.5)
  tapes <- buildsmdtape("sim", "sqrt", "sph", "ppi",
                        samp3[1, ], intheta,
                        bdryw = "minsq",
                        acut = acut)
  SE <- sqrt(diag(sme_estvar(tapes$smdtape, drop(estimate1all), samp3)))

  #3*SE bounds
  expect_absdiff_lte_v(estimate1all, theta[is.na(intheta)], 3 * SE)
})
