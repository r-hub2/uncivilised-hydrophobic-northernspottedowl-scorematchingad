####model in Section 2.3 with beta=(-0.8,-0.8,-0.5) here a_c is set to 20%####
# from 'modelA.R' in the original code

test_that("Score1ac estimator works on highly concentrated data, with some components close to the boundary", {
  #dimension
  p=3

  #sample size
  n=100

  #parameters for the PPI model
  muL=matrix(0,p-1,1)
  muL[1:sum(p,-1)]=0.12
  aa=matrix(1/500,p-1,1)
  D=diag(as.vector(aa))
  SigA=D
  SigA[1,1]=SigA[1,1]*2
  cor=0.5
  SigA[1,2]=cor*sqrt(SigA[1,1]*SigA[2,2])
  SigA[2,1]=SigA[1,2]
  ALs=-0.5*solve(SigA)
  bL=solve(SigA)%*%muL
  beta0=matrix(-0.8,p,1)
  beta0[p]=-0.5

  #simulate sample from PPI model
  set.seed(31654)
  samp3=rppi(n,beta=beta0,AL=ALs,bL=bL,maxden=4)

  ####Score1ac estimator##

  #a_c for h function:
  acut=1.6e-02

  #calculate scoring estimate for full model (only beta[p] fixed at -0.5):
  estimator=estimatorall1(samp3,acut,-0.5)
  estimate1all=estimator$estimator1
  # use CppAD for SE
  cppadest <- ppi(Y = samp3, paramvec = ppi_paramvec(p = p, betap = -0.5), trans = "sqrt", bdryw="minsq", acut = acut)
  expect_equal(cppadest$est$paramvec, c(estimate1all, -0.5), ignore_attr = "names")
  # use SE estimates as if beta0 was fixed at the estimate (not estimated)
  SE <- cppadest$SE$paramvec
  # within 2*SE
  expect_absdiff_lte_v(c(estimate1all, -0.5),  c(diag(ALs), ALs[upper.tri(ALs)], bL, beta0), 2*SE)
  #plus invented bounds for beta0 estimates for now
  expect_true(all(abs(beta0[-p] - estimate1all[6:7]) <= 2*3/sqrt(n)))

  #calculate scoring estimate with beta fixed at beta0:
  estimator=estimator1(samp3,acut,1, beta0, computeSE = TRUE)
  estimate1=estimator$est$paramvec
  std1=estimator$SE$paramvec
  # check
  #2*SE bounds
  expect_absdiff_lte_v(estimate1, c(diag(ALs), ALs[upper.tri(ALs)], bL, beta0), 2*std1)
})

