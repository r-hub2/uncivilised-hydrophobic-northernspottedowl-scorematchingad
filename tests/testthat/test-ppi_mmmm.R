skip_on_cran()
# test of the 'model 1' described in section A.11 of the supplementary material.
# This is a gGamma model with β = (−0.80, −0.85, 0, −0.2, 0), b = 0, and A_L = Table 3. Beta and b considered fixed.

#### Simulate Model ####
#dimension
p=3

#sample size
n=92

#set seed
# set.seed(1)

#parameters for the PPI model:
ALs=matrix(0,p-1,p-1)
bL=matrix(0,p-1,1)
ALs[1,1]= -127480.0929
ALs[1,2]= 14068.39057
ALs[2,1]= ALs[1,2]
ALs[2,2]= -8191.17253
beta0=matrix(0,p,1)
beta0[1]=-0.80
beta0[2]=-0.85
beta0[3]=0

#simulate sample from the PPI model
samp3=rppi(n,beta=beta0,AL=ALs,bL=bL,maxden=0)

#### Estimate from Simulated Sample ####
test_that("ppi_mmmm gives numerical non-NA values", {
  #simulate sample from the multinomial PPI model:
  ni=matrix(2000,n,1)
  ni=as.vector(ni)
  x=matrix(0,n,p)
  for (j in 1:n)
  {
    x[j,]=rmultinom(1,ni[j],prob=samp3[j,])
  }

  mult=ppi_mmmm(x, ni, beta0)
  expect_length(mult, 3)
  expect_true(sum(mult) != 0)
})
