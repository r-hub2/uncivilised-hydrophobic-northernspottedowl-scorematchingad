# @title Score matching estimate of the PPI model using product-like Hyvärinen weight function
# @description Estimates \eqn{A_L} and \eqn{b_L} of the PPI model using score matching and a product-like Hyvärinen weight function. \eqn{\beta_0}{beta0} is fixed.
# @param prop compositional data (each row is a sample, each column corresponds to a component)
# @param  acut \eqn{a_c} for the weighting function \eqn{h}.
# @param incb if `incb=1` then \eqn{b_L} is estimated, otherwise \eqn{b_L} is fixed at zero
# @param beta0 The (fixed) beta0 of the model.
# @details The PPI model is given in equation 3 of (Scealy and Wood, 2021). The matrices \eqn{A_L} and \eqn{b_L} must be estimated.
# This function implements the score matching estimator,
# \deqn{\hat{W}^{-1}\hat{d},}{W^{-1}d,}
# using a product-like Hyvärinen weight function
# \deqn{\tilde{h}(z)^2 = \min(\prod_{j=1}^{p} z_j^2, a_c^2).}{h(z)^2 = min(z1^2 * z2^2 * ... * zp^2, a_c^2),}
# where \eqn{z} is a point in the positive orthant of the p-dimensional unit sphere
# and \eqn{z_j}{zj} is the jth component of z.
# For details of the score matching estimator see equations 16 - 19 in (Scealy and Wood, 2021).
# If \eqn{a_c} is greater than or equal to 1 then this product-like Hyvärinen weight function corresponds to (Scealy and Wood, 2021; eqn 8), if it is less than 1 then it corresponds to (Scealy and Wood, 2021; eqn 9).
# For more on the Hyvärinen weight (see equation 7 and Section 3.2 of (Scealy and Wood, 2021)).
# @return A vector of the estimates for individual entries in the matrices \eqn{A} and \eqn{b}, and the estimated \eqn{\hat{W}}{W}. The former first contains the diagonal of \eqn{A} (except the last entry that is always zero for identifiability in the PPI model), then the upper triangle of \eqn{A} without the last column (again for identifiability) and finally the elements of \eqn{b} (except the last element, which is always 0 due to identifiability also) if `incb=1`.

# @export
estimator2 <- function(prop,acut,incb, beta0, w = rep(1, nrow(prop)))
{
  n=nrow(prop)
  p=ncol(prop)

	#response on sphere scale
	z=sqrt(prop)

	#h function
	h=matrix(1,n,1)
	for (j in 1:p)
	{
		h=h*z[,j]
	}
	for (j in 1:n)
	{
		if (h[j] > acut){h[j]=acut}
	}

	################### ##calculate W ##################
	sp <- p - 1
	ind_qind <- indexcombinations(sp)
	ind <- ind_qind$ind
	qind <- ind_qind$qind
  W <- calcW11(p, z, h, ind, qind, w = w)
	################### ##calculate d(6) ##################
  ev <- calcd6_fixedbeta(p, sp, z, h, ind, qind, beta0, w=w)
	################### ##calculate d(1) ##################
  d <- calcd1A(p, sp, z, h, ind, qind, w=w)
	################### ##calculate d(2) ##################

	ind2=matrix(1,n,1)
	for (j in 1:n)
	{
		if (h[j] == acut) {ind2[j]=0}
	}

	dv_A=matrix(0,sp,1)
	for (j in 1:sp)
	{
		dv_A[j]=weighted.mean(-2*4*ind2*(h^2)*(z[,j]^2*(1-p*z[,j]^2)), w=w)
	}

	dv_B=matrix(0,qind,1)
	for (j in 1:qind)
	{
		dv_B[j]=2*weighted.mean((h^2)*ind2*(8*p*z[,ind[1,j]]^2*z[,ind[2,j]]^2-4*z[,ind[1,j]]^2-4*z[,ind[2,j]]^2), w=w)
	}

	dv_C=matrix(0,sp,1)
	for (j in 1:sp)
	{
		dv_C[j]=weighted.mean(4*ind2*(h^2)*(p*z[,j]^2-1), w=w)
	}

	dv=rbind(dv_A,dv_B,dv_C)

	################### ##calculate d total and scoring estimate ##################

	d=d+dv+ev


	if (incb==1)
	{
		#include bL in the model
		num1=sp+qind+sp
	}
	else
	{
		#omit bL from the model
		num1=sp+qind
	}





	#save W
	W_est=W

	#scoring estimator
	quartic_sphere=solve(W[1:num1,1:num1])%*%t(t(d[1:num1]))

        #convert to PPI param vec
        if (incb==1){
          theta <- c(quartic_sphere, beta0)
        } else {
          theta <- c(quartic_sphere, rep(0, p-1), beta0)
        }
	return(list(estimator2=quartic_sphere,W_est=W_est, theta = theta))
}

# cleaner PPI interface for estimator2
ppi_sqrt_prodh_zerob <- function(Y, acut, beta, w){
  rawout <- estimator2(Y,acut = acut,incb = 0,
                       beta0 = beta,
                       w= w)
  rawSE <- estimator2SE(prop = Y, acut = acut,
                        estimate2 = rawout$estimator2,
                        W_est = rawout$W_est,
                        incb = 0,
                        beta0 = beta,
                        w = w
  )
  fit <- list(est = list(
    paramvec = ppi_paramvec(AL = tosmatrix(rawout$estimator2),
                             bL = 0,
                             beta = beta)
  ))
  fit$est <- c(fit$est, ppi_parammats(fit$est$paramvec))
  fit$SE <- list(paramvec = ppi_paramvec(AL = tosmatrix(rawSE),
                                              bL = 0, beta = 0))
  fit$SE <- c(fit$SE, ppi_parammats(fit$SE$paramvec))
  return(fit)
}

ppi_sqrt_prodh <- function(Y, acut, beta, w){
  rawout <- estimator2(Y,acut = acut,incb = 1,
                       beta0 = beta,
                       w= w)
  rawSE <- estimator2SE(prop = Y, acut = acut,
                        estimate2 = rawout$estimator2,
                        W_est = rawout$W_est,
                        incb = 1,
                        beta0 = beta,
                        w = w
  )
  estparamvec <- c(rawout$estimator2, beta)
  fit <- list()
  fit$est <- c(list(paramvec = estparamvec),
               ppi_parammats(estparamvec))
  SEparamvec <- c(rawSE, 0 * beta)
  fit$SE <- c(list(paramvec = SEparamvec),
                 ppi_parammats(SEparamvec))
  fit$SE <- c(fit$SE, ppi_parammats(fit$SE$paramvec))
  return(fit)
}
