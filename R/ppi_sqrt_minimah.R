# @title Score matching estimate of the PPI model using a minima-based Hyvärinen weight function
# @description Estimates \eqn{A_L} and \eqn{b_L} of the PPI model using score matching and a minimum-like Hyvärinen weight function. \eqn{\beta} is fixed.
# @param Y compositional data (each row is a sample, each column corresponds to a component)
# @param  acut \eqn{a_c} for the weighting function \eqn{h}.
# @param incb if `incb=1` then \eqn{b_L} is estimated, otherwise \eqn{b_L} is fixed at zero
# @param beta The (fixed) Dirichlet exponents (\eqn{\beta}) of the model.
# @param computeSE Computes the standard error using [estimator1SE()]
# @details The PPI model is given in equation 3 of (Scealy and Wood, 2021). The matrices \eqn{A_L} and \eqn{b_L} must be estimated.
# This function implements the score matching estimator,
# \deqn{\hat{W}^{-1}\hat{d},}{W^{-1}d,}
# using a minima-based Hyvärinen weight function
# \deqn{\tilde{h}(z)^2 = \min(z_1^2, z_2^2, ..., z_p^2, a_c^2).}{h(z)^2 = min(z1^2, z2^2, ..., zp^2, a_c^2),}
# where \eqn{z} is a point in the positive orthant of the p-dimensional unit sphere
# and \eqn{z_j}{zj} is the jth component of z.
# For details of the score matching estimator see equations 16 - 19 in (Scealy and Wood, 2021).
# If \eqn{a_c} is greater than or equal to 1 then this Hyvärinen weight function corresponds to (Scealy and Wood, 2021; eqn 11), if it is less than 1 then it corresponds to (Scealy and Wood, 2021; eqn 12).
# For more on the Hyvärinen weight (see equation 7 and Section 3.2 of (Scealy and Wood, 2021)).
# @return A vector of the estimates for individual entries in the matrices \eqn{A} and \eqn{b}, and the estimated \eqn{\hat{W}}{W}. The former first contains the diagonal of \eqn{A} (except the last entry that is always zero for identifiability in the PPI model), then the upper triangle of \eqn{A} without the last column (again for identifiability) and finally the elements of \eqn{b} (except the last element, which is always 0 due to identifiability also) if `incb=1`.

ppi_usertheta_estimator1_compatible_incb <- function(usertheta){
  p <- ppiltheta2p(length(usertheta))
  d_isfixed <- ppi_paramvec(p, AL = FALSE, bL = FALSE, beta = TRUE)
  isfixed <- t_u2i(usertheta)
  if (all(d_isfixed == isfixed)){return(TRUE)}
  else {return(FALSE)}
}

ppi_usertheta_estimator1_compatible_zerob <- function(usertheta){
  p <- ppiltheta2p(length(usertheta))
  d_isfixed <- ppi_paramvec(p, AL = FALSE, bL = TRUE, beta = TRUE)
  isfixed <- t_u2i(usertheta)
  if (!all(d_isfixed == isfixed)){return(FALSE)}

  mats <- ppi_parammats(usertheta)
  if (!all(mats$bL == 0)){return(FALSE)}
  else {return(TRUE)}
}

# @export
estimator1 <- function(Y,acut,incb, beta, w=rep(1, nrow(Y)), computeSE = FALSE)
{
  n <- nrow(Y) #number of samples
  p <- ncol(Y) #number of dimensions, although what happens when the beta need to be estimated?
	#response on sphere scale
	z=sqrt(Y)

	#h is the h function without taking the bound acut into account
	h=matrix(1,n,1)
	for (j in 1:p)
	{
		h=h*z[,j]
	}
	indh=matrix(0,n,1)
	# in the following indh is either 0, 1, ..,p for each row
  # 0 means the acut contraint hit and h2 is given acut
	# 1, .., p means the minimum h occured at that component
	h2=h
	for (j in 1:n)
	{
		indh[j]=1
		zmin=z[j,1]

		for (i in 2:p)
		{
			zmin_prev=zmin
			zmin=min(zmin,z[j,i])
			if (zmin_prev > zmin){indh[j]=i}
		}
		zmin_prev=zmin
		zmin=min(zmin,acut)
		if (zmin_prev > zmin){indh[j]=0}
		h2[j]=zmin

	}
	h=h2


	################### ##calculate W ##################
	sp <- p - 1
	ind_qind <- indexcombinations(sp)
	ind <- ind_qind$ind
	qind <- ind_qind$qind
	W <- calcW11(p, z, h, ind, qind, w=w)

	################### ##calculate d(6) ##################
	ev <- calcd6_fixedbeta(p, sp, z, h, ind, qind, beta, w=w)

	################### ##calculate d(1) ##################
	d <- calcd1A(p, sp, z, h, ind, qind, w=w)

	################### ##calculate d(2) ##################
	dv <- calcd2A_minimah(sp, n, z, ind, qind, indh, w=w)

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
          theta <- c(quartic_sphere, beta) 
        } else {
          theta <- c(quartic_sphere, rep(0, p-1), beta) 
        }

        # compute SE
        if (computeSE){
          SE <- try(estimator1SE(Y, acut, quartic_sphere, W_est, incb, beta, w))
          if (length(SE) > 1){
            if (incb==1){
              SE <- c(SE, 0*beta)
            } else {
              SE <- c(SE, rep(0, p-1), 0*beta) 
            }
            SE <- c(list(paramvec = SE), ppi_parammats(SE))
          }
        } else {
          SE <- "Not calculated."
        }

	return(list(est = c(list(paramvec=theta), ppi_parammats(theta)),
                    SE = SE,
                    info = list(W=W_est)))
}

