#' @noRd
#' @title Score matching estimates for the Dirichlet distribution
#' @description Score matching estimates for the Dirichlet distribution, which is the PPI with \eqn{A_L=0} and \eqn{b_L=0}. The parameters to be estimates are the \eqn{\beta}{beta}. Standard errors are not calculated.
#' @param Y Compositional data (n by p matrix)
#' @param acut \eqn{a_c} for \eqn{h} function.
#' @return A vector of estimates.
#' @describeIn estimator_dir The score matching estimator using the product-based Hyvärinen weight
#' \deqn{\tilde{h}(z)^2 = \min(\prod_{j=1}^{p} z_j^2, a_c^2).}{h(z)^2 = min(z1^2 * z2^2 * ... * zp^2, a_c^2).}
dir_sqrt_prodh <- function(Y,acut, w = rep(1, nrow(Y)))
{
  n=nrow(Y)
  p=ncol(Y)
	z=sqrt(Y)

	h=matrix(1,n,1)
	for (j in 1:p)
	{
		h=h*z[,j]
	}
	for (j in 1:n)
	{
		if (h[j] > acut) {h[j]=acut}
	}



	sp=p-1

        indqind <- indexcombinations(sp)
	ind=indqind$ind
	qind=indqind$qind

	h4m <- h2onz2_mean(p, n, z, h, indh = NULL, hstyle = "product", w = w)#note that indh
  W <- calcW22(p, h, h4m, w = w)

	d1=t(((p-2)*weighted.mean(h^2, w=w)+h4m))

	#ind2 indicates whether the acut constraint was hit
	ind2=matrix(1,n,1)
	for (j in 1:n)
	{
		if (h[j] == acut) {ind2[j]=0}
	}

  # the following is h4m as above, except that h4s is zero
	# whenever the acut constraint is hit
	# this is the first half of the equation for d(2)B at the end of
	# page 22 in the appendix
	#homit matrix. For a given component, k, homit is the multiple
	# of the z entries in the other components. I.e. prod(z_j)/z_k when z_k is not zero.
	homit=matrix(1,n,p)
	for (k in 1:p)
	{
		for (j in 1:p)
		{
			if (k != j){homit[,k]=homit[,k]*z[,j]}
		}

	}
	h4s=matrix(0,n,p)
	h4m=matrix(0,1,p)
	for (j in 1:p)
	{
		for (i in 1:n)
		{
			if (h[i] > 0){h4s[i,j]=ind2[i]*(h[i]^2)/z[i,j]^2}
			else {h4s[i,j]=ind2[i]*homit[i,j]^2}

		}
		h4m[j]=weighted.mean(h4s[,j], w=w)
	}



	d2=t(2*h4m-2*p*weighted.mean(ind2*h^2, w = w))

	d=d1-d2

	pp=solve(W)%*%d


	estimate2=(pp-1)/2

	return(estimate2)

}


