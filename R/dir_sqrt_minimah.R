#' @noRd
#' @describeIn estimator_dir The score matching estimator using the minima-based Hyvärinen weight function
#' \deqn{\tilde{h}(z)^2 = \min(z_1^2, z_2^2, ..., z_p^2, a_c^2).}{h(z)^2 = min(z1^2, z2^2, ..., zp^2, a_c^2).}
dir_sqrt_minimah <- function(Y,acut, w=rep(1, nrow(Y)))
{
  n=nrow(Y)
  p=ncol(Y)

	z=sqrt(Y)

	h=matrix(1,n,1)
	for (j in 1:p)
	{
		h=h*z[,j]
	}
	indh=matrix(0,n,1)
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

	sp=p-1

        indqind <- indexcombinations(sp)
	ind=indqind$ind
	qind=indqind$qind

  h4m <- h2onz2_mean(p, n, z, h, indh, w = w)
  W <- calcW22(p, h, h4m, w=w)

	d1=t(((p-2)*weighted.mean(h^2, w)+h4m))


	d2=matrix(0,n,p)

	for (i in 1:n)
	{
		for (j in 1:p)
		{
			if (indh[i]==j){d2[i,j]=2*(1-z[i,j]^2)}
			else if (indh[i]==0){d2[i,j]=0}
			else {d2[i,j]=-2*(z[i,indh[i]]^2)}
		}
	}


	d3=matrix(0,1,p)
	for (j in 1:p)
	{
		d3[j]=weighted.mean(d2[,j], w = w)
	}


	d=d1-t(d3)

	pp=solve(W)%*%d
	estimate1=(pp-1)/2

	return(estimate1)
}

# function that determins if a ppi usertheta matches a dirichlet model
ppi_usertheta_for_dir_sqrt_minimah <- function(usertheta){
  p <- ppiltheta2p(length(usertheta))
  d_utheta <- ppi_paramvec(p, AL=0, bL = 0)
  if (isTRUE(all(d_utheta[!is.na(d_utheta)] ==
             usertheta[!is.na(d_utheta)])) &&
      all(is.na(usertheta[is.na(d_utheta)]))){
        return(TRUE)
  } else {return(FALSE)}
}
