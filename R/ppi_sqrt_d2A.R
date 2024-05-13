#' @noRd
#' @description Calculated d(2)A for when the weight function h is a minimum version
calcd2A_minimah <- function(sp, n, z, ind, qind, indh, w = rep(1, nrow(z))){
	dv_A=matrix(0,n,sp)
	for (i in 1:n)
	{
		for (j in 1:sp)
		{
			if (indh[i]==j){dv_A[i,j]=8*z[i,j]^4*(1-z[i,j]^2)}
			else if (indh[i]==0){dv_A[i,j]=0}
			else {dv_A[i,j]=-8*(z[i,j]^4*z[i,indh[i]]^2)}
		}
	}
	dv_A_mean=matrix(0,1,sp)
	for (j in 1:sp)
	{
		dv_A_mean[j]=weighted.mean(dv_A[,j], w = w)
	}


	dv_B=matrix(0,n,qind)
	for (i in 1:n)
	{
		for (j in 1:qind)
		{
			if (indh[i]==ind[1,j]){dv_B[i,j]=8*z[i,ind[2,j]]^2*z[i,ind[1,j]]^2*(1-z[i,ind[1,j]]^2)-8*z[i,ind[1,j]]^4*z[i,ind[2,j]]^2}
			else if (indh[i]==ind[2,j]){dv_B[i,j]=8*z[i,ind[2,j]]^2*z[i,ind[1,j]]^2*(1-z[i,ind[2,j]]^2)-8*z[i,ind[1,j]]^2*z[i,ind[2,j]]^4}
			else if (indh[i]==0){dv_B[i,j]=0}
			else {dv_B[i,j]=-16*(z[i,ind[1,j]]^2*z[i,ind[2,j]]^2*z[i,indh[i]]^2)}
		}
	}
	dv_B_mean=matrix(0,1,qind)
	for (j in 1:qind)
	{
		dv_B_mean[j]=weighted.mean(dv_B[,j], w=w)
	}

	dv_C=matrix(0,n,sp)
	for (i in 1:n)
	{
		for (j in 1:sp)
		{
			if (indh[i]==j){dv_C[i,j]=4*z[i,j]^2*(1-z[i,j]^2)}
			else if (indh[i]==0){dv_C[i,j]=0}
			else {dv_C[i,j]=-4*(z[i,j]^2*z[i,indh[i]]^2)}
		}
	}
	dv_C_mean=matrix(0,1,sp)
	for (j in 1:sp)
	{
		dv_C_mean[j]=weighted.mean(dv_C[,j], w = w)
	}

	dv=-t(cbind(dv_A_mean,dv_B_mean,dv_C_mean))
  return(dv)
}
