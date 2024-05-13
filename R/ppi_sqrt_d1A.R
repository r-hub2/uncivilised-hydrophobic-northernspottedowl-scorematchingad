
calcd1A <- function(p, sp, z, h, ind, qind, w = rep(1, nrow(z))){
	lambda4=4*(4+p-2)
	lambda2=2*(2+p-2)

	d_A=matrix(0,sp,1)
	for (j in 1:sp)
	{
		d_A[j]=weighted.mean((h^2)*(lambda4*z[,j]^4-12*z[,j]^2), w = w)
	}

	d_B=matrix(0,qind,1)
	for (j in 1:qind)
	{
		d_B[j]=2*weighted.mean((h^2)*(lambda4*z[,ind[1,j]]^2*z[,ind[2,j]]^2-2*z[,ind[1,j]]^2-2*z[,ind[2,j]]^2), w = w)
	}

	d_C=matrix(0,sp,1)
	for (j in 1:sp)
	{
		d_C[j]=weighted.mean((h^2)*(lambda2*z[,j]^2-2), w = w)
	}

	d=rbind(d_A,d_B,d_C)
  return(d)
}
