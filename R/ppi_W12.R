
calcd6_fixedbeta <- function(p, sp, z, h, ind, qind, beta0, w = rep(1, nrow(z))){
  Wnew <- calcW12(p, sp, z, h, ind, qind, w=w)
	pi2=1+2*beta0

	ev=-1*Wnew%*%pi2
	return(ev)
}

calcW12 <- function(p, sp, z, h, ind, qind, w = rep(1, nrow(z))){
	D1=matrix(0,1,sp)
	for (j in 1:sp)
	{
		D1[j]=weighted.mean(4*(h^2)*z[,j]^2, w=w)
	}
	AD1=diag(as.vector(D1))
	AD=matrix(0,sp,p)
	AD[1:sp,1:sp]=AD1



	BD=matrix(0,qind,p)
	for (i in 1:qind)
	{
		for (j in 1:p)
		{
			if (j==ind[1,i]){BD[i,j]=weighted.mean(4*h^2*z[,ind[2,i]]^2, w=w)}
			else if (j==ind[2,i]){BD[i,j]=weighted.mean(4*h^2*z[,ind[1,i]]^2, w=w)}
		}
	}


	D2=matrix(0,1,sp)
	for (j in 1:sp)
	{
		D2[j]=weighted.mean(2*(h^2), w = w)
	}
	CD0=diag(as.vector(D2))
	CD=matrix(0,sp,p)
	CD[1:sp,1:sp]=CD0


	AD2=AD
	for (i in 1:sp)
	{
		for (j in 1:p)
		{
			AD2[i,j]=4*weighted.mean((h^2)*z[,i]^4, w = w)
		}
	}

	BD2=BD
	for (i in 1:qind)
	{
		for (j in 1:p)
		{
			BD2[i,j]=8*weighted.mean((h^2)*(z[,ind[1,i]]^2*z[,ind[2,i]]^2), w=w)
		}
	}

	CD2=CD
	for (i in 1:sp)
	{
		for (j in 1:p)
		{
			CD2[i,j]=2*weighted.mean((h^2)*z[,i]^2, w=w)
		}
	}


	Wnew1=rbind(AD,BD,CD)

	Wnew2=rbind(AD2,BD2,CD2)

	Wnew=Wnew1-Wnew2
  return(Wnew)
}

