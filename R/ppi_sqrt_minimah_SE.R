# @title Standard error estimates for Score2 estimates
# @description Estimates the standard errors of Score2 estimates. See (Scealy and Wood, 2021; Sec. 3.2).
# @param prop compositional data (n by p matrix)
# @param acut \eqn{a_c} for the Hyvärinen weight function \eqn{h}
# @param estimate1 the value of [estimator1()] estimates
# @param W_est The \eqn{\hat{W}}{W} matrix estimated by [estimator1()]
# @param incb if `incb=1` then \eqn{b_L} is estimated otherwise \eqn{b_L} is fixed at zero.
# @param beta0 The fixed \eqn{beta_0}{beta0}.
# @return A vector of standard errors corresponding to each entry of the estimate by [estimator1()].
# @export
estimator1SE <- function(prop,acut,estimate1,W_est,incb, beta0, w = rep(1, nrow(prop)))
{
  n<-nrow(prop)
  p<-ncol(prop)


	#response on sphere scale
	z=sqrt(prop)

	#h function
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



	AD_big=array(0, dim=c(sp, p, n))
	for (i in 1:n)
	{
		D1=matrix(0,1,sp)
		for (j in 1:sp)
		{
			D1[j]=4*(h[i]^2)*z[i,j]^2
		}
		AD1=diag(as.vector(D1))
		AD=matrix(0,sp,p)
		AD[1:sp,1:sp]=AD1
		AD_big[1:sp,1:p,i]=AD
	}


	BD_big=array(0, dim=c(qind, p, n))
	for (k in 1:n)
	{
		BD=matrix(0,qind,p)
		for (i in 1:qind)
		{
			for (j in 1:p)
			{
				if (j==ind[1,i]){BD[i,j]=4*h[k]^2*z[k,ind[2,i]]^2}
				else if (j==ind[2,i]){BD[i,j]=4*h[k]^2*z[k,ind[1,i]]^2}
			}
		}
		BD_big[1:qind,1:p,k]=BD
	}

	CD_big=array(0, dim=c(sp, p, n))
	for (k in 1:n)
	{
		D2=matrix(0,1,sp)
		for (j in 1:sp)
		{
			D2[j]=2*(h[k]^2)
		}
		CD0=diag(as.vector(D2))
		CD=matrix(0,sp,p)
		CD[1:sp,1:sp]=CD0
		CD_big[1:sp,1:p,k]=CD
	}


	AD2_big=array(0, dim=c(sp, p, n))
	AD2=AD
	for (k in 1:n)
	{
		for (i in 1:sp)
		{
			for (j in 1:p)
			{
				AD2[i,j]=4*((h[k]^2)*z[k,i]^4)
			}
		}
		AD2_big[1:sp,1:p,k]=AD2
	}



	BD2_big=array(0, dim=c(qind, p, n))
	BD2=BD
	for (k in 1:n)
	{
		for (i in 1:qind)
		{
			for (j in 1:p)
			{
				BD2[i,j]=8*(h[k]^2)*(z[k,ind[1,i]]^2*z[k,ind[2,i]]^2)
			}
		}
		BD2_big[1:qind,1:p,k]=BD2
	}


	CD2_big=array(0, dim=c(sp, p, n))
	CD2=CD
	for (k in 1:n)
	{
		for (i in 1:sp)
		{
			for (j in 1:p)
			{
				CD2[i,j]=2*(h[k]^2)*z[k,i]^2
			}
		}
		CD2_big[1:sp,1:p,k]=CD2
	}

	pi2=1+2*beta0
	evm=matrix(0,sum(sp,sp,qind),n)
	for (k in 1:n)
	{

		Wnew1=rbind(AD_big[1:sp,1:p,k],BD_big[1:qind,1:p,k],CD_big[1:sp,1:p,k])
		Wnew2=rbind(AD2_big[1:sp,1:p,k],BD2_big[1:qind,1:p,k],CD2_big[1:sp,1:p,k])
		Wnew=Wnew1-Wnew2


		evm[,k]=-1*Wnew%*%pi2
	}



	lambda4=4*(4+p-2)
	lambda2=2*(2+p-2)

	d_A=matrix(0,sp,n)
	for (j in 1:sp)
	{
		d_A[j,]=((h^2)*(lambda4*z[,j]^4-12*z[,j]^2))
	}

	d_B=matrix(0,qind,n)
	for (j in 1:qind)
	{
		d_B[j,]=2*((h^2)*(lambda4*z[,ind[1,j]]^2*z[,ind[2,j]]^2-2*z[,ind[1,j]]^2-2*z[,ind[2,j]]^2))
	}

	d_C=matrix(0,sp,n)
	for (j in 1:sp)
	{
		d_C[j,]=((h^2)*(lambda2*z[,j]^2-2))
	}

	dm=rbind(d_A,d_B,d_C)

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


	dvm=t(cbind(dv_A,dv_B,dv_C))

	dm=dm-dvm+evm


	AA_big=array(0, dim=c(sp, sp, n))
	for (k in 1:n)
	{
		A1=matrix(0,1,p-1)
		for (j in 1:sp)
		{
			A1[j]=(16*(h[k]^2)*z[k,j]^6)
		}
		AA=diag(as.vector(A1))
		AA_big[1:sp,1:sp,k]=AA
	}



	AB_big=array(0, dim=c(sp, qind, n))
	for (k in 1:n)
	{
		AB=matrix(0,sp,qind)
		for (i in 1:sp)
		{
			for (j in 1:qind)
			{
				if (i==ind[1,j]){AB[i,j]=(16*z[k,ind[1,j]]^4*z[k,ind[2,j]]^2*h[k]^2)}
				if (i==ind[2,j]){AB[i,j]=(16*z[k,ind[1,j]]^2*z[k,ind[2,j]]^4*h[k]^2)}
			}
		}
		AB_big[1:sp,1:qind,k]=AB
	}


	BB_big=array(0, dim=c(qind, qind, n))
	for (k in 1:n)
	{
		BB=matrix(0,qind,qind)
		for (i in 1:qind)
		{
			for (j in 1:qind)
			{
				if (ind[1,i]==ind[1,j] & ind[2,i]==ind[2,j]){BB[i,j]=(16*z[k,ind[1,j]]^4*z[k,ind[2,j]]^2*h[k]^2)+(16*z[k,ind[1,j]]^2*z[k,ind[2,j]]^4*h[k]^2)}
				else if (ind[1,i]==ind[1,j]){BB[i,j]=((h[k]^2)*16*z[k,ind[2,j]]^2*z[k,ind[2,i]]^2*z[k,ind[1,j]]^2)}
				else if (ind[2,i]==ind[2,j]){BB[i,j]=((h[k]^2)*16*z[k,ind[1,j]]^2*z[k,ind[1,i]]^2*z[k,ind[2,j]]^2)}
				else if (ind[1,i]==ind[2,j]){BB[i,j]=((h[k]^2)*16*z[k,ind[2,i]]^2*z[k,ind[1,j]]^2*z[k,ind[2,j]]^2)}
				else if (ind[2,i]==ind[1,j]){BB[i,j]=((h[k]^2)*16*z[k,ind[1,i]]^2*z[k,ind[2,j]]^2*z[k,ind[1,j]]^2)}

			}
		}
		BB_big[1:qind,1:qind,k]=BB
	}



	AC_big=array(0, dim=c(sp, sp, n))
	for (k in 1:n)
	{

		C1=matrix(0,1,p-1)
		for (j in 1:sp)
		{
			C1[j]=(8*(h[k]^2)*z[k,j]^4)
		}
		AC=diag(as.vector(C1))
		AC_big[1:sp,1:sp,k]=AC
	}




	BC_big=array(0, dim=c(qind, sp, n))
	for (k in 1:n)
	{

		BC=matrix(0,qind,sp)
		for (i in 1:qind)
		{
			for (j in 1:sp)
			{
				if (j==ind[1,i]){BC[i,j]=(h[k]^2*z[k,ind[2,i]]^2*z[k,ind[1,i]]^2*8)}
				else if (j==ind[2,i]){BC[i,j]=(h[k]^2*z[k,ind[1,i]]^2*z[k,ind[2,i]]^2*8)}
			}
		}
		BC_big[1:qind,1:sp,k]=BC
	}


	CC_big=array(0, dim=c(sp, sp, n))
	for (k in 1:n)
	{

		C2=matrix(0,1,p-1)
		for (j in 1:sp)
		{
			C2[j]=(4*(h[k]^2)*z[k,j]^2)
		}
		CC=diag(as.vector(C2))
		CC_big[1:sp,1:sp,k]=CC
	}


	AA2_big=array(0, dim=c(sp, sp, n))
	for (k in 1:n)
	{

		AA2=AA
		for (i in 1:sp)
		{
			for (j in 1:sp)
			{
				AA2[i,j]=16*((h[k]^2)*z[k,i]^4*z[k,j]^4)
			}
		}
		AA2_big[1:sp,1:sp,k]=AA2
	}



	AB2_big=array(0, dim=c(sp, qind, n))
	for (k in 1:n)
	{

		AB2=AB
		for (i in 1:sp)
		{
			for (j in 1:qind)
			{
				AB2[i,j]=32*((h[k]^2)*z[k,i]^4*z[k,ind[1,j]]^2*z[k,ind[2,j]]^2)
			}
		}
		AB2_big[1:sp,1:qind,k]=AB2
	}


	BB2_big=array(0, dim=c(qind, qind, n))
	for (k in 1:n)
	{
		BB2=BB
		for (i in 1:qind)
		{
			for (j in 1:qind)
			{
				BB2[i,j]=64*((h[k]^2)*(z[k,ind[1,i]]^2*z[k,ind[2,i]]^2)*(z[k,ind[1,j]]^2*z[k,ind[2,j]]^2))
			}
		}
		BB2_big[1:qind,1:qind,k]=BB2
	}



	BC2_big=array(0, dim=c(qind, sp, n))
	for (k in 1:n)
	{

		BC2=BC
		for (i in 1:qind)
		{
			for (j in 1:sp)
			{
				BC2[i,j]=16*((h[k]^2)*z[k,j]^2*(z[k,ind[1,i]]^2*z[k,ind[2,i]]^2))
			}
		}
		BC2_big[1:qind,1:sp,k]=BC2
	}


	AC2_big=array(0, dim=c(sp, sp, n))
	for (k in 1:n)
	{
		AC2=AC
		for (i in 1:sp)
		{
			for (j in 1:sp)
			{
				AC2[i,j]=8*((h[k]^2)*z[k,i]^4*z[k,j]^2)
			}
		}
		AC2_big[1:sp,1:sp,k]=AC2
	}



	CC2_big=array(0, dim=c(sp, sp, n))
	for (k in 1:n)
	{
		CC2=AC
		for (i in 1:sp)
		{
			for (j in 1:sp)
			{
				CC2[i,j]=4*((h[k]^2)*z[k,i]^2*z[k,j]^2)
			}
		}
		CC2_big[1:sp,1:sp,k]=CC2
	}


	diff=matrix(0,sum(sp,sp,qind),n)
	for (k in 1:n)
	{

		AA[,]=AA_big[,,k]
		AB[,]=AB_big[,,k]
		AC[,]=AC_big[,,k]

		AB[,]=AB_big[,,k]
		BB[,]=BB_big[,,k]
		BC[,]=BC_big[,,k]

		AC[,]=AC_big[,,k]
		BC[,]=BC_big[,,k]
		CC[,]=CC_big[,,k]


		W1=cbind(AA,AB,AC)
		W2=cbind(t(AB),BB,BC)
		W3=cbind(t(AC),t(BC),CC)
		W0=rbind(W1,W2,W3)

		AA2[,]=AA2_big[,,k]
		AB2[,]=AB2_big[,,k]
		AC2[,]=AC2_big[,,k]

		AB2[,]=AB2_big[,,k]
		BB2[,]=BB2_big[,,k]
		BC2[,]=BC2_big[,,k]

		AC2[,]=AC2_big[,,k]
		BC2[,]=BC2_big[,,k]
		CC2[,]=CC2_big[,,k]


		W12=cbind(AA2,AB2,AC2)
		W22=cbind(t(AB2),BB2,BC2)
		W32=cbind(t(AC2),t(BC2),CC2)
		W2=rbind(W12,W22,W32)

		W=W0-W2

		diff[1:num1,k]=t(t(dm[1:num1,k]))-W[1:num1,1:num1]%*%estimate1
	}
	Sig0=(diff[1:num1,]%*%(t(diff[1:num1,]) * w))/sum(w) #this is variability matrix
	Gam0=W_est[1:num1,1:num1] #W is the double derivative of the score matching discrepancy, and thus the gradient of the 'score', and so W_est is the sensitivity matrix.
	var0=solve(Gam0)%*%Sig0%*%solve(Gam0)/sum(w)


	std=sqrt(diag(var0))



	return(std)
}

