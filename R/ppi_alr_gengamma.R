# @title Score matching estimate of the generalised gamma form of the PPI model(b_L = 0) using an alr transformation
# @description Estimates \eqn{A_L} and \eqn{\beta_0}{beta0} of the PPI model using score matching in Euclidean space. Data is transformed to the Euclidean space using the additive log-ration transformation.
# The final element of \eqn{\beta_0}{beta0} is fixed at a chosen value.
# For estimating the elements of \eqn{A_L} more precisely, we recommend the score-matching estimators that transform the simplex on the sphere.
# @param Y compositional data (each row is a sample, each column corresponds to a component)
# @param betap The (fixed) final element of beta0 for the model.
# @param w A vector of weights to apply to the measurements in `Y`. The length of `w` must equal the number of rows of `Y`.
# @export
ppi_alr_gengamma <- function(Y,betap,w)
{
        p=ncol(Y)
	sp=p-1
        n=nrow(Y)

        indqind <- indexcombinations(sp)
        ind <- indqind$ind
        qind <- indqind$qind

	V1S_sum=matrix(0,sum(sp,sp,qind),1)
	V2S_sum=matrix(0,sum(sp,sp,qind),1)
	V3S_sum=matrix(0,sum(sp,sp,qind),sum(sp,sp,qind))


	for (b in 1:n)
	{

		V1A=matrix(0,sp,sp)
		for (j in 1:sp)
		{
			for (i in 1:sp)
			{
				V1A[i,j]=-2*(Y[b,i]^2)*Y[b,j]+6*(Y[b,i]^2)*(Y[b,j]^2)
				V1A[j,i]=-2*(Y[b,j]^2)*Y[b,i]+6*(Y[b,j]^2)*(Y[b,i]^2)
			}
			V1A[j,j]=4*Y[b,j]^2-10*Y[b,j]^3+6*Y[b,j]^4
		}


		V1B=matrix(0,qind,sp)
		for (j in 1:qind)
		{

			for (k in 1:sp)
			{
				if (ind[1,j]==k){V1B[j,k]=2*Y[b,ind[1,j]]*Y[b,ind[2,j]]-12*(Y[b,ind[1,j]]^2)*Y[b,ind[2,j]]+12*(Y[b,ind[1,j]]^3)*Y[b,ind[2,j]]}
				else if (ind[2,j]==k){V1B[j,k]=2*Y[b,ind[2,j]]*Y[b,ind[1,j]]-12*(Y[b,ind[2,j]]^2)*Y[b,ind[1,j]]+12*(Y[b,ind[2,j]]^3)*Y[b,ind[1,j]]}
				else {V1B[j,k]=-4*Y[b,k]*(1-Y[b,k])*Y[b,ind[1,j]]*Y[b,ind[2,j]]+8*(Y[b,k]^2)*Y[b,ind[1,j]]*Y[b,ind[2,j]]}

			}

		}


		V1C=matrix(0,sp,sp)
		for (j in 1:sp)
		{

			V1C[,j]=-Y[b,j]*(1-Y[b,j])
		}



		V2A=matrix(0,sp,sp)
		for (j in 1:sp)
		{

			for (i in 1:sp)
			{
				V2A[i,j]=-2*(Y[b,i]^2)*Y[b,j]
				V2A[j,i]=-2*(Y[b,j]^2)*Y[b,i]
			}
			V2A[j,j]=2*(Y[b,j]^2)*(1-Y[b,j])
		}


		V2B=matrix(0,qind,sp)
		for (j in 1:qind)
		{

			for (k in 1:sp)
			{
				if (ind[1,j]==k){V2B[j,k]=2*Y[b,ind[1,j]]*Y[b,ind[2,j]]*(1-2*Y[b,ind[1,j]])}
				else if (ind[2,j]==k){V2B[j,k]=2*Y[b,ind[2,j]]*Y[b,ind[1,j]]*(1-2*Y[b,ind[2,j]])}
				else {V2B[j,k]=-4*Y[b,k]*Y[b,ind[1,j]]*Y[b,ind[2,j]]}

			}

		}


		V2C=matrix(0,sp,sp)
		for (j in 1:sp)
		{

			for (i in 1:sp)
			{
				V2C[i,j]=-1*Y[b,j]
				V2C[j,i]=-1*Y[b,i]
			}
			V2C[j,j]=(1-Y[b,j])
		}


		V1=rbind(V1A,V1B,V1C)
		V2=rbind(V2A,V2B,V2C)

		V1S=matrix(0,sum(sp,sp,qind),1)
		V2S=matrix(0,sum(sp,sp,qind),1)
		V3S=matrix(0,sum(sp,sp,qind),sum(sp,sp,qind))

		for (j in 1:sp)
		{
			V1S=V1S+V1[,j]
			V2S=V2S+V2[,j]*Y[b,j]
			V3S=V3S+t(t(V2[,j]))%*%t(V2[,j])
		}

		V1S_sum=V1S_sum+V1S*w[b]
		V2S_sum=V2S_sum+V2S*w[b]
		V3S_sum=V3S_sum+V3S*w[b]


	}

	V1S_sum=V1S_sum/n
	V2S_sum=V2S_sum/n
	V3S_sum=V3S_sum/n

	alpha_p=betap+1

	d=((alpha_p*V2S_sum)-V1S_sum)
	W=V3S_sum
	ppi=solve(W)%*%d

	tot=sum(sp,qind,sp)
	ppi[sum(sp,qind,1):tot]=ppi[sum(sp,qind,1):tot]-1

	#ppi is the score matching estimate
        #for consistence add bL and betap
  ppiAll <- c(ppi[seq(1, p-1 + (p-2)*(p-1)/2)], #ALs
              rep(0,sp), #bL
              ppi[p-1 + (p-2)*(p-1)/2 + seq(1, sp)], #betaL
              betap)
	return(list(est = c(list(paramvec = ppiAll), ppi_parammats(ppiAll)),
                    SE = "Not calculated.",
                    info = list(W=W,d=d)))


}

usertheta_ppi_alr_gengamma_compatible <- function(usertheta){
  if (is.na(utils::tail(ppi_parammats(usertheta)$beta, 1))){return(FALSE)}
  p <- ppiltheta2p(length(usertheta))
  d_utheta <- ppi_paramvec(p, bL = 0, betap = utils::tail(ppi_parammats(usertheta)$beta, 1))
  if (isTRUE(all((d_utheta == usertheta)[!is.na(d_utheta)])) &&
      all(is.na(usertheta[is.na(d_utheta)])) ){return(TRUE)}
  else (return(FALSE))
}
