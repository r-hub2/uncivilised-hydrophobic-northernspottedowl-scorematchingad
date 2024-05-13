#' @title A PPI Score-Matching Marginal Moment Matching Estimator (dimension=3 only)
#' @description Computes a marginal moment matching estimator \insertCite{@Section 6.2, @scealy2023sc}{scorematchingad}, which assumes \eqn{\beta} is a known vector with the same value in each element, and \eqn{b_L = 0}.
#' Only \eqn{A_L} is estimated. 
#' @details
#' \eqn{\beta=\beta_0} is fixed and not estimated. \eqn{b_L} is fixed at zero.
#' See \insertCite{@Section 6.2 and A.8 of @scealy2023sc}{scorematchingad}.
#' The boundary weight function in the score matching discrepancy is the unthresholded product weight function
#' \deqn{h(z)^2 = \min\left(\prod_{j=1}^{p} z_j^2, a_c^2\right).}{h(z)^2 = min(z1^2 * z2^2 * ... * zp^2).}
#' @param Y Count data, each row is a multivariate observation.
#' @param ni The total for each sample (sum across rows)
#' @param beta0 \eqn{\beta=\beta_0}{beta=beta0} is the same for each component.
#' @param w Weights for each observation. Useful for weighted estimation in [`Windham()`].
#' @return A vector of estimates for \eqn{A_L} entries (diagonal and off diagonal).
#' @references
#' \insertAllCited{}
#' @export
ppi_mmmm <- function(Y, ni, beta0, w = rep(1, nrow(Y)))
{

	X1=Y
	X2=Y*(Y-1)
	X3=Y*(Y-1)*(Y-2)
	X4=Y*(Y-1)*(Y-2)*(Y-3)
	X5=Y*(Y-1)*(Y-2)*(Y-3)*(Y-4)
	X6=Y*(Y-1)*(Y-2)*(Y-3)*(Y-4)*(Y-5)

	fac1=ni
	fac2=ni*(ni-1)
	fac3=ni*(ni-1)*(ni-2)
	fac4=ni*(ni-1)*(ni-2)*(ni-3)
	fac5=ni*(ni-1)*(ni-2)*(ni-3)*(ni-4)
	fac6=ni*(ni-1)*(ni-2)*(ni-3)*(ni-4)*(ni-5)
	fac7=ni*(ni-1)*(ni-2)*(ni-3)*(ni-4)*(ni-5)*(ni-6)

	p3p1=weighted.mean(X3[,1]*X1[,2]/fac4, w=w)
	p4p1=weighted.mean(X4[,1]*X1[,2]/fac5, w=w)
	p3p2=weighted.mean(X3[,1]*X2[,2]/fac5, w=w)
	p2p1=weighted.mean(X2[,1]*X1[,2]/fac3, w=w)
	p2p2=weighted.mean(X2[,1]*X2[,2]/fac4, w=w)
	p1p3=weighted.mean(X1[,1]*X3[,2]/fac4, w=w)
	p1p4=weighted.mean(X1[,1]*X4[,2]/fac5, w=w)
	p2p3=weighted.mean(X2[,1]*X3[,2]/fac5, w=w)
	p1p2=weighted.mean(X1[,1]*X2[,2]/fac3, w=w)
	p2p2=weighted.mean(X2[,1]*X2[,2]/fac4, w=w)
	p3p3=weighted.mean(X3[,1]*X3[,2]/fac6, w=w)
	p1p1=weighted.mean(X1[,1]*X1[,2]/fac2, w=w)
	p5p1=weighted.mean(X5[,1]*X1[,2]/fac6, w=w)
	p4p2=weighted.mean(X4[,1]*X2[,2]/fac6, w=w)
	p6p1=weighted.mean(X6[,1]*X1[,2]/fac7, w=w)
	p5p2=weighted.mean(X5[,1]*X2[,2]/fac7, w=w)
	p3p4=weighted.mean(X3[,1]*X4[,2]/fac7, w=w)
	p4p3=weighted.mean(X4[,1]*X3[,2]/fac7, w=w)
	p2p4=weighted.mean(X2[,1]*X4[,2]/fac6, w=w)
	p1p5=weighted.mean(X1[,1]*X5[,2]/fac6, w=w)
	p2p5=weighted.mean(X2[,1]*X5[,2]/fac7, w=w)
	p1p6=weighted.mean(X1[,1]*X6[,2]/fac7, w=w)



	d=matrix(0,5,1)
	d[1]=32*p3p1-20*p4p1-20*p3p2-12*p2p1+12*p2p2
	d[2]=32*p1p3-20*p1p4-20*p2p3-12*p1p2+12*p2p2
	d[3]=48*p2p2-40*p3p2-40*p2p3-4*p1p2+4*p1p3-4*p2p1+4*p3p1
	d[4]=8*p2p1-6*p3p1-2*p1p1+2*p1p2-6*p2p2
	d[5]=8*p1p2-6*p2p2-6*p1p3-2*p1p1+2*p2p1



	dv=matrix(0,5,1)
	dv[1]=2*(4*p2p1-16*p3p1-4*p2p2+12*p4p1+12*p3p2)
	dv[2]=2*(4*p1p2-4*p2p2-16*p1p3+12*p2p3+12*p1p4)
	dv[3]=2*(4*p1p2+4*p2p1-32*p2p2-4*p1p3-4*p3p1+24*p3p2+24*p2p3)
	dv[4]=2*(2*p1p1-8*p2p1-2*p1p2+6*p3p1+6*p2p2)
	dv[5]=2*(2*p1p1-8*p1p2-2*p2p1+6*p1p3+6*p2p2)


	ev=matrix(0,5,1)
	ev[1]=12*(p3p1-p4p1-p3p2)-4*(p2p1-p3p1-p2p2)
	ev[2]=12*(p1p3-p2p3-p1p4)-4*(p1p2-p2p2-p1p3)
	ev[3]=24*(p2p2-p3p2-p2p3)-4*(p2p1-p3p1-p2p2)-4*(p1p2-p2p2-p1p3)
	ev[4]=6*(p2p1-p3p1-p2p2)-2*(p1p1-p2p1-p1p2)
	ev[5]=6*(p1p2-p2p2-p1p3)-2*(p1p1-p2p1-p1p2)



	d=d-dv+ev*(1+2*beta0[1])



	W=matrix(0,5,5)
	W[1,1]=16*p4p1-32*p5p1-16*p4p2+16*p6p1+16*p5p2
	W[1,2]=-16*p3p3+16*p3p4+16*p4p3
	W[1,3]=16*p3p2-48*p4p2-16*p3p3+32*p5p2+32*p4p3
	W[1,4]=8*p3p1-16*p4p1-8*p3p2+8*p5p1+8*p4p2
	W[1,5]=-8*p3p2+8*p4p2+8*p3p3
	W[2,1]=W[1,2]
	W[3,1]=W[1,3]
	W[4,1]=W[1,4]
	W[5,1]=W[1,5]
	W[2,2]=16*p1p4-16*p2p4-32*p1p5+16*p2p5+16*p1p6
	W[2,3]=16*p2p3-16*p3p3-48*p2p4+32*p3p4+32*p2p5
	W[2,4]=-8*p2p3+8*p3p3+8*p2p4
	W[2,5]=8*p1p3-8*p2p3-16*p1p4+8*p2p4+8*p1p5
	W[3,2]=W[2,3]
	W[4,2]=W[2,4]
	W[5,2]=W[2,5]
	W[3,3]=16*p2p3-96*p3p3-16*p2p4-16*p4p2+16*p3p2+64*p4p3+64*p3p4
	W[3,4]=8*p2p2-24*p3p2-8*p2p3+16*p4p2+16*p3p3
	W[3,5]=8*p2p2-8*p3p2-24*p2p3+16*p3p3+16*p2p4
	W[4,3]=W[3,4]
	W[5,3]=W[3,5]
	W[4,4]=4*p2p1-8*p3p1-4*p2p2+4*p4p1+4*p3p2
	W[4,5]=-4*p2p2+4*p3p2+4*p2p3
	W[5,4]=W[4,5]
	W[5,5]=4*p1p2-4*p2p2-8*p1p3+4*p2p3+4*p1p4

	pp=solve(W[1:3,1:3])%*%t(t(d[1:3]))

	mult=pp

	return(mult)

}


