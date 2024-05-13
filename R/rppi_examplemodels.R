# @title Simulate from an example PPI model
# @family PPI model tools
#' @return `rppi_egmodel` returns a list: 
#'  * `sample` A matrix of the simulated samples (`n` rows)
#'  * `p` The number of components of the model
#'  * `theta` The PPI parameter vector
#'  * `AL` The \eqn{A_L} parameter matrix
#'  * `bL` The \eqn{b_L} parameter vector
#'  * `beta` The \eqn{\beta} parameter vector 
#' @examples
#' rppi_egmodel(1000)
#' @describeIn rppi 
#' Simulates the 3-component PPI model from \insertCite{@Section 2.3, @scealy2023sc}{scorematchingad} and returns both simulations and model parameters. 
#' @order 2
#' @references
#' \insertAllCited{}
#' @export
rppi_egmodel <- function(n, maxden = 4){
  mats <- pars_sec2dot3model(3)
  mats$beta = mats$beta

  #simulate sample from PPI model
  samp3=rppi(n,beta = mats$beta,AL = mats$AL,bL = mats$bL, maxden = maxden)

  out <- c(list(
    sample = samp3,
    p = 3,
    theta = ppi_paramvec(AL = mats$AL, bL = mats$bL, beta = mats$beta)
    ),
    mats)
  return(out)
}

# @rdname rppi_egmodel
# @export
rppi_egmodel_p4 <- function(n, maxden = 8){
  mats <- pars_sec2dot3model(4)

  #simulate sample from PPI model
  samp3=rppi(n,beta = mats$beta,AL = mats$AL,bL = mats$bL, maxden = maxden)

  out <- c(list(
    sample = samp3,
    p = 4,
    theta = ppi_paramvec(AL = mats$AL, bL = mats$bL, beta = mats$beta)
    ),
    mats)
  return(out)
}


# internal function for building Section 2.3-like model for any number of components
# it pretty poor for anything but p = 3
# returns the parameter matrices/vectors
pars_sec2dot3model <- function(p){
  #parameters for the PPI model
  muL=matrix(0,p-1,1)
  muL[1:sum(p,-1)]=0.12
  aa=matrix(1/500,p-1,1)
  D=diag(as.vector(aa))
  SigA=D
  SigA[1,1]=SigA[1,1]*2
  cor=0.5
  SigA[1,2]=cor*sqrt(SigA[1,1]*SigA[2,2])
  SigA[2,1]=SigA[1,2]
  ALs=-0.5*solve(SigA)
  bL=solve(SigA)%*%muL
  beta0=matrix(-0.8,p,1)
  beta0[p]=-0.5
  return(list(
    AL = ALs,
    bL = bL,
    beta = beta0
  ))
}
