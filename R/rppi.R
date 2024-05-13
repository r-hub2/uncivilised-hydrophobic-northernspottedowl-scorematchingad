#' @order 1
#' @title Simulate from a PPI Model
#' @family PPI model tools
#' @description Given parameters of the PPI model, generates independent samples.
#' @param n Number of samples to generate
#' @param paramvec The PPI parameter vector, created easily using [`ppi_paramvec()`] and also returned by [`ppi()`]. Use `paramvec` instead of `...`.
#' @param maxden This is the constant \eqn{log(C)} in \insertCite{@Appendix A.1.3 @scealy2023sc}{scorematchingad}.
#' @param maxmemorysize Advanced use. The maximum size, in bytes, for matrices containing simulated Dirchlet samples. The default of `1E5` corresponds to 100 mega bytes.
#' @return A matrix with `n` rows and `p` columns. Each row is an independent draw from the specified PPI distribution.
#' @inheritDotParams ppi_paramvec
#' @details
#' We recommend running `rppi()` a number of times to ensure the choice of `maxden` is good. `rppi()` will error when `maxden` is too low.
#' 
#' The simulation uses a rejection-sampling algorithm with Dirichlet proposal \insertCite{@Appendix A.1.3 @scealy2023sc}{scorematchingad}.
#' Initially `n` Dirichlet proposals are generated. After rejection there are fewer samples remaining, say \eqn{n^*}{n*}.
#' The ratio \eqn{n^*/n}{n*/n} is used to guess the number of new Dirichlet proposals to generate until `n` samples of the PPI model are reached. 
#'
#' Advanced use: The number of Dirichlet proposals created at a time is limited such that the matrices storing the Dirchlet proposals are always smaller than `maxmemorysize` bytes (give or take a few bytes for wrapping). 
#' Larger `maxmemorysize` leads to faster simulation so long as `maxmemorysize` bytes are reliably contiguously available in RAM.
#' @examples
#' beta0=c(-0.8, -0.8, -0.5)
#' AL = diag(nrow = 2)
#' bL = c(2, 3)
#' samp <- rppi(100,beta=beta0,AL=AL,bL=bL)
#' @export
rppi <- function(n, ..., paramvec = NULL, maxden = 4, maxmemorysize = 1E5){
  ellipsis::check_dots_used()
  #process inputs
  if (is.null(paramvec)){
    paramvec <- ppi_paramvec(...)
  }
  if (any(is.na(paramvec))){stop("All elements of paramvec must be non-NA. Did you forget to specify AL, bL or beta?")}
  parammats <- ppi_parammats(paramvec)
  beta <- parammats$beta
  AL <- parammats$AL
  bL <- parammats$bL
  p <- length(beta)
  
  # a warning if maxden is high
  if (maxden > 10){
    mess <- paste(sprintf("'maxden' of %0.2f is higher than 10.", maxden),
                                "When rppi() requires a high 'maxden' it could mean that",
                                "PPI density is hugely different from the Dirichlet component of the density.")
    if (any(beta < 0)){mess <- paste(mess, "This could mean that the concentrations on the boundary from the Dirichlet component will be too narrow to be represented in simulatad samples.")}
    warning(mess)
  }


  maxdenin <- maxden
  propaccepted <- 1 #start at 100 acceptance rate
  samples <- matrix(NA_real_, nrow = n, ncol = p) #create empty sample matrix
  nsamples <- 0
  nproposals <- 0
  maxden <- maxdenin
  maxblockrows <- floor(maxmemorysize/(p*8))
  stopifnot(maxblockrows > 1)
  while (nsamples < n){
    blocksize <- min(ceiling((n - nsamples)/propaccepted), maxblockrows) #choose so that the matrices of Dirichlet samples never get bigger than maxmemorysize bytes
    newsamples <- rppi_block(blocksize,p,beta,AL,bL,maxden)
    if (newsamples$maxden > maxden){warning(sprintf("Input maxden (i.e. the log(C) maximum) was %0.2f, but sampler suggests higher than %0.2f is required. You will have to call rppi() with a higher maxden.", maxdenin, newsamples$maxden))} 
    maxden <- newsamples$maxden
    nproposals <- nproposals + blocksize
    propaccepted <- (nsamples + nrow(newsamples$accepted))/nproposals
    if (nrow(newsamples$accepted)>0){
      newsampleskept <- newsamples$accepted[1:min(nrow(newsamples$accepted), n-nsamples), , drop = FALSE] #keep up to the n requested samples
      samples[nsamples+(1:nrow(newsampleskept)), ] <- newsampleskept
      nsamples <- nsamples + nrow(newsampleskept)
    }
    # continue until n or more samples accepted
  } 

  #maxden is the constant log(C) in Appendix A.1.3. Need to run the sampler
  #a few times to check that it is an appropriate upper bound.
  if (maxden > maxdenin){stop(sprintf("Input maxden (i.e. the log(C) maximum) was %0.2f, but sampler suggests higher than %0.2f is required.", maxdenin, maxden))}

  stopifnot(all(is.finite(samples)))

  return(samples)
}


# simulates samples one at a time
rppi_singly <- function(n,p,beta,AL,bL,maxden)
{

	alpha=beta+1
	coun=0
	samp1=matrix(0,1,p)
	count2=0
	while (coun < n)
	{

		count2=count2+1
		Uni=MCMCpack::rdirichlet(1, alpha)
		u=stats::runif(1,0,1)
		Uni_nop <- Uni[1:(p-1)]
		tUni_nop <- t(Uni[1:(p-1)])
    num <- ppi_uAstaru(matrix(Uni_nop, nrow = 1), AL, bL) - maxden #uT * AL * u + t(bL) * u - maxden
		if (num > 0){maxden=num+maxden}
		#print(maxden)
		if (u < exp(num)){samp1=rbind(samp1,Uni);coun=coun+1}
		#print(coun)
	}
	samp3=samp1[2:length(samp1[,1]),]

	return(list(samp3=samp3,maxden=maxden))
}



# try generating samples without a 'while' loop
rppi_block <- function(n,p,beta,AL,bL,maxden){
  Uni <- MCMCpack::rdirichlet(n, beta+1)
  Uni_nop <- Uni[, -p]
  nums <- ppi_uAstaru(Uni_nop, AL, bL) - maxden #uT * AL * u + t(bL) * u - maxden

  # update maxdens
  max_num <- max(nums)
  if (max_num > 0) {maxden=max_num+maxden}

  # accept some of points
  unif <- stats::runif(n,0,1)
  accepted <- Uni[unif < exp(nums), , drop = FALSE]

  return(list(accepted = accepted, maxden = maxden))
}

#' @title Improper Log-Density of the PPI Model
#' @family PPI model tools
#' @description Compute the __natural logarithm__ of the improper density for the PPI model for the given matrix of measurements `Y`. Rows with negative values or with a sum that differs from `1` by more than `1E-15` are assigned a value of `-Inf`.
#' @param Y A matrix of measurements in the simplex. Each row is a multivariate measurement.
#' @inheritParams rppi
#' @param ... Arguments passed on to [`ppi_paramvec()`].
#' @inheritDotParams ppi_paramvec
#' @details The value calculated by `dppi` is
#' \deqn{z_L^TA_Lz_L + b_L^Tz_L + \beta^T \log(z),}
#' where \eqn{z} is the multivariate observation (i.e. a row of `Y`), and \eqn{z_L} ommits the final element of \eqn{z}.
#' @examples
#' m <- rppi_egmodel(10)
#' dppi(m$sample, paramvec = m$theta)
#' @export
dppi <- function(Y, ..., paramvec = NULL){
  ellipsis::check_dots_used()
  #process inputs
  if (is.null(paramvec)){
    paramvec <- ppi_paramvec(...)
  }
  if (any(is.na(paramvec))){stop("All elements of paramvec must be non-NA. Did you forget to specify AL, bL or beta?")}
  parammats <- ppi_parammats(paramvec)
  beta <- parammats$beta
  AL <- parammats$AL
  bL <- parammats$bL
  
  p <- ncol(Y)
  if (!("matrix" %in% class(bL))){bL <- as.matrix(bL, ncol = 1)}
  stopifnot(isTRUE(ncol(bL) == 1))

  uAstaru <- ppi_uAstaru(Y[,-p, drop = FALSE], AL, bL) #result is a vector
  if (all(beta == 0)){return(as.vector(uAstaru))} #skip the computation below when beta0 is zero

  if (!("matrix" %in% class(beta))){beta <- as.matrix(beta, ncol = 1)}
  logprop <- log(Y)
  # define u^0 as 1 when u goes to -Inf
  logprop[, beta == 0] <- 0
  logdirichlet <- logprop %*% beta
  logdensity <- as.vector(uAstaru + logdirichlet)

  # set points outside the simplex to 0
  negatives <- rowSums(Y < 0) > 0
  sumisnt1 <- abs(rowSums(Y) -1) > 1E-15
  logdensity[negatives|sumisnt1] <- -Inf
  return(logdensity)
}

# below is function for uT * ALs * u + t(bL) * u
# prop_nop is the measurements WITHOUT the last column
ppi_uAstaru <- function(prop_nop, ALs, bL){
  stopifnot(ncol(prop_nop) == ncol(ALs))
  stopifnot(ncol(ALs) == nrow(bL))
  stopifnot(ncol(ALs) == nrow(ALs))
  nums <- rowSums((prop_nop %*% ALs) * prop_nop) + #uT * ALs * u
    as.vector(prop_nop %*% bL) #+ u * bL
  return(nums)
}
