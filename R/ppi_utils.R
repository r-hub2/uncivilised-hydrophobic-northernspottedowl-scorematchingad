#' @noRd
#' @title To or From vector form of parameters for PPI
#' purely for testing hardcoded ppi estimators - returns canonical parameters
toPPIcannparam <- function(ALs, bL, beta, manifold = "sphere"){
  if (manifold == "sphere"){out <- ppi_paramvec(AL = ALs, bL = bL, beta = 1 + 2 * beta)}
  else if (manifold == "simplex"){out <- ppi_paramvec(AL = ALs, bL = bL, beta = beta)}
  else {stop("manifold not supported")}
  return(out)
}

#' @noRd
#' @title returns the length of the parameter vector for given dimension p
ppithetalength <- function(p){
  p + #the beta
  (p-1) + #the diagonal of AL
  (p-2) * (p-1)/2 + #the upper triangle of AL
  (p-1) #the bL
}

#' @noRd
#' @title From the length of a ppi parameter vector, get the number of components
ppiltheta2p <- function(ltheta){#ltheta is length of theta
  #ltheta = p + (p-1) + (p-2) * (p-1)/2 + (p-1)
  # = 3p - 2 + (p^2 - 3p + 2)/2
  # 0 = p^2/2 + 1.5p -1-ltheta
  # 0 = p^2 + 3p - (2+2ltheta)
  #p = (-3 pm sqrt(9 + 4(2+2ltheta)) / 2
  #p = (-3 + sqrt(9 + 8 + 8ltheta)) /2
  p <- (-3 + sqrt(9 + 8 + 8 * ltheta))/2
  return(p)
}

#' @rdname ppi_param_tools
#' @order 3
#' @title Convert a PPI Parameter Vector to AL, bL and beta
#' @param paramvec A PPI parameter vector, typically created by [`ppi_paramvec()`] or as an output of [`ppi()`].
#' @return `ppi_parammats()`: A named list of \eqn{A_L}, \eqn{b_L}, and \eqn{\beta}.
#' @examples
#' vec <- ppi_paramvec(AL = rsymmetricmatrix(2), beta = c(-0.8, -0.7, 0))
#' ppi_parammats(vec)
#' @export
ppi_parammats <- function(paramvec){
  calcp <- ppiltheta2p(length(paramvec))
  p <- calcp
  AL <- tosmatrix(paramvec[1:((p-1) + (p-1)*(p-2)/2)])
  bL <- paramvec[p - 1 + ((p-2) * (p-1)/2) + 1:(p-1)]
  beta <- paramvec[(p - 1 + ((p-2) * (p-1)/2) + (p-1) + 1):length(paramvec)]
  return(list(
    AL = AL,
    bL = bL,
    beta = beta
  ))
}

#' @noRd
#' @title Functions to prepare microbiome data for unit tests
#' @description 
#' Prepares the micobiome data with or without two outliers.
#' Two different subsets of the measurement components are available.
#' @description Cleaned. TM7, Cyanobacteria/Chloroplast, Actinobacteria, Proteobacteria, other

ppi_microbiomedata_cleaned_TCAP <- function(){
  microdata <- scorematchingad::microbiome[(scorematchingad::microbiome$Year == 2008) & scorematchingad::microbiome$Helminth, ]
  microdata <- microdata[!microdata$IndividualID %in% c(2079, 2280), ] #remove two outlying measurements
  countdata=as.matrix(microdata[,12:31])

  #sample size
  n=92

  #dimension
  p=20

  #calculate totals
  tot=matrix(0,n,1)
  for (j in 1:p)
  {
   tot=tot+countdata[,j]
  }
  tot=as.vector(tot)

  #proportion data
  prop=countdata
  for (j in 1:n)
  {
    	prop[j,]=countdata[j,]/tot[j]
  }

  ####Reduce dimensions to p=5####


  #calculate 5D dataset
  comb=matrix(0,n,5)
  comb[,1]=prop[,"TM7"]
  comb[,2]=prop[,"Cyanobacteria/Chloroplast"]
  comb[,3]=prop[,"Actinobacteria"]
  comb[,4]=prop[,"Proteobacteria"]
  comb[,5]=abs(1-comb[,1]-comb[,2]-comb[,4]-comb[,3])
  propreal=comb
  colnames(propreal) <- c("TM7", "Cyanobacteria/Chloroplast", "Actinobacteria", "Proteobacteria", "pool")

  #dimension
  p=5

  #set beta (this is fixed here)
  beta0=matrix(0,p,1)
  beta0[1]=-0.8
  beta0[2]=-0.85
  beta0[3]=0
  beta0[4]=-0.2
  beta0[5]=0
  return(list(
    propreal = propreal,
    beta0 = beta0,
    p = p
  ))
}

#' @noRd
#' @description Not cleaned. TM7, Cyanobacteria/Chloroplast, Actinobacteria, Proteobacteria, other
ppi_microbiomedata_TCAP <- function(){
  microdata <- scorematchingad::microbiome[(scorematchingad::microbiome$Year == 2008) & scorematchingad::microbiome$Helminth, ]
  countdata=as.matrix(microdata[,12:31])
  
  #sample size
  n=94
  
  #dimension
  p=20
  
  #calculate totals
  tot=matrix(0,n,1)
  for (j in 1:p)
  {
   tot=tot+countdata[,j]
  }
  tot=as.vector(tot)
  
  #proportion data
  prop=countdata
  for (j in 1:n)
  {
  	prop[j,]=countdata[j,]/tot[j]
  }
  
  ###Reduce dimensions to p=5
  
  
  #calculate 5D dataset
  comb=matrix(0,n,5)
  comb[,1]=prop[,"TM7"]
  comb[,2]=prop[,"Cyanobacteria/Chloroplast"]
  comb[,3]=prop[,"Actinobacteria"]
  comb[,4]=prop[,"Proteobacteria"]
  comb[,5]=abs(1-comb[,1]-comb[,2]-comb[,4]-comb[,3])
  propreal=comb
  colnames(propreal) <- c("TM7", "Cyanobacteria/Chloroplast", "Actinobacteria", "Proteobacteria", "pool")
  
  
  #dimension
  p=5
  
  return(list(
    propreal = propreal,
    p = p
  ))
}


#' @noRd
#' @description Not cleaned. Spirochates, Verrucomicrobia, Cyanobacteria/Chloroplast, TM7 and pooled
ppi_microbiomedata_SVCTP <- function(){
  microdata <- scorematchingad::microbiome[(scorematchingad::microbiome$Year == 2008) & scorematchingad::microbiome$Helminth, ]
  countdata=as.matrix(microdata[,12:31])

  #sample size
  n=94

  #dimension
  p=20

  #calculate totals
  tot=matrix(0,n,1)
  for (j in 1:p)
  {
    tot=tot+countdata[,j]
  }
  tot=as.vector(tot)

  #proportion data
  prop=countdata
  for (j in 1:n)
  {
    prop[j,]=countdata[j,]/tot[j]
  }

  ##Reduce dimensions to p=5


  #dimension
  p=5

  #calculate 5D dataset
  comb=matrix(0,n,p)
  comb[,1]=prop[,"Spirochaetes"]
  comb[,2]=prop[,"Verrucomicrobia"]
  comb[,3]=prop[,"Cyanobacteria/Chloroplast"]
  comb[,4]=prop[,"TM7"]
  for (j in 1:sum(p,-1))
  {
    comb[,p]=comb[,p]+comb[,j]
  }
  comb[,p]=1-comb[,p]
  colnames(comb) <- c("Spirochaetes", "Verrucomicrobia", "Cyanobacteria/Chloroplast", "TM7", "pool")

  #save data
  propreal=comb

  return(list(
    propreal = propreal,
    p = ncol(propreal)
  ))
}
