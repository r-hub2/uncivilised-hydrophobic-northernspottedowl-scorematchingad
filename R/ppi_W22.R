
calcW22 <- function(p, h, h4m, w = rep(1, length(h))){

	DD=matrix(0,p,p)
	for (i in 1:p)
	{
		for (j in 1:p)
		{
			if (i==j){DD[i,j]=h4m[i]}
		}
	}

	DD2=matrix(weighted.mean(h^2, w = w),p,p)


	W=DD-DD2
  return(W)
}

# compute the first half of of the sample moment of wij(22) in section A.5.
# This is the sample moment of h(z)^2 * t(mu(22)_i) * mu(22)_j.
# When i=j this is non-zero, and valued at the sample moment of h(z)^2 / z_j^2
# The below calculates this for all j, for each measurement and averages.
# When h(z)^2 = 0, then the limit of h(z)^2/z_j^2 as z_j goes to zero is used,
# For the minima h(z) form, this limit is 1 when z_j is the minimum
# Note: indh is not used for the product style boundary weight function (hstyle)
h2onz2_mean <- function(p, n, z, h, indh, hstyle = "minima", w = rep(1, nrow(z))){
  if (hstyle=="minima"){
    limit <- matrix(0, n, p)
    for (i in 1:n){
      if (indh[i] > 0){
        limit[i, indh[[i]]] <- 1
      }
    }
  }
  if (hstyle=="product"){
    limit=matrix(1,n,p)
    for (k in 1:p)
    {
      for (j in 1:p)
      {
        if (k != j){limit[,k]=limit[,k]*z[,j]}
      }
    }
    limit = limit^2
  }

	h4s=matrix(0,n,p)
	h4m=matrix(0,1,p)
	for (j in 1:p)
	{
		for (i in 1:n)
		{
			if (h[i] > 0){h4s[i,j]=(h[i]^2)/z[i,j]^2}
		  else {h4s[i,j]=limit[i,j]}
		  # for Score1 estimators the indh gives the component that is smallest (or 0 if acut constraint hit)
      # so h4s[i,j] above is set to 1 whenever h is zero for that row and j is the first component that is zero.
		}
		h4m[j]=weighted.mean(h4s[,j], w = w)
	}
	return(h4m)
}
