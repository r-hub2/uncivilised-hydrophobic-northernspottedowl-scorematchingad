#' @name buildsmdtape
# Generates a `CppAD` tape of the empirical score matching discrepancy function \eqn{A(z) + B(z) + C(z)} ([`scorematchingtheory`]). Requires a tape of an improper log-likelihood function and the corresponding manifold with the transformation. 
# @param tranobj A transform object (of type `Rcpp_transform_ad`), typically created by [`manifoldtransform()`].
# @param man A transform object (of type `Rcpp_man_ad`) matching `tranobj`, and typically created by [`manifoldtransform()`].
# @param lltape Tape of the log-likelihood function constructed using [`tapell()`].
#' @param bdryw The name of the boundary weight function. "ones" for manifolds without boundary. For the simplex and positive orthant of the sphere, "prodsq" and "minsq" are possible - see [`ppi()`] for more information on these.
#' @param acut A parameter passed to the boundary weight function `bdryw`. Ignored for `bdryw = "ones"`.
#' @param verbose If `TRUE` more details are printed when taping. These details are for debugging and will likely be comprehensible only to users familiar with the source code of this package.
NULL

tapesmd <- function(lltape,
                   tranobj,
                   man,
                   bdryw,
                   acut = 1,
                   verbose = FALSE){
  inherits(tranobj, "Rcpp_transform_ad")
  inherits(man, "Rcpp_man_ad")
  stopifnot(is.numeric(acut))
  
  smdtape <- ptapesmd(attr(lltape, "ytape"),
                      lltape$dyntape,
                      lltape$ptr,
                      tranobj,
                      man,
                      bdryw, 
                      acut, 
                      verbose = verbose)
  out <- ADFun$new(ptr = smdtape,
                   name = paste(lltape$name, "smd", sep = "-"),
                   xtape = as.numeric(lltape$dyntape),
                   dyntape = as.numeric(attr(lltape, "ytape")),
                   usertheta = as.numeric(lltape$usertheta))
  return(out)
}
