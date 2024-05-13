ppi_alrsqrt <- function(Y, paramvec = ppi_paramvec(p = ncol(Y), bL = 0, betap = 0), acut, bdryw = "minsq", w = rep(1, nrow(Y))){
  est1 <- ppi(Y, paramvec = paramvec,
            trans = "alr",
            w = w)
  paramspec <- ppi_parammats(paramvec)
  paramspec$beta <- est1$est$beta
  newparamvec <- do.call(ppi_paramvec, paramspec)
  est2 <- ppi(Y, paramvec = newparamvec,
              trans = "sqrt",
              bdryw = bdryw,
              acut = acut,
              w = w)
  if (is.list(est1$SE) && is.list(est2$SE)){ #as in SEs were estimated
    est2$SE$beta <- est1$SE$beta
    est2$SE$paramvec[-(ncol(Y)-1):0 + length(est2$SE$paramvec)] <- est1$SE$beta
  }
  return(est2)
}

ppi_alrsqrt_robust <- function(Y, cW, ...){
  ldenfun <- function(Y, theta){ #here theta is the usual parameters of PPI model from
    return(drop(dppi(Y, paramvec = theta)))
  }
  
  est <- Windham(Y = Y,
                    estimator = ppi_alrsqrt,
                    ldenfun = ldenfun,
                    cW = cW,
                    ...)
  out <- list(
    est = c(list(paramvec = est$paramvec), ppi_parammats(est$paramvec)),
    SE = "Not calculated.",
    info = est$optim
  )

  return(out)
}
