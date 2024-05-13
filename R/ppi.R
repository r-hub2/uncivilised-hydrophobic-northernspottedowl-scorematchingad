#' @title Estimation of Polynomially-Tilted Pairwise Interaction (PPI) Model
#' @order 1
#' @family PPI model tools
#' @description 
#' Estimates the parameters of the Polynomially-Tilted Pairwise Interaction (PPI) model \insertCite{scealy2023sc}{scorematchingad} for compositional data.
#' By default `ppi()` uses [`cppad_closed()`] to find estimate.
#' For many situations a hard-coded implementation of the score matching estimator is also available.

#' @section PPI Model:
#' The PPI model density is proportional to
#' \deqn{\exp(z_L^TA_Lz_L + b_L^Tz_L)\prod_{i=1}^p z_i^{\beta_i},}
#' where \eqn{p} is the dimension of a compositional measurement \eqn{z}, and \eqn{z_L} is \eqn{z} without the final (\eqn{p}th) component.
#' \eqn{A_L} is a \eqn{p-1 \times p-1} symmetric matrix that controls the covariance between components.
#' \eqn{b_L} is a \eqn{p-1} vector that controls the location within the simplex.
#' The \eqn{i}th component \eqn{\beta_i} of \eqn{\beta} controls the concentration of density when \eqn{z_i} is close to zero: when \eqn{\beta_i \geq 0} there is no concentration and \eqn{\beta_i} is hard to identify; when \eqn{\beta_i < 0} then the probability density of the PPI model increases unboundedly as \eqn{z_i} approaches zero, with the increasing occuring more rapidly and sharply the closer \eqn{\beta_i} is to \eqn{-1}.
#' 
#' @details
#' Estimation may be performed via transformation of the measure in Hyvärinen divergence from Euclidean space to the simplex (inverse of the additive log ratio transform), from a hyperplane to the simplex (inverse of the centred log ratio transform), from the positive quadrant of the sphere to the simplex (inverse of the square root transform), or without any transformation. In the latter two situations there is a boundary and *weighted Hyvärinen divergence* \insertCite{@Equation 7, @scealy2023sc}{scorematchingad} is used.
#' Properties of the estimator using the square root transform were studied by \insertCite{scealy2023sc;textual}{scorematchingad}.
#' Properties of the estimator using the additive log ratio transfrom were studied by \insertCite{scealy2024ro;textual}{scorematchingad}.
#'
#' There are three boundary weight functions available:
#' * The function "ones" applies no weights and should be used whenever the manifold does not have a bounday.
#' * The function "minsq" is the minima-based boundary weight function for the PPI model \insertCite{@Equation 12, @scealy2023sc}{scorematchingad}
#' \deqn{\tilde{h}(z)^2 = \min(z_1^2, z_2^2, ..., z_p^2, a_c^2).}{h(z)^2 = min(z1^2, z2^2, ..., zp^2, a_c^2),}
#' where \eqn{z} is a point in the positive orthant of the p-dimensional unit sphere
#' and \eqn{z_j}{zj} is the jth component of z.
#' * The function "prodsq" is the product-based \insertCite{@Equation 9, @scealy2023sc}{scorematchingad}
#' \deqn{\tilde{h}(z)^2 = \min(\prod_{j=1}^{p} z_j^2, a_c^2).}{h(z)^2 = min(z1^2 * z2^2 * ... * zp^2, a_c^2),}
#' where \eqn{z} is a point in the positive orthant of the p-dimensional unit sphere
#' and \eqn{z_j}{zj} is the jth component of z.
#'
#' Scealy and Wood \insertCite{@Theorem 1, @scealy2023sc}{scorematchingad} prove that minimising the weighted Hyvärinen Divergence is equivalent to minimising \eqn{\psi(f, f_0)} (See [`scorematchingtheory`])
#' when the boundary weight function is smooth or for the functions "minsq" and "prodsq"  above when the manifold is the simplex or positive orthant of a sphere.
#'
#' Hard-coded estimators are available for the following situations:
#'  + Square root transformation ("sqrt") with the "minsq" boundary weight function:
#'    + full parameter vector (`paramvec` not provided)
#'    + `paramvec` fixes only the final element of \eqn{\beta}
#'    + `paramvec` fixes all elements of \eqn{\beta} 
#'    + `paramvec` fixes \eqn{b_L = 0} and provides fixed values of \eqn{\beta}
#'    + `paramvec` fixes \eqn{A_L=0} and \eqn{b_L=0}, leaving \eqn{\beta} to be fitted. 
#'  + Square root transformation ("sqrt") with the "prodsq" boundary weight function:
#'    + `paramvec` fixes all elements of \eqn{\beta} 
#'    + `paramvec` fixes \eqn{b_L = 0} and provides fixed values of \eqn{\beta}
#'    + `paramvec` fixes \eqn{A_L=0} and \eqn{b_L=0}, leaving \eqn{\beta} to be fitted. 
#'  + The additive log ratio transformation ("alr") using the final component on the denominator, with \eqn{b_L=0} and fixed final component of \eqn{\beta}.

#' @param trans The name of the transformation of the manifold in Hyvärinen divergence (See [`scorematchingtheory`]): "clr" (centred log ratio), "alr" (additive log ratio), "sqrt" or "none".
#' @param Y A matrix of measurements. Each row is a compositional measurement (i.e. each row sums to 1 and has non-negative elmenents).
#' @param paramvec Optionally a vector of the PPI models parameters. `NA`-valued elements of this vector are estimated and non-`NA` values are fixed. Generate `paramvec` easily using [ppi_paramvec()].  If `NULL` then all elements of \eqn{A_L}, \eqn{b_L} and \eqn{\beta} are estimated.
#' @param bdryw The boundary weight function for down weighting measurements as they approach the manifold boundary. Either "ones", "minsq" or "prodsq". See details.
#' @param acut The threshold \eqn{a_c} in `bdryw` to avoid over-weighting measurements interior to the simplex
#' @param control `iterative` only. Passed to [`optimx::Rcgmin()`] to control the iterative solver.
#' @param bdrythreshold `iterative` or `closed` methods only. For measurements within `bdrythreshold` of the simplex boundary a Taylor approximation is applied by shifting the measurement `shiftsize` towards the center of the simplex.
#' @param shiftsize `iterative` or `closed` methods only. For measurements within `bdrythreshold` of the simlex boundary a Taylor approximation is applied by shifting the measurement `shiftsize` towards the center of the simplex.
#' @param approxorder `iterative` or `closed` methods only. Order of the Taylor approximation for measurements on the boundary of the simplex.
#' @param method `"closed"` uses `CppAD` to solve in closed form the a quadratic score matching discrepancy using [`cppad_closed()`]. `"hardcoded"` uses hardcoded implementations. "iterative" uses [`cppad_search()`] (which uses `CppAD` and [`optimx::Rcgmin()`]) to iteratively find the minimum of the weighted Hyvärinen divergence.
#' @param paramvec_start `iterative` method only. The starting guess for `Rcgmin`. Generate `paramvec_start` easily using [`ppi_paramvec()`].
#' @param w Weights for each observation, if different observations have different importance. Used by [`Windham()`] and [`ppi_robust()`] for robust estimation.
#' @param constrainbeta If `TRUE`, elements of \eqn{\beta} that are less than `-1` are converted to `-1 + 1E-7`.
#' @return `ppi()` returns:
#' A list of `est`, `SE` and `info`.
#'  * `est` contains the estimates in vector form, `paramvec`, and as \eqn{A_L}, \eqn{b_L} and \eqn{\beta}.
#'  * `SE` contains estimates of the standard errors if computed. See [`cppad_closed()`].
#'  * `info` contains a variety of information about the model fitting procedure and results.

#' @references
#' \insertAllCited{}
#' @examples
#' model <- rppi_egmodel(100)
#' estalr <- ppi(model$sample,
#'               paramvec = ppi_paramvec(betap = -0.5, p = ncol(model$sample)),
#'               trans = "alr")
#' estsqrt <- ppi(model$sample,
#'               trans = "sqrt",
#'               bdryw = "minsq", acut = 0.1)
#' @export
ppi <- function(Y, paramvec = NULL,
                trans, method = "closed", w = rep(1, nrow(Y)),
                constrainbeta = FALSE,
                bdryw = "ones", acut = NULL, #specific to some methods
                bdrythreshold = 1E-10, shiftsize = bdrythreshold, approxorder = 10, control = list(tol = 1E-15, checkgrad = TRUE), paramvec_start = NULL#specific to cppad methods
                ){
  ###### process inputs #####
  stopifnot("matrix" %in% class(Y))
  p = ncol(Y)
  numneg <- sum(Y < 0)
  if (numneg > 0){
     warning(sprintf("Y contains %i negative values.", numneg))
  }
  sum_m1 <- max(abs(rowSums(Y) - 1))
  if (sum_m1 > 1E-15){
     warning(sprintf("Y contains measurement that don't add to 1. Largest discrepancy is %s.", sum_m1))
  }

  stopifnot(trans %in% c("alr", "sqrt", "clr", "none"))
  man <- switch(trans,
           alr = "Euc",
           clr = "Hn111",
           sqrt = "sph",
           none = "sim")
  if (trans == "none"){trans <- "identity"}

  if (man %in% c("sim", "sph")){
    if (bdryw == "ones"){stop("Manifold supplied has a boundary - set bdryw to something that isn't 'ones'")}
  } else {
    if (bdryw != "ones"){warning("Manifold supplied has no boundary. Setting bdryw to 'ones'.")}
  }
  if (bdryw == "ones"){
    if (!is.null(acut)){warning("The value of 'acut' is ignored for bdryw == 'ones'")}
    acut <- 1 #set just for passing to CppAD
  }
  if (is.null(paramvec)){usertheta <- rep(NA, ppithetalength(p))}
  else {usertheta <- paramvec}

  stopifnot(length(usertheta) == ppithetalength(p))

  fitobj <- list()
  fitfun <- NA

  ########## DO THE FITTING ###############

  
  # hardcoded methods:
  if (method == "hardcoded"){
    if (trans == "alr"){
      if (usertheta_ppi_alr_gengamma_compatible(usertheta)){
        fitobj <- ppi_alr_gengamma(Y, betap = utils::tail(usertheta, 1), w = w) #any theta is fine
        fitfun <- "ppi_alr_gengamma"
      }
    }
    if (trans == "sqrt"){ # a number of methods implemented
      if (bdryw == "minsq"){
        if (ppi_usertheta_for_dir_sqrt_minimah(usertheta)){
          betaest <- as.vector(dir_sqrt_minimah(Y, acut = acut, w = w))
          estparamvec <- t_fu2t(betaest, usertheta)
          fitobj$est <- c(list(paramvec = estparamvec),
                            ppi_parammats(estparamvec))
          fitobj$SE <- "Not calculated."
          fitfun <- "dir_sqrt_minimah"
        } else if (ppi_usertheta_estimator1_compatible_zerob(usertheta)){
          fitobj <- estimator1(Y,acut = acut,incb = 0,
                            beta = ppi_parammats(usertheta)$beta,
                            w= w, computeSE = TRUE)
          fitfun <- "estimator1_zerob"
        } else if (ppi_usertheta_estimator1_compatible_incb(usertheta)){
          fitobj <- estimator1(Y,acut = acut,incb = 1,
                            beta = ppi_parammats(usertheta)$beta,
                            w= w, computeSE = TRUE)
          fitfun <- "estimator1_incb"
        } else if (utheta_estimatorall1_betap_compatible(usertheta)){
          fitobj <- ppi_sqrt_minimah_full(Y, acut, utils::tail(ppi_parammats(usertheta)$beta, 1),
                                            w)
          fitfun <- "ppi_sqrt_minimah_betap"
        } else if (utheta_estimatorall1_full_compatible(usertheta)){
          fitobj <- ppi_sqrt_minimah_full(Y, acut, betap = NULL,
                                            w)
          fitfun <- "ppi_sqrt_minimah_full"
        }
      }

      if (bdryw == "prodsq"){
        if (ppi_usertheta_for_dir_sqrt_minimah(usertheta)){
          betaest <- dir_sqrt_prodh(Y, acut = acut, w = w)
          fitfun <- "dir_sqrt_prodh"
          estparamvec <- t_fu2t(betaest, usertheta)
          fitobj$est <- c(list(paramvec = estparamvec),
                            ppi_parammats(estparamvec))
          fitobj$SE <- "Not calculated."
        } else if (ppi_usertheta_estimator1_compatible_zerob(usertheta)){
          fitobj <- ppi_sqrt_prodh_zerob(Y, acut, beta = ppi_parammats(usertheta)$beta, w)
          fitfun <- "ppi_sqrt_prodh_zerob"
        } else if (ppi_usertheta_estimator1_compatible_incb(usertheta)){
          fitobj <- ppi_sqrt_prodh(Y, acut, beta = ppi_parammats(usertheta)$beta, w)
          fitfun <- "ppi_sqrt_prodh"
        }
      }
    }
    if (is.na(fitfun)){
      warning("No hard-coded estimator exists for parameter set. Using cppad closed-form solution.")
      method <- "closed"
    }
  }

  # cppad methods
  if (method %in% c("iterative", "closed")){
    if (is.null(paramvec_start)){stheta <- t_u2s_const(usertheta, 0.2)}
    else {stheta <- t_us2s(usertheta, paramvec_start)}
    isfixed <- t_u2i(usertheta)
    # prepare tapes
    tapes <- buildsmdtape(
       start = "sim",
       tran = trans,
       end = man,
       ll = "ppi",
       ytape =  rep(1/p, p),
       usertheta = t_si2u(stheta, isfixed),
       bdryw = bdryw,
       acut = acut,
       verbose = FALSE)
    smdtape <- tapes$smdtape
  
    # find boundary points and their approximation centres
    isbdry <- simplex_isboundary(Y, bdrythreshold)
    Yapproxcentres <- Y 
    Yapproxcentres[!isbdry, ] <- NA
    Yapproxcentres[isbdry, ] <- simplex_boundaryshift(Y[isbdry, , drop = FALSE], shiftsize = shiftsize)

    fitobj <- list()
    if (method == "iterative"){ 
      rawout <- cppad_search(smdtape = smdtape,
                   theta = t_si2f(stheta, isfixed),
                   Y = Y,
                   Yapproxcentres = Yapproxcentres,
                   w = w,
                   approxorder = approxorder,
                   control = control)
  
      fitobj$info$smval <- rawout$value
      fitfun <- "iterative"
    }
    if (method == "closed"){
      rawout <- cppad_closed(smdtape = smdtape,
                   Y = Y,
                   Yapproxcentres = Yapproxcentres,
                   w = w,
                   approxorder = approxorder)
      fitfun <- "closed"
    }
    #process the theta and SE
    thetaest <- t_sfi2u(rawout$est, stheta, isfixed)
                   
    #refactor results to fit with ppi() standard output
    fitobj$est <- c(list(paramvec = thetaest),
                      ppi_parammats(thetaest))
    if (isa(rawout$SE, "numeric")){
      SE <- t_sfi2u(rawout$SE, rep(0, length(stheta)), isfixed)
      fitobj$SE <- c(list(paramvec = SE),
                      ppi_parammats(SE))
    } else {fitobj$SE <- rawout$SE}
    fitobj$info <- rawout 
    fitobj$info$boundarypoints <- sum(isbdry)
  }

  if (length(fitobj) == 0){stop("Could not find fitting method.")}

  if (constrainbeta & (any(fitobj$est$beta <= -1))){
    warning("Beta estimates of -1 or lower, replacing these with -1 + 1E-7")
    fitobj$est$beta[fitobj$est$beta <= -1] <- -1 + 1E-7
    fitobj$est$paramvec <- do.call(ppi_paramvec, fitobj$est[c("AL", "bL", "beta")])
  }

  fitobj$info$method <- fitfun
  return(fitobj)
}



# Tests of the ppi() structure
expect_ppi_str <- function(object, p){
  act <- testthat::quasi_label(rlang::enquo(object), arg = "object")

  #check
  correctfirstnames <- isTRUE(setequal(names(object), c("est", "SE", "info")))
  if (!correctfirstnames){
    testthat::fail("Incorrect slot names")
  }

  correctest <- correcteststr(object$est, p)
  if ((length(object$SE) == 1) && (object$SE == "Not calculated.")){
    correctSE = list(SEnotcalculated = TRUE)
  } else {
    correctSE <- correcteststr(object$SE, p)
  }

  correctmethodstr <- isTRUE((class(object$info$method) == "character") & (nchar(object$info$method) > 0))
  passsummary <- list(
    slotnames = correctfirstnames,
    est_structure = correctest,
    SE_structure = correctSE,
    info_method_exists = correctmethodstr
  )
  if (all(unlist(passsummary))){
    testthat::succeed()
    return(invisible(act$val))
  }
  else {
    fails <- names(which(!unlist(passsummary)))
    fails <- gsub(".", "$", fails, fixed = TRUE)
    testthat::fail(
      sprintf("The following conditions fail: %s", paste(fails, collapse = ", ")))
  }
}

correcteststr <- function(est, p){
  correctestnames <- setequal(names(est), c("paramvec", "AL", "bL", "beta"))
  correctestdims <- list()
  if (correctestnames){
    correctestdims$paramvec <- (length(est$paramvec) == ppithetalength(p))
    correctestdims$AL <- isTRUE(all(dim(est$AL) == c(p-1, p-1)))
    correctestdims$bL <- (length(est$bL) == p-1)
    correctestdims$beta <- (length(est$beta) == p)
  }
  return(c(Names = correctestnames, correctestdims))
}


