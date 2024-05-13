#' @title Windham Robustification of Point Estimators for Exponential Family Distributions
#' @family generic score matching tools
#' @description Performs a generalisation of Windham's robustifying method \insertCite{windham1995ro}{scorematchingad} for exponential models with natural parameters that are a linear function of the parameters for estimation.
#' Estimators must solve estimating equations of the form
#' \deqn{\sum_{i = 1}^n U(z_i; \theta) = 0.}
#' The estimate is found iteratively through a fixed point method as suggested by \insertCite{windham1995ro;textual}{scorematchingad}.
#'

#' @param Y A matrix of measurements. Each row is a measurement, each component is a dimension of the measurement.
#' @param estimator A function that estimates parameters from weighted observations.
#' It must have arguments `Y` that is a matrix of measurements and `w` that are weights associated with each row of `Y`. If it accepts arguments `paramvec` or `paramvec_start` then these will be used to specify fixed elements of the parameter vector and the starting guess of the parameter vector, respectively. The estimated parameter vector, including any fixed elements, must be the returned object, or the first element of a returned list, or as the `paramvec` slot within the `est` slot of the returned object.
#' @param ldenfun A function that returns a vector of values propotional to the log-density for a matrix of observations `Y` and parameter vector `theta`.
#' @param ... Arguments passed to `estimator`.
#' @param paramvec_start
#' Initially used to check the function `estimator`. If `estimator` accepts a `paramvec_start`, then the current estimate of the parameter vector is passed as `paramvec_start` to `estimator` in each iteration.
#' @param cW A vector of robustness tuning constants. When computing the weight for an observation the parameter vector is multiplied element-wise with `cW`. For the PPI model, generate `cW` easily using [ppi_cW()] and [ppi_cW_auto()].
#' @param fpcontrol A named list of control arguments to pass to [FixedPoint::FixedPoint()] for the fixed point iteration.
#' @param alternative_populationinverse The default is to use [`Windham_populationinverse()`]. If TRUE an alternative implementation in [`Windham_populationinverse_alternative()`] is used. So far we have not seen any difference between the results.


#' @details
#' For any family of models with density \eqn{f(z; \theta)}, Windham's method finds the parameter set \eqn{\hat\theta} such that the estimator applied to observations weighted by \eqn{f(z; \hat\theta)^c} returns an estimate that matches the theoretical effect of weighting the full population of the model.
#' When \eqn{f} is proportional to \eqn{\exp(\eta(\theta) \cdot T(z))} and \eqn{\eta(\theta)} is linear, these weights are equivalent to \eqn{f(z; c\hat\theta)} and the theoretical effect of the weighting on the full population is to scale the parameter vector \eqn{\theta} by \eqn{1+c}.
#'
#' The function `Windham()` assumes that \eqn{f} is proportional to \eqn{\exp(\eta(\theta) \cdot T(z))} and \eqn{\eta(\theta)} is linear. It allows a generalisation where \eqn{c} is a vector so the weight for an observation \eqn{z} is \deqn{f(z; c \circ \theta),} where \eqn{\theta} is the parameter vector, \eqn{c} is a vector of tuning constants, and \eqn{\circ} is the element-wise product (Hadamard product).
#'
#' The solution is found iteratively \insertCite{windham1995ro}{scorematchingad}. 
#' Given a parameter set \eqn{\theta_n}, `Windham()` first computes weights \eqn{f(z; c \circ \theta_n)} for each observation \eqn{z}.
#' Then, a new parameter set \eqn{\tilde{\theta}_{n+1}} is estimated by `estimator` with the computed weights.
#' This new parameter set is element-wise-multiplied by the (element-wise) reciprical of \eqn{1+c} to obtain an adjusted parameter set \eqn{\theta_{n+1}}.
#' The estimate returned by `Windham()` is the parameter set \eqn{\hat{\theta}} such that \eqn{\theta_n \approx \theta_{n+1}}.
#' @family Windham robustness functions
#' @return
#' A list:
#' * `paramvec` the estimated parameter vector
#' * `optim` information about the fixed point iterations and opimisation process. Including a slot `finalweights` for the weights in the final iteration.
#' @examples
#' if (requireNamespace("movMF")){
#'   Y <- movMF::rmovMF(1000, 100 * c(1, 1) / sqrt(2))
#' } else {
#'   Y <- matrix(rnorm(1000 * 2, sd = 0.01), ncol = 2)
#'   Y <- Y / sqrt(rowSums(Y^2))
#' }
#' Windham(Y = Y,
#'    estimator = vMF,
#'    ldenfun = function(Y, theta){ #here theta is km
#'      return(drop(Y %*% theta))
#'    },
#'    cW = c(0.01, 0.01),
#'    method = "Mardia")
#' @export
Windham <- function(Y, estimator, ldenfun, cW, ..., fpcontrol = list(Method = "Simple", ConvergenceMetricThreshold = 1E-10), paramvec_start = NULL, alternative_populationinverse = FALSE){#... earlier so that fpcontrol and paramvec_start can only be passed by being named
  extraargs <- list(...)
  ellipsis::check_dots_used()
  # assuming estimator has arguments: Y, paramvec, w, and optionally paramvec_start.
  # and assume that the return vector can be extracted using `extract_paramvec()` and similar and that this is the full model parameter vector, including the fixed elements (this is important for computing density).
  estargs <- c(list(Y = Y), extraargs)
  estargs$paramvec_start <- paramvec_start #adding this slot this way so that it is omitted if NULL

  #assess the passes estimator
  assessment <- do.call(Windham_assess_estimator, c(list(estimator = estimator), estargs[intersect(names(formals(estimator)), names(estargs))]))


  # extract start vector from a paramvec and paramvec_start
  if (!is.null(paramvec_start)){
    if (!is.null(extraargs$paramvec)){
      starttheta <- t_us2s(extraargs$paramvec, paramvec_start)
    } else {
      starttheta <- paramvec_start
    }
  } else { #use estimator for the start values
    starttheta <- assessment$est
  }

  # calculate isfixed
  if (!is.null(extraargs$paramvec)){isfixed = t_u2i(extraargs$paramvec)}
  else {isfixed <- rep(FALSE, length(starttheta))}

  # cW checks
  stopifnot(length(cW) == length(starttheta))
  stopifnot(is.numeric(cW))
  if (any((cW * starttheta)[isfixed] != 0)){stop("Elements of cW corresponding to fixed non-zero parameters should be zero")}

  # Correction of parameter preparation
  # use the WindhamCorrection(), the alternative is Scealy's original additive method in the draft paper
  if (alternative_populationinverse){
   if (length(cW) > 1){ if (stats::var(cW[cW > 1E-10]) > (1E-10)^2){ #require constant cW (or zero) because I'm not sure what Scealy's correction method should be in the presence of a different tuning constants per value
     stop("Non-zero cW values vary, which is not supported by 'additive' correction of the parameter estimate.")
   }}
   inWW <- (cW > 1E-10)
   cWav <- mean(cW[cW > 1E-10]) #note that cW ~~ inWW * cWav
   thetaadjuster <- Windham_populationinverse_alternative
  } else {
    tauinv <- Windham_populationinverse(cW)
    cWav <- NULL  #not relevant to this correction method
    thetaadjuster <- function(newtheta, previoustheta = NULL, cW = NULL, cWav = NULL){tauinv %*% newtheta}
  }

  # functions for adding paramvec_start or otherwise to estimator arguments
  if (!assessment$paramvec_start){
    additionalargsbuilder <- function(extraargs = list(), paramvec_start = NULL){
      extraargs$paramvec_start <- NULL
      return(extraargs) #paramvec passed as part of extraargs
    }
  }
  if (assessment$paramvec_start){
    additionalargsbuilder <- function(extraargs = list(), paramvec_start = NULL){
      extraargs$paramvec_start <- paramvec_start #overwrites or adds a new element to the argument list
      return(extraargs) #paramvec passed as part of extraargs
    }
  }

  # define the function that extracts the estimated parameter value
  getparamfun <- extract_paramvec_fun(assessment$estlocation)

  #############
  # build the function that takes a theta and returns a new theta, depending on assessment results
  #############

  fpiterator <- function(fitted){
      stopifnot(length(fitted) == sum(!isfixed))
      fulltheta <- t_sfi2u(fitted, starttheta, isfixed) #including fitted and non-fitted parameter elements
      previous <- fulltheta
      weight_vec <- Windham_weights(ldenfun = ldenfun, Y = Y,
                 theta = fulltheta, cW)

      #calculate estimate:
      args = c(list(Y = Y, w = weight_vec), additionalargsbuilder(extraargs, fulltheta)) #paramvec passed
      estobj <- do.call(estimator, args = args)
      estparamvec <- getparamfun(estobj) #extract result
      #### adjust the estimates (Step 4 in Notes5.pdf)
      estparamvec <- thetaadjuster(estparamvec, previous, cW, cWav) #for WindhamCorrections() only estparamvec is used
      fitted <- t_si2f(estparamvec, isfixed)
      return(fitted)
  }

  #try fpiterator once
  tmp <- fpiterator(starttheta[!isfixed])

  # do the main computation, first check control args
  noncontrolargs <- setdiff(names(fpcontrol), methods::formalArgs(FixedPoint::FixedPoint)) 
  if (length(noncontrolargs) > 0){
    stop(paste("The follow arguments are not accepted by FixedPoint():", paste0(noncontrolargs, collapse = ", ")))
  }
  if (!isTRUE(fpcontrol$Method == "Simple")){warning("You have chosen to use a fixed point search that isn't the standard 'Simple'")}
  est <- do.call(FixedPoint::FixedPoint, 
                 c(list(Function = fpiterator, Inputs = starttheta[!isfixed]),
                    fpcontrol))

  # process results
  nevals <- ncol(est$Inputs)
  theta <- starttheta
  theta[!isfixed] <- est$FixedPoint

  # get weights corresponding the final iteration
  thetaprevious <- t_sfi2u(est$Inputs[,ncol(est$Inputs)], starttheta, isfixed) #including fitted and non-fitted parameter elements
  weight_vec <- Windham_weights(ldenfun = ldenfun, Y = Y,
                                theta = thetaprevious, cW)

  if (mean(weight_vec < 1E-10 / nrow(Y)) > 0.1){
    warning("More than 10% of weights are extremely small (smaller than 1E-10 / nrow(Y)) and are being treated like outliers. The fixed point search may have gone in an extreme direction. Are there too many parameters in the model?")
  }

  return(list(paramvec = theta,
           optim = c(est, list(finalweights = weight_vec))))
}

#' @title Inverse Transform for the Population Parameters Under Windham Weights
#' @description Returns the matrix which reverses the effect of weights on a population for certain models.
#' @param cW A vector of tuning constants for the Windham robustification method performed by [`Windham()`].
#' @return A diagonal matrix with the same number of columns as `cW`.
#' @details 
#' In the Windham robustification method ([`Windham()`]) the effect of weighting a population plays a central role.
#' When the 
#' the model density is proportional to \eqn{\exp(\eta(\theta) \cdot T(u))},
#' where \eqn{T(u)} is a vector of sufficient statistics for a measurement \eqn{u},
#' and \eqn{\eta} is a *linear* function,
#' Then weights proportional to 
#' \eqn{\exp(\eta(c \circ \theta) \cdot t(u))},
#' where \eqn{c} is a vector of tuning constants and \eqn{\circ} is the Hadamard (element-wise) product,
#' have a very simple effect on the population parameter vector \eqn{\theta}:
#' the weighted population follows a density of the same form, but with a parameter vector of 
#' \eqn{(1 + c) \circ \theta}.
#' The inverse of this change to the parameter vector is then a matrix multiplication by a diagonal matrix with elements \eqn{1/(1+c_i)}, with \eqn{c_i} denoting the elements of \eqn{c}.
#' @name Windham_populationinverse 
NULL

#' @describeIn Windham_populationinverse The matrix with diagonal elements \eqn{1/(1+c_i)} 
#' @export
Windham_populationinverse <- function(cW){
  tauinv <- diag(1/(1 + cW), nrow = length(cW)) #matrix that converts theta to the new theta*cW based on inclusion/exclusion  #klh: the extra argument nrow = length(cW) forces diag() to use the cW values on the diagonal, rather than treat them as the size of the matrix desired - useful when cW is legitimately length 1
  return(tauinv)
}

#' @describeIn Windham_populationinverse The transform implemented as described by \insertCite{scealy2024ro;textual}{scorematchingad}. It is mathematically equivalent to multiplication by the result of `Windham_populationinverse()` in the situation in \insertCite{scealy2024ro;textual}{scorematchingad}.
#' @export
#' @param newtheta The parameter vector most recently estimated
#' @param previoustheta The parameter vector estimated in the previous step
#' @param cWav The value of the non-zero elements of `cW`. That is `cW` have elements that are zero or equal to `cWav`.
Windham_populationinverse_alternative <- function(newtheta, previoustheta, cW, cWav){ #cW is a vector, cWav is the average of the non-zero elements of cW
  # generate the tuning constants dbeta, dA
  dtheta <- previoustheta * (cW - cWav) # = theta * cWav * (inWW - 1) = theta * cWav * -1 * !inWW = -cW * theta * !inWW
  return((newtheta - dtheta)/(cWav + 1))
}


