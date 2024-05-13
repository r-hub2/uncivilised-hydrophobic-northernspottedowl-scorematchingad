#' @title Windham Robustness for Scealy et al 2024
#' @rdname ppi_robust
#' @order 2
#' @description
#' `ppi_robust_alrgengamma()` performs the Windham robustification algorithm exactly as described in \insertCite{scealy2024ro;textual}{scorematchingad} for score matching via log-ratio transform of the PPI model with \eqn{b_L = 0}. This function calls the more general [`Windham()`] and [`ppi()`].
#' @inheritParams Windham
#' @inheritParams ppi_robust
#' @inherit ppi_robust return
#' @param ... Passed to [`Windham()`] and on to [`ppi()`]. 
#' @details
#' `ppi_robust_alrgengamma()`: must fit a PPI model via additive-log ratio transform of the simplex with \eqn{b_L=0} fixed and the final element of \eqn{\beta} fixed.
#' The default convergence metric and threshold are different to the default for [`ppi_robust()`] to match the implementation in \insertCite{scealy2024ro}{scorematchingad}: convergence is measured by the change in the first element of \eqn{\beta}, and convergence is reached when the change is smaller than `1E-6`. Override this behaviour by specifying the elements `ConvergenceMetric` and `ConvergenceMetricThreshold` in a list passed as `fpcontrol`.
#' [`Windham()`] is called with `alternative_populationinverse = TRUE`.
#' @references
#' \insertAllCited{}
#' @examples
#' model <- rppi_egmodel(100)
#' ppi_robust_alrgengamma(model$sample,
#'    cW = ppi_cW_auto(0.01, model$sample),
#'    paramvec = ppi_paramvec(betap = -0.5, p = ncol(model$sample)))
#' @export
ppi_robust_alrgengamma <- function(Y, cW, ..., fpcontrol = list(Method = "Simple", ConvergenceMetricThreshold = 1E-10)){
  ldenfun <- function(Y, theta){ #here theta is the usual parameters of PPI model from
    return(drop(dppi(Y, paramvec = theta)))
  }

  JSConvergenceMetricThreshold = 1E-6 #used by Janice Scealy in original code
  # The following convergence metric matches the metric used by Janice Scealy in her original code: a threshold on the change in the first element of the beta estimate
  JSConvergenceMetric <- function(residuals){
    p <- (-1 + sqrt(1 + 8 + 8*length(residuals)))/2 #from quadratic formula: qstar = p(p-1)/2 + p-1
    return(abs(residuals[p*(p-1)/2 + 1])) #the first element after the AL is the first element of beta - this is to math
  }

  #construct fpcontrol list
  fpcontrol = c(fpcontrol, 
                      list(ConvergenceMetric = JSConvergenceMetric,
                      ConvergenceMetricThreshold = JSConvergenceMetricThreshold))
  # remove duplicates - removes the second element
  fpcontrol <- fpcontrol[!duplicated(names(fpcontrol))]

  est = Windham(Y = Y,
                    estimator = ppi,
                    ldenfun = ldenfun,
                    cW = cW,
                    trans = "alr",
                    ...,
                    fpcontrol = fpcontrol,
                    alternative_populationinverse = TRUE)

  #make results nicer and consistent with ppi()
  out <- list(
    est = c(list(paramvec = est$paramvec), ppi_parammats(est$paramvec)),
    SE = "Not calculated.",
    info = est$optim
  )

  return(out)
}
