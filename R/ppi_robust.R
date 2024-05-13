#' @title Robustly Estimate Parameters of the PPI Distribution
#' @order 1
#' @family PPI model tools
#' @description `ppi_robust()` uses [`Windham()`] and [`ppi()`] to estimate a PPI distribution robustly.
#' There are many arguments to the [`ppi()`] function and we highly recommend trialling your arguments on [`ppi()`] first before running `ppi_robust()`.
#' @param Y A matrix of measurements. Each row is a measurement, each component is a dimension of the measurement.
#' @param cW A vector of robustness tuning constants. Easy to build using [`ppi_cW()`] and [`ppi_cW_auto()`]. See [`Windham()`] for more details on `cW`.
#' @param ... Passed to [`Windham()`] then to [`ppi()`].
#' @family Windham robustness functions
#' @return A list:
#'  * `est` The estimated parameters in vector form (`paramvec`) and as `AL`, `bL` and `beta`.
#'  * `SE` "Not calculated." Returned for consistency with other estimators.
#'  * `info` Information returned in the `optim` slot of [`Windham()`]. Includes the final weights in `finalweights`.
#' @examples
#' model <- rppi_egmodel(100)
#' estsqrt <- ppi_robust(model$sample,
#'   cW = ppi_cW_auto(0.01, model$sample),
#'   trans = "sqrt", bdryw = "minsq", acut = 0.1)
#' @export
ppi_robust <- function(Y, cW, ...){
  ellipsis::check_dots_used()
  ldenfun <- function(Y, theta){ #here theta is the usual parameters of PPI model from
    return(drop(dppi(Y, paramvec = theta)))
  }

  est = Windham(Y = Y,
                    estimator = ppi,
                    ldenfun = ldenfun,
                    cW = cW,
                    ...)

  #make results nicer and consistent with ppi()
  out <- list(
    est = c(list(paramvec = est$paramvec), ppi_parammats(est$paramvec)),
    SE = "Not calculated.",
    info = est$optim
  )

  return(out)
}
