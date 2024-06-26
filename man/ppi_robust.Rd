% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppi_robust.R, R/ppi_robust_alrgengamma.R
\name{ppi_robust}
\alias{ppi_robust}
\alias{ppi_robust_alrgengamma}
\title{Robustly Estimate Parameters of the PPI Distribution}
\usage{
ppi_robust(Y, cW, ...)

ppi_robust_alrgengamma(
  Y,
  cW,
  ...,
  fpcontrol = list(Method = "Simple", ConvergenceMetricThreshold = 1e-10)
)
}
\arguments{
\item{Y}{A matrix of measurements. Each row is a measurement, each component is a dimension of the measurement.}

\item{cW}{A vector of robustness tuning constants. Easy to build using \code{\link[=ppi_cW]{ppi_cW()}} and \code{\link[=ppi_cW_auto]{ppi_cW_auto()}}. See \code{\link[=Windham]{Windham()}} for more details on \code{cW}.}

\item{...}{Passed to \code{\link[=Windham]{Windham()}} and on to \code{\link[=ppi]{ppi()}}.}

\item{fpcontrol}{A named list of control arguments to pass to \code{\link[FixedPoint:FixedPoint]{FixedPoint::FixedPoint()}} for the fixed point iteration.}
}
\value{
A list:
\itemize{
\item \code{est} The estimated parameters in vector form (\code{paramvec}) and as \code{AL}, \code{bL} and \code{beta}.
\item \code{SE} "Not calculated." Returned for consistency with other estimators.
\item \code{info} Information returned in the \code{optim} slot of \code{\link[=Windham]{Windham()}}. Includes the final weights in \code{finalweights}.
}
}
\description{
\code{ppi_robust()} uses \code{\link[=Windham]{Windham()}} and \code{\link[=ppi]{ppi()}} to estimate a PPI distribution robustly.
There are many arguments to the \code{\link[=ppi]{ppi()}} function and we highly recommend trialling your arguments on \code{\link[=ppi]{ppi()}} first before running \code{ppi_robust()}.

\code{ppi_robust_alrgengamma()} performs the Windham robustification algorithm exactly as described in \insertCite{scealy2024ro;textual}{scorematchingad} for score matching via log-ratio transform of the PPI model with \eqn{b_L = 0}. This function calls the more general \code{\link[=Windham]{Windham()}} and \code{\link[=ppi]{ppi()}}.
}
\details{
\code{ppi_robust_alrgengamma()}: must fit a PPI model via additive-log ratio transform of the simplex with \eqn{b_L=0} fixed and the final element of \eqn{\beta} fixed.
The default convergence metric and threshold are different to the default for \code{\link[=ppi_robust]{ppi_robust()}} to match the implementation in \insertCite{scealy2024ro}{scorematchingad}: convergence is measured by the change in the first element of \eqn{\beta}, and convergence is reached when the change is smaller than \code{1E-6}. Override this behaviour by specifying the elements \code{ConvergenceMetric} and \code{ConvergenceMetricThreshold} in a list passed as \code{fpcontrol}.
\code{\link[=Windham]{Windham()}} is called with \code{alternative_populationinverse = TRUE}.
}
\examples{
model <- rppi_egmodel(100)
estsqrt <- ppi_robust(model$sample,
  cW = ppi_cW_auto(0.01, model$sample),
  trans = "sqrt", bdryw = "minsq", acut = 0.1)
model <- rppi_egmodel(100)
ppi_robust_alrgengamma(model$sample,
   cW = ppi_cW_auto(0.01, model$sample),
   paramvec = ppi_paramvec(betap = -0.5, p = ncol(model$sample)))
}
\references{
\insertAllCited{}
}
\seealso{
Other PPI model tools: 
\code{\link{dppi}()},
\code{\link{ppi}()},
\code{\link{ppi_param_tools}},
\code{\link{rppi}()}

Other Windham robustness functions: 
\code{\link{Windham}()},
\code{\link{vMF_robust}()}
}
\concept{PPI model tools}
\concept{Windham robustness functions}
