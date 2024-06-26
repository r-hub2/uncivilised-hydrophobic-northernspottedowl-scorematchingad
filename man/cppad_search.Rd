% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cppad_search.R
\name{cppad_search}
\alias{cppad_search}
\title{Iterative Score Matching Estimator Using Conjugate-Gradient Descent}
\usage{
cppad_search(
  smdtape,
  theta,
  Y,
  Yapproxcentres = NA * Y,
  w = rep(1, nrow(Y)),
  approxorder = 10,
  control = list(tol = 1e-15, checkgrad = TRUE)
)
}
\arguments{
\item{smdtape}{A \code{CppAD} tape of a score matching discrepancy function that has \emph{quadratic form}. Test for quadratic form using \code{\link[=testquadratic]{testquadratic()}}.
The \code{smdtape}'s independent variables are assumed to be the model parameters to fit
and the \code{smdtape}'s dynamic parameter is a (multivariate) measurement.}

\item{theta}{The starting parameter set}

\item{Y}{A matrix of multivariate observations. Each row is an observation.}

\item{Yapproxcentres}{A matrix of Taylor approximation centres for rows of Y that require approximation. \code{NA} for rows that do not require approximation.}

\item{w}{Weights for each observation.}

\item{approxorder}{The order of Taylor approximation to use.}

\item{control}{Control parameters passed to \code{\link[optimx:Rcgmin]{optimx::Rcgmin()}}}
}
\description{
Uses conjugate gradient descent to search for a vector of parameters such that gradient of the score matching discrepancy is within tolerance of zero.
Also estimates standard errors and covariance.
}
\details{
The score matching discrepancy function and gradient of the score matching function are passed to \code{\link[optimx:Rcgmin]{optimx::Rcgmin()}}.
The call to \code{\link[optimx:Rcgmin]{optimx::Rcgmin()}} uses the \emph{sum} of observations (as opposed to the mean) to reduce floating point inaccuracies. This has implications for the meaning of the control parameters passed to \code{Rcgmin()} (e.g. \code{tol}). The results are converted into averages so the use of sums can be ignored when not setting control parameters, or studying the behaviour of Rcgmin.

Standard errors use the Godambe information matrix (aka sandwich method) and are only computed when the weights are constant.
The estimate of the sensitivity matrix \eqn{G} is
the negative of the average over the Hessian of \code{smdtape} evaluated at each observation in \code{Y}.
The estimate of the variability matrix \eqn{J} is then
the sample covariance (denominator of \eqn{n-1}) of the gradiant of \code{smdtape} evaluated at each of the observations in \code{Y} for the estimated \eqn{\theta}.
The variance of the estimator is then estimated as
\eqn{G^{-1}JG^{-1}/n,}
where \code{n} is the number of observations.

Taylor approximation is available because boundary weight functions and transformations of the measure in Hyvärinen divergence can remove singularities in the model log-likelihood, however evaluation at these singularities may still involve computing intermediate values that are unbounded.
If the singularity is ultimately removed, then Taylor approximation from a nearby location will give a very accurate evaluation at the removed singularity.
}
\section{Warning}{
There appears to be floating point issues with evaluation of the gradient leading to slow or no convergence in many situations. Tweaking the convergence tolerance \code{tol} may be appropriate. If the closed form solution exists please use \code{cppad_closed()} instead.
}

\examples{
smdtape <- buildsmdtape("sim", "sqrt", "sph", "ppi",
              ytape = rep(1/3, 3),
              usertheta = ppi_paramvec(p = 3),
              bdryw = "minsq", acut = 0.01,
              verbose = FALSE
              )$smdtape
Y <- rppi_egmodel(100)
cppad_search(smdtape, 0.9 * Y$theta, Y$sample)
sum((smvalues_wsum(smdtape, Y$sample, Y$theta)$grad/nrow(Y$sample))^2)
}
\seealso{
Other generic score matching tools: 
\code{\link{Windham}()},
\code{\link{cppad_closed}()}
}
\concept{generic score matching tools}
