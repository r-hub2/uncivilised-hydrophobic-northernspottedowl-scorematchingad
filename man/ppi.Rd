% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppi.R, R/ppi_smvalues.R
\name{ppi}
\alias{ppi}
\alias{ppi_smvalues}
\title{Estimation of Polynomially-Tilted Pairwise Interaction (PPI) Model}
\usage{
ppi(
  Y,
  paramvec = NULL,
  trans,
  method = "closed",
  w = rep(1, nrow(Y)),
  constrainbeta = FALSE,
  bdryw = "ones",
  acut = NULL,
  bdrythreshold = 1e-10,
  shiftsize = bdrythreshold,
  approxorder = 10,
  control = list(tol = 1e-15, checkgrad = TRUE),
  paramvec_start = NULL
)

ppi_smvalues(
  Y,
  paramvec = NULL,
  evalparam,
  trans,
  method = "closed",
  w = rep(1, nrow(Y)),
  bdryw = "ones",
  acut = NULL,
  bdrythreshold = 1e-10,
  shiftsize = bdrythreshold,
  approxorder = 10,
  average = TRUE
)
}
\arguments{
\item{Y}{A matrix of measurements. Each row is a compositional measurement (i.e. each row sums to 1 and has non-negative elmenents).}

\item{paramvec}{Optionally a vector of the PPI models parameters. \code{NA}-valued elements of this vector are estimated and non-\code{NA} values are fixed. Generate \code{paramvec} easily using \code{\link[=ppi_paramvec]{ppi_paramvec()}}.  If \code{NULL} then all elements of \eqn{A_L}, \eqn{b_L} and \eqn{\beta} are estimated.}

\item{trans}{The name of the transformation of the manifold in Hyvärinen divergence (See \code{\link{scorematchingtheory}}): "clr" (centred log ratio), "alr" (additive log ratio), "sqrt" or "none".}

\item{method}{\code{"closed"} uses \code{CppAD} to solve in closed form the a quadratic score matching discrepancy using \code{\link[=cppad_closed]{cppad_closed()}}. \code{"hardcoded"} uses hardcoded implementations. "iterative" uses \code{\link[=cppad_search]{cppad_search()}} (which uses \code{CppAD} and \code{\link[optimx:Rcgmin]{optimx::Rcgmin()}}) to iteratively find the minimum of the weighted Hyvärinen divergence.}

\item{w}{Weights for each observation, if different observations have different importance. Used by \code{\link[=Windham]{Windham()}} and \code{\link[=ppi_robust]{ppi_robust()}} for robust estimation.}

\item{constrainbeta}{If \code{TRUE}, elements of \eqn{\beta} that are less than \code{-1} are converted to \code{-1 + 1E-7}.}

\item{bdryw}{The boundary weight function for down weighting measurements as they approach the manifold boundary. Either "ones", "minsq" or "prodsq". See details.}

\item{acut}{The threshold \eqn{a_c} in \code{bdryw} to avoid over-weighting measurements interior to the simplex}

\item{bdrythreshold}{\code{iterative} or \code{closed} methods only. For measurements within \code{bdrythreshold} of the simplex boundary a Taylor approximation is applied by shifting the measurement \code{shiftsize} towards the center of the simplex.}

\item{shiftsize}{\code{iterative} or \code{closed} methods only. For measurements within \code{bdrythreshold} of the simlex boundary a Taylor approximation is applied by shifting the measurement \code{shiftsize} towards the center of the simplex.}

\item{approxorder}{\code{iterative} or \code{closed} methods only. Order of the Taylor approximation for measurements on the boundary of the simplex.}

\item{control}{\code{iterative} only. Passed to \code{\link[optimx:Rcgmin]{optimx::Rcgmin()}} to control the iterative solver.}

\item{paramvec_start}{\code{iterative} method only. The starting guess for \code{Rcgmin}. Generate \code{paramvec_start} easily using \code{\link[=ppi_paramvec]{ppi_paramvec()}}.}

\item{evalparam}{The parameter set to evaluate the score matching values.
This is different to \code{paramvec}, which specifies which parameters to estimate.
All elements of \code{evalparam} must be non-NA, and any parameters fixed by \code{paramvec} must have the same value in \code{evalparam}.}

\item{average}{If TRUE return the (weighted average) of the measurements, otherwise return the values for each measurement.}
}
\value{
\code{ppi()} returns:
A list of \code{est}, \code{SE} and \code{info}.
\itemize{
\item \code{est} contains the estimates in vector form, \code{paramvec}, and as \eqn{A_L}, \eqn{b_L} and \eqn{\beta}.
\item \code{SE} contains estimates of the standard errors if computed. See \code{\link[=cppad_closed]{cppad_closed()}}.
\item \code{info} contains a variety of information about the model fitting procedure and results.
}

\code{ppi_smvalues()} returns a list of
\itemize{
\item \code{obj} the score matching discrepancy value
\item \code{grad} the gradient of the score matching discrepancy
\item \code{hess} the Hessian of the score matching discrepancy
\item \code{offset} gradient offset (see \code{\link[=quadratictape_parts]{quadratictape_parts()}})
}
}
\description{
Estimates the parameters of the Polynomially-Tilted Pairwise Interaction (PPI) model \insertCite{scealy2023sc}{scorematchingad} for compositional data.
By default \code{ppi()} uses \code{\link[=cppad_closed]{cppad_closed()}} to find estimate.
For many situations a hard-coded implementation of the score matching estimator is also available.

For a given parameter vector \code{evalparam}, \code{ppi_smvalues()} computes the score matching discrepancy, the gradient and the Hessian of the score matching discrepancy (see \code{\link[=smvalues]{smvalues()}}) and the gradient offset of the score matching discrepancy (see \code{\link[=quadratictape_parts]{quadratictape_parts()}} and \code{\link[=tapeGradOffset]{tapeGradOffset()}}).
}
\details{
Estimation may be performed via transformation of the measure in Hyvärinen divergence from Euclidean space to the simplex (inverse of the additive log ratio transform), from a hyperplane to the simplex (inverse of the centred log ratio transform), from the positive quadrant of the sphere to the simplex (inverse of the square root transform), or without any transformation. In the latter two situations there is a boundary and \emph{weighted Hyvärinen divergence} \insertCite{@Equation 7, @scealy2023sc}{scorematchingad} is used.
Properties of the estimator using the square root transform were studied by \insertCite{scealy2023sc;textual}{scorematchingad}.
Properties of the estimator using the additive log ratio transfrom were studied by \insertCite{scealy2024ro;textual}{scorematchingad}.

There are three boundary weight functions available:
\itemize{
\item The function "ones" applies no weights and should be used whenever the manifold does not have a bounday.
\item The function "minsq" is the minima-based boundary weight function for the PPI model \insertCite{@Equation 12, @scealy2023sc}{scorematchingad}
\deqn{\tilde{h}(z)^2 = \min(z_1^2, z_2^2, ..., z_p^2, a_c^2).}{h(z)^2 = min(z1^2, z2^2, ..., zp^2, a_c^2),}
where \eqn{z} is a point in the positive orthant of the p-dimensional unit sphere
and \eqn{z_j}{zj} is the jth component of z.
\item The function "prodsq" is the product-based \insertCite{@Equation 9, @scealy2023sc}{scorematchingad}
\deqn{\tilde{h}(z)^2 = \min(\prod_{j=1}^{p} z_j^2, a_c^2).}{h(z)^2 = min(z1^2 * z2^2 * ... * zp^2, a_c^2),}
where \eqn{z} is a point in the positive orthant of the p-dimensional unit sphere
and \eqn{z_j}{zj} is the jth component of z.
}

Scealy and Wood \insertCite{@Theorem 1, @scealy2023sc}{scorematchingad} prove that minimising the weighted Hyvärinen Divergence is equivalent to minimising \eqn{\psi(f, f_0)} (See \code{\link{scorematchingtheory}})
when the boundary weight function is smooth or for the functions "minsq" and "prodsq"  above when the manifold is the simplex or positive orthant of a sphere.

Hard-coded estimators are available for the following situations:
\itemize{
\item Square root transformation ("sqrt") with the "minsq" boundary weight function:
\itemize{
\item full parameter vector (\code{paramvec} not provided)
\item \code{paramvec} fixes only the final element of \eqn{\beta}
\item \code{paramvec} fixes all elements of \eqn{\beta}
\item \code{paramvec} fixes \eqn{b_L = 0} and provides fixed values of \eqn{\beta}
\item \code{paramvec} fixes \eqn{A_L=0} and \eqn{b_L=0}, leaving \eqn{\beta} to be fitted.
}
\item Square root transformation ("sqrt") with the "prodsq" boundary weight function:
\itemize{
\item \code{paramvec} fixes all elements of \eqn{\beta}
\item \code{paramvec} fixes \eqn{b_L = 0} and provides fixed values of \eqn{\beta}
\item \code{paramvec} fixes \eqn{A_L=0} and \eqn{b_L=0}, leaving \eqn{\beta} to be fitted.
}
\item The additive log ratio transformation ("alr") using the final component on the denominator, with \eqn{b_L=0} and fixed final component of \eqn{\beta}.
}
}
\section{PPI Model}{

The PPI model density is proportional to
\deqn{\exp(z_L^TA_Lz_L + b_L^Tz_L)\prod_{i=1}^p z_i^{\beta_i},}
where \eqn{p} is the dimension of a compositional measurement \eqn{z}, and \eqn{z_L} is \eqn{z} without the final (\eqn{p}th) component.
\eqn{A_L} is a \eqn{p-1 \times p-1} symmetric matrix that controls the covariance between components.
\eqn{b_L} is a \eqn{p-1} vector that controls the location within the simplex.
The \eqn{i}th component \eqn{\beta_i} of \eqn{\beta} controls the concentration of density when \eqn{z_i} is close to zero: when \eqn{\beta_i \geq 0} there is no concentration and \eqn{\beta_i} is hard to identify; when \eqn{\beta_i < 0} then the probability density of the PPI model increases unboundedly as \eqn{z_i} approaches zero, with the increasing occuring more rapidly and sharply the closer \eqn{\beta_i} is to \eqn{-1}.
}

\examples{
model <- rppi_egmodel(100)
estalr <- ppi(model$sample,
              paramvec = ppi_paramvec(betap = -0.5, p = ncol(model$sample)),
              trans = "alr")
estsqrt <- ppi(model$sample,
              trans = "sqrt",
              bdryw = "minsq", acut = 0.1)
}
\references{
\insertAllCited{}
}
\seealso{
Other PPI model tools: 
\code{\link{dppi}()},
\code{\link{ppi_param_tools}},
\code{\link{ppi_robust}()},
\code{\link{rppi}()}
}
\concept{PPI model tools}
