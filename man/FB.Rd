% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FB.R
\name{FB}
\alias{FB}
\title{Estimate the Fisher-Bingham Distribution}
\usage{
FB(Y, km = NULL, A = NULL)
}
\arguments{
\item{Y}{A matrix of multivariate observations in Cartesian coordinates. Each row is a multivariate measurement (i.e. each row corresponds to an individual).}

\item{km}{Optional. A vector of the same length as the dimension, representing the parameter vector for the von Mises-Fisher component (i.e. the \eqn{\kappa \mu} see \code{\link[=vMF]{vMF()}}).
If supplied, the non-NA elements are fixed.}

\item{A}{Optional. The Bingham matrix. If supplied the non-NA elements of the Bingham matrix are fixed.
The final element of the diagonal of \code{A} must be NA as the software calculates this value to ensure the trace of the Bingham matrix is zero.}
}
\description{
Estimates parameters for the Fisher-Bingham distribution using score-matching.
}
\details{
The density of the Fisher-Bingham distribution is proportional to
\deqn{\exp(z^TAz + \kappa\mu^Tz),}
where \eqn{A} is a matrix as in the Bingham distribution, and
\eqn{\kappa} and \eqn{\mu} are the concentration and mean direction, respectively, as in the von Mises-Fisher distribution.
}
\section{Warning: Slow Convergence with Sample Size}{
Score matching estimates of all elements of \eqn{A} and \eqn{\kappa\mu} converge slowly with sample size.
Even with a million simulated measurements,
the gradient of the score matching discrepancy at the true parameters can have size (L2 Euclidean norm) more than 0.001, which is substantially non-zero.
}

\examples{
p <- 3
A <- rsymmetricmatrix(p, -10, 10)
A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
m <- runif(p, -10, 10)
m <- m / sqrt(sum(m^2))
if (requireNamespace("simdd")){
  Y <- simdd::rFisherBingham(1000, 2 * m, A)
  FB(Y)
}
}
\seealso{
Other directional model estimators: 
\code{\link{Bingham}()},
\code{\link{vMF}()},
\code{\link{vMF_robust}()}
}
\concept{directional model estimators}
