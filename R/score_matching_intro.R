#' @title Introduction to Score Matching
#' @name scorematchingtheory
#' 
#' @description
#' This package includes score matching estimators for particular distributions and a general capacity to implement additional score matching estimators.
#' Score matching is a popular estimation technique when normalising constants for the proposed model are difficult to calculate or compute.
#' Score matching was first developed by \insertCite{hyvarinen2005es;textual}{scorematchingad} and was further developed for subsets of Euclidean space \insertCite{hyvarinen2007ex,yu2019ge,yu2020ge,liu2019es}{scorematchingad}, Riemannian manifolds \insertCite{mardia2016sc,mardia2018ne}{scorematchingad},
#'  and Riemannian manifolds with boundary \insertCite{scealy2023sc}{scorematchingad}.
#' In this help entry we briefly describe score matching estimation.
#' 
#' @details
#' # Score Matching in General
#' In the most general form (Riemannian manifolds with boundary) score matching minimises the weighted Hyvärinen divergence \insertCite{@Equation 7, @scealy2023sc}{scorematchingad} 
#' \deqn{
#' \phi(f, f_0) =  \frac{1}{2} \int_M f_0(z)h(z)^2 \left\lVert P(z)\Big(\nabla_z \log(f) - \nabla_z \log(f_0)\Big)\right\rVert^2 dM(z),
#' }
#'  where
#' 
#'  + \eqn{M} is the manifold, isometrically embedded in Euclidean space, and \eqn{dM(z)} is the unnormalised uniform measure on \eqn{M}.
#'  + \eqn{P(z)} is the matrix that projects points onto the tangent space of the manifold at \eqn{z}, which is closely related to to Riemannian metric of \eqn{M}.
#'  + \eqn{f_0} is the density of the data-generating process, defined with respect to \eqn{dM(z)}.
#'  + \eqn{f} is the density of a posited model, again defined with respect to \eqn{dM(z)}.
#'  + \eqn{h(z)} is a function, termed the *boundary weight function*, that is zero on the boundary of \eqn{M} and smooth \insertCite{@Section 3.2, @scealy2023sc}{scorematchingad} or potentially piecewise smooth.
#'  + \eqn{\nabla_z} is the Euclidean gradient operator that differentiates with respect to \eqn{z}.
#'  + \eqn{\lVert \cdot \rVert} is the Euclidean norm.
#'  
#'  Note that, because \eqn{P(z)} is the projection matrix, \eqn{\left\lVert P(z)\Big(\nabla_z \log(f) - \nabla_z \log(f_0)\Big)\right\rVert^2} is the natural inner product of the gradient of the log ratio of \eqn{f} and \eqn{f_0}.
#' 
#'  When the density functions \eqn{f} and \eqn{f_0} are smooth and positive inside \eqn{M},
#'  and the boundary weight function is smooth or of particular forms for specific manifolds \insertCite{@Section 3.2, @scealy2023sc}{scorematchingad},
#'  then minimising the weighted Hyvärinen divergence \eqn{\phi(f, f_0)} is equivalent to minimising the score matching discrepancy \insertCite{@Theorem 1, @scealy2023sc}{scorematchingad}
#' \deqn{
#' \psi(f, f_0) = \int f_0(z)\big(A(z) + B(z) + C(z)\big)dM(z),
#' }
#'  where 
#'  \deqn{A(z) = \frac{1}{2} h(z)^2 \left(\nabla_z \log(f)\right)^T P(z) \left(\nabla_z \log(f)\right),}
#'  \deqn{B(z) = h(z)^2\Delta_z\log(f),}
#'  \deqn{C(z) = \left(\nabla_z h(z)^2)\right)^T P(z) \left(\nabla_z \log(f)\right),}
#' and \eqn{\Delta} is the Laplacian operator.
#' We term \deqn{A(z) + B(z) + C(z)} the *score matching discrepancy function*.
#' 
#' We suspect that  \insertCite{@Theorem 1, @scealy2023sc}{scorematchingad} holds more generally for nearly all realistic continuous and piecewise-smooth boundary weight functions, although no proof exists to our knowledge.
#' 
#'  When \eqn{n} independent observations from \eqn{f_0} are available, the integration in \eqn{\psi(f, f_0)} can be approximated by an average over the observations, 
#'   \deqn{\psi(f, f_0) \approx \hat\psi(f, f_0) = \frac{1}{n} \sum_{i = 1}^n A(z_i) + B(z_i) + C(z_i).}
#' 
#' If we parameterise a family of models \eqn{f_\theta} according to a vector of parameters \eqn{\theta}, then the *score matching estimate* is the \eqn{\theta} that minimises \eqn{\hat\psi(f_\theta, f_0)}.
#' In general, the score matching estimate must be found via numerical optimisation techniques, such as in the function `cppad_search()`.
#' However, when the family of models is a canonical exponential family then often \eqn{\hat\psi(f_\theta, f_0)} and the score matching discrepancy function is a quadratic function of \eqn{\theta} \insertCite{mardia2018ne}{scorematchingad} and the minimum has a closed-form solution found by `cppad_closed()`.
#' 
#' Note that when \eqn{M} has a few or more dimensions, the calculations of \eqn{A(z)}, \eqn{B(z)} and \eqn{C(z)} can become cumbersome. This package uses `CppAD` to automatically compute \eqn{A(z)}, \eqn{B(z)} and \eqn{C(z)}, and the quadratic simplification if it exists.
#' 
#' # Transformations
#' Hyvärinen divergence \eqn{\phi(f, f_0)} is sensitive to transformations of the measure \eqn{dM(z)}, including transforming the manifold.
#' That is, transforming the manifold \eqn{M} changes the divergence between distributions and changes the minimum of \eqn{\hat\psi(f_\theta, f_0)}.
#' The transformation changes measure \eqn{dM(z)}, the divergence and the estimator but does *not* transform the data.
#' 
#' For example, many different transformations of the simplex (i.e. compositional data) are possible \insertCite{@Appendix A.3, @scealy2024ro}{scorematchingad}.
#' Hyvärinen divergences that use the sphere, obtained from the simplex by a square root, have different behaviour to Hyvärinen divergence using Euclidean space obtained from the simplex using logarithms \insertCite{scealy2024ro}{scorematchingad}.
#' The estimator for the latter does not apply logarithms to the observations, in fact the estimator involves only polynomials of the observed compositions \insertCite{scealy2024ro}{scorematchingad}.
#' 
#' The variety of estimator behaviour available through different transformations was a major motivator for this package as each transformation has different \eqn{A(z)}, \eqn{B(z)} and \eqn{C(z)}, and without automatic differentiation, implementation of the score matching estimator in each case would require a huge programming effort.
#' 
#' @references \insertAllCited{}
NULL

