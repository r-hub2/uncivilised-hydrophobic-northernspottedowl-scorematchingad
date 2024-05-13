#' @importFrom stats weighted.mean
#' @importFrom Rdpack reprompt
#' @useDynLib scorematchingad
#' @details
#' This package's main features are
#'  * A general capacity to implement score matching estimators that use algorithmic differentiation to avoid tedious manual algebra.
#' The package uses `CppAD` and `Eigen` to differentiate model densities and compute the score matching discrepancy function (see [`scorematchingtheory`]).
#' The score matching discrepancy is usually minimised by solving a quadratic equation, but a method for solving numerically (through [`optimx::Rcgmin()`]) is also included.
#' On Linux platforms using the `gcc` compiler new models can be fitted with the help of [`customll()`], in a similar fashion to models in the `TMB` package.
#' New manifolds or new transforms require small alterations to the source code of this package.
#'  * Score matching estimators for the Polynomially-Tilted Pairwise Interaction (PPI) model \insertCite{scealy2023sc,scealy2024ro}{scorematchingad}. See function [`ppi()`].
#'  * Score matching and hybrid score matching estimators for von Mises Fisher, Bingham and Fisher-Bingham directional distributions \insertCite{mardia2016sc}{scorematchingad}. See [`vMF()`], [`Bingham()`] and [`FB()`].
#'  * Implementation of a modification of Windham's robustifying method \insertCite{windham1995ro}{scorematchingad} for many exponential family distributions. See [`Windham()`].
#' For some models the density approaches infinity at some locations, creating difficulties for the weights in Windham's original method \insertCite{scealy2024ro}{scorematchingad}.
#' \insertNoCite{*}{scorematchingad}
#' # Acknowledgements
#' Colleagues Andrew T. A. Wood and John T. Kent played important roles in developing the statistical ideas and theory for score matching estimation for the PPI model \insertCite{scealy2024ro}{scorematchingad}.
#' 
#' We developed this package on Ngunnawal and Ngambri Country. We thank the Country for its influence.
#' @references
#' \insertAllCited{}
"_PACKAGE"

