#' @noRd 
inlineCxxPlugin <- function(...) {
    plugin <-
        Rcpp::Rcpp.plugin.maker(
                  include.before = c('#include <scorematchingad.h>'),
                  include.after = paste("#include", c("<utils/PrintFor.hpp>", "<utils/wrapas.hpp>", "<likelihoods/likelihoods.hpp>"), collapse = "\n"),
                  package        = "scorematchingad",
                  LinkingTo = c("RcppEigen", "scorematchingad", "Rcpp")
              )
    settings <- plugin()
    settings
}

