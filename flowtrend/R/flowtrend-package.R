#' flowtrend
#'
#' This package implements the `flowtrend` method for automatic gating of flow cytometry data using trend filtering.
#' It was proposed in <https://arxiv.org/abs/2504.12287>.  To learn
#' more about this package, please visit its website
#' <https://sangwon-hyun/flowtrend-project>.
#' 
#' @docType package
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import RcppEigen
#' @importFrom Rcpp sourceCpp
#' @useDynLib flowtrend, .registration = TRUE
## usethis namespace: end
NULL
