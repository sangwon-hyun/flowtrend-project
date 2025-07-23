# Generated from _main.Rmd: do not edit by hand

#' a
#'
#' @param TT
etilde_mat <- function(TT){

  mats <- lapply(1:TT, FUN = function(t){
    e_vec <- rep(0, TT)
    e_vec[t] <- 1
    (e_vec - 1/TT) %*% t(e_vec - 1/TT)
  })
  Reduce('+', mats)
}
