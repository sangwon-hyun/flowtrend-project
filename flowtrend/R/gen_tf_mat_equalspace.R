# Generated from _main.Rmd: do not edit by hand

#' Creates trend filtering regression matrix using Lemma 2 of Ryan Tibshirani's
#' trendfilter paper (2014); works on equally spaced data only.
#'
#' @param n Number of data points.
#' @param k Order of trend filter. 0 is fused lasso, and so on.
#' 
#' @examples 
#' ord = 1
#' H_tf <- gen_tf_mat_equalspace(n = 100, k = ord)
#' H_tf[,1] * 100
#' H_tf[,2] * 100
#' H_tf[,3] * 100
#' H_tf[,4] * 100
#' H_tf[,99] * 100
#' H_tf[,100] * 100
#' 
#' @return n by n matrix.
gen_tf_mat_equalspace <- function(n, k){
  nn = n
  kk = k
  ##' Connects kk to kk-1.
  sigm <- function(ii, kk){
    if(kk == 0) return(rep(1, ii))
    cumsum(sigm(ii, kk-1))
  }

  mat = matrix(NA, ncol = nn, nrow = nn)
  for(ii in 1:nn){
    for(jj in 1:nn){
      if(jj <= kk+1) mat[ii,jj] = ii^(jj-1) / nn^(jj-1)
      if(ii <= jj-1 & jj >= kk+2) mat[ii, jj] = 0
      if(ii > jj-1 & jj >= kk+2){
        mat[ii, jj] = (sigm(ii-jj+1, kk) %>% tail(1)) * factorial(kk) / nn^kk
      }
    }
  }
  return(mat)
}
