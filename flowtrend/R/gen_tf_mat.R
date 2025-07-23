# Generated from _main.Rmd: do not edit by hand

#' A lasso regressor matrix H that can be used to solve an equivalent problem as the trend filtering of the \code{k}'th degree.
#'
#' @param n Total number of time points.
#' @param k Degree of trend filtering for cluster probabilities. $k=0$ is fused lasso, $k=1$ is linear trend filtering, and so on.
#' @param x Time points
#' 
#' @return $n$ by $n$ matrix.
#' 
#' @export
gen_tf_mat <- function(n, k, x = NULL){

  if(is.null(x) ){
    x = (1:n)/n
  }
  if(!is.null(x)){
    stopifnot(length(x) == n)
  }

  ## For every i,j'th entry, use this helper function (from eq 25 of Tibshirani
  ## (2014)).
  gen_ij <- function(x, i, j, k){
    xi <- x[i]
    if(j %in% 1:(k+1)){
      return(xi^(j-1))
    }
    if(j %in% (k+2):n){

      ## Special handling for k==0, See lemma 4 eq 25
      if(k == 0){
        prd = 1 
        ind = j
      }
      if(k >= 1){
        ind = j - (1:k)
        prd = prod(xi - x[ind]) 

      }
      return(prd * ifelse(xi >= x[max(ind)], 1, 0))
    }
  }

  ## Generate the H matrix, entry-wise.
  H <- matrix(nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      H[i,j] <- gen_ij(x, i,j, k)
    }
  }
  return(H)
}
