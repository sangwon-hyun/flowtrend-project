# Generated from _main.Rmd: do not edit by hand

#' Solve a fused lasso problem for the W update. Internally, a fused lasso dynamic
#' programming solver \code{prox()} (which calls \code{prox_dp()} written in C)
#' is used.
#'
W_update_fused <- function(l, TT, mu, uw, rho, lambda, Dl){

  # modified lambda for fused lasso routine
  mod_lam <- lambda/rho

  # generating pseudo response xi
  if( l < 0 ){
    stop("l should never be /below/ zero!")
  } else if( l == 0 ){
    xi <- mu + 1/rho * uw  ## This is faster
  } else {
    xi <- Dl %*% mu + 1/rho * uw
    if(any(is.nan(xi))) browser()

    ## l = 2 is quadratic trend filtering
    ## l = 1 is linear trend filtering
    ## l = 0 is fused lasso
    ## D^{(1)} is first differences, so it correponds to l=0
    ## Dl = gen_diff_mat(n = TT, l = l, x = x)  <---  (T-l) x T matrix 
  }

  ## Running the fused LASSO
  ## which solves min_zhat 1/2 |z-zhat|_2^2 + lambda |D^{(1)}zhat|
  fit <- prox_dp(z = xi, lam = mod_lam) 

  return(fit)
}

Z_update  <- function(m, Uz, C, rho){
  mat = m + Uz/rho
  Z = projCmat(mat, C)
  return(Z)
}

#' Project rows of matrix \code{mat} into a ball of size \code{C}.
#' 
#' @param mat Matrix whose rows will be projected into a C-sized ball.
#' @param C radius
#'
#' @return Projected matrix.
projCmat <- function(mat, C){
  if(!is.null(C)){
    vlens = sqrt(rowSums(mat * mat))
    inds = which(vlens > C)
    if(length(inds) > 0){
      mat[inds,] = mat[inds,] * C / vlens[inds]
    }
  }
  return(mat)
}
