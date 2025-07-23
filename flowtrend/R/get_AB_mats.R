# Generated from _main.Rmd: do not edit by hand

#' a
#'
#' @param y
get_AB_mats <- function(y, resp, Sigma_inv, e_mat, Dlsqrd, N, Dlp1, Dl, rho, z, w, uz, uw){

  # A matrix
  A <- 1/N * Sigma_inv

  # B matrix
  sum_resp <- sapply(resp, sum)
  #B <- rho*(t(Dl)%*%Dl + e_mat)%*% diag(1/unlist(sum_resp))
  B <- rho*(Dlsqrd + e_mat)
  B <- B/sum_resp[col(B)]
  #B <- rho*(Dlsqrd + e_mat)  %*% diag(1/unlist(sum_resp))

  return(list(A = A, B = B))
}
