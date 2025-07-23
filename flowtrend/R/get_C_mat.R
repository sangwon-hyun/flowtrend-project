# Generated from _main.Rmd: do not edit by hand

#' a
#' @param C1
get_C_mat <- function(C1, resp_sum, TT, dimdat, Sigma_inv, e_mat, N, Dlp1, Dl,
                      rho, z, w, uz, uw, l){

  C2 <- t(uz - rho*z)

  # averaging
  C2 <- C2 - rowMeans(C2)

  # third component
  C3 <- do.call(rbind, lapply(1:dimdat, FUN = function(j){
    ((uw[,j, drop=TRUE] - rho * w[,j,drop=TRUE]) %*% Dl) %>% as.numeric()
  }))

  # combining
  C <-  (-1/N * C1 + C2 + C3)
  C <- C/resp_sum[col(C)]
  return(C)
}
