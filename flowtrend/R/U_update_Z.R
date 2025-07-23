# Generated from _main.Rmd: do not edit by hand

#' @param U (T x dimdat) matrix.
U_update_Z <- function(U, rho, mu, Z, TT){
 # return(U + rho * (scale(mu, scale = F) - Z))
  stopifnot(nrow(U) == TT)

  centered_mu = sweep(mu, 2, colMeans(mu)) 
  ## stopifnot(all(abs(colMeans(centered_mu))<1E-8))
  Unew = U + rho * (centered_mu - Z)

  ## Expect a (T-l) x dimdat matrix.
  stopifnot(all(dim(U) == dim(Unew))) 
  stopifnot(nrow(U) == TT)
  return(Unew)
}
