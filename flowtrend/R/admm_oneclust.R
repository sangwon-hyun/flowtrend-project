# Generated from _main.Rmd: do not edit by hand

#' One cluster's ADMM for a fixed step size (\code{rho}).
#' 
#' 
#' @export
#' @param iclust 
#' @param niter 
#' @param y 
#' @param Dl 
#' @param tDl 
#' @param Dlp1 
#' @param l Default: NULL
#' @param TT 
#' @param N 
#' @param dimdat 
#' @param maxdev 
#' @param rho 
#' @param rhoinit Default: rho
#' @param Xinv 
#' @param schurA 
#' @param schurB 
#' @param sigmainv 
#' @param lambda 
#' @param resp 
#' @param resp_sum 
#' @param ylist 
#' @param err_rel Default: 1e-3
#' @param err_abs Default: 0
#' @param zerothresh 
#' @param mu 
#' @param z 
#' @param w 
#' @param uw 
#' @param uz 
#' @param first_iter 
#' @param em_iter 
#' @param outer_iter 
#' @param local_adapt 
#' @param sigma 
#' @param sigma_eig_by_clust 
#'
#' @return A list containing the updated parameters for the cluster.
#' @seealso \code{\link{la_admm_oneclust}} which is a wrapper over this function.
admm_oneclust <- function(iclust = 1, niter, y,
                          Dl, tDl, Dlp1, l = NULL,
                          TT, N, dimdat, maxdev,
                          rho,
                          rhoinit = rho,
                          Xinv,
                          schurA,
                          schurB,
                          sigmainv,
                          lambda,
                          resp,
                          resp_sum,
                          ylist, err_rel = 1e-3, err_abs = 0,
                          zerothresh,
                          mu, 
                          z,
                          w,
                          uw,
                          uz,
                          first_iter,## Not used 
                          em_iter,
                          outer_iter,
                          local_adapt,
                          sigma,
                          sigma_eig_by_clust){

  ## Initialize the variables ###
  ## resid_mat = matrix(NA, nrow = ceiling(niter/5), ncol = 4)
  ## colnames(resid_mat) = c("primresid", "primerr", "dualresid", "dualerr")
  resid_mat = matrix(NA, nrow = ceiling(niter/5), ncol = 6)
  colnames(resid_mat) = c("prim1", "prim2", "primresid", "primerr", "dualresid", "dualerr")
  rhofac = rho / rhoinit 

  ## This doesn't change over iterations
  schurB = myschur(schurB$orig * rhofac) ## In flowmix, this is done on A. Here, it's done on B (in AX + XB + C = 0).
  TA = schurA$T 
  TB = schurB$T
  UA = schurA$Q
  UB = schurB$Q
  tUA = schurA$tQ
  tUB = schurB$tQ

  ## This also doesn't change over iterations
  C1 <- do.call(cbind, lapply(1:TT, FUN = function(tt){
      multmat <- apply(y[[tt]], FUN = function(yy) yy * resp[[tt]], MARGIN = 2)
      sigmainv %*% colSums(multmat)
    }))

  
  for(iter in 1:niter){
    syl_C <- get_C_mat(C1 = C1, resp_sum = resp_sum, TT = TT, dimdat = dimdat,
                       Sigma_inv = sigmainv, N = N, Dl = Dl, rho = rho,
                       z = z, w = w, uz = uz, uw = uw, l=l)
    FF =  (-1) * tUA %*% syl_C %*% UB
    mu = UA %*% matrix_function_solve_triangular_sylvester_barebonesC2(TA, TB, FF) %*% tUB
    mu = t(mu)
    stopifnot(nrow(mu) == TT)
    stopifnot(ncol(mu) == dimdat)

    centered_mu = sweep(mu, 2, colMeans(mu)) 
    z <- Z_update(centered_mu, Uz = uz, C = maxdev, rho = rho)

    if(any(abs(mu)>1E2)){
      stop("mu is blowing up!")
      ## break
    }
    wlist = lapply(1:dimdat, function(j){
      W_update_fused(l = l, TT = TT, mu = mu[, j, drop = TRUE],
                     rho = rho, lambda = lambda,
                     uw = uw[,j,drop=TRUE],
                     Dl = Dl)})
    w <- do.call(cbind, wlist)
    stopifnot(nrow(w) == TT-l)
    stopifnot(ncol(w) == dimdat)


    uz = U_update_Z(uz, rho, mu, z, TT)


    ## uw = U_update_W(uw, rho, mu, w, l, TT)
    uw = U_update_W(uw, rho, mu, w, l, Dl, TT)


    ## Check convergence
    if( iter > 1  & iter %% 5 == 0){## & !local_adapt){

      ## Calculate convergence criterion
      obj = check_converge(mu, rho,
                           w, z,
                           w_prev, z_prev,
                           uw, uz,
                           Dl, tDl, err_rel = err_rel, err_abs = err_abs)

      jj = (iter/ 5)
      resid_mat[jj,] = c(
          norm(obj$primal_resid1, "F"), ## Temp
          norm(obj$primal_resid2, "F"), ## Temp
          norm(obj$primal_resid, "F"),
                         obj$primal_err,
                         norm(obj$dual_resid,"F"),
                         obj$dual_err)

      if(is.na(obj$converge)){
        obj$converge = converge = FALSE
        warning("Convergence was NA")
      }
      if(obj$converge){
        converge = TRUE
        break
      } else {
        converge = FALSE
      }
    }
    w_prev = w
    z_prev = z
  }

  if(FALSE){
    ## Calculate optimization objective values for this cluster.
    obj.value <- objective_per_cluster(y = y, mu = mu, resp = resp,
                                       Sigma_inv = sigmainv, TT = TT,
                                       d = dimdat, Dlp1 = Dlp1, Dl = Dl, l = l,
                                       maxdev = maxdev, lambda = lambda,
                                       rho = rho, N = N)
    ## This is very expensive to do within each ADMM iteration, so it's commented out for now.
  }
  obj.value = NA

  return(list(mu = mu,
              resid_mat = resid_mat,
              converge = converge,
              ## Other variables to return.
              Z = z,
              W = w,
              uz = uz,
              uw = uw,
              inner.iter = iter,
              single_admm_objective = obj.value))
}
