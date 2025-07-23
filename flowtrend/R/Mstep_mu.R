# Generated from _main.Rmd: do not edit by hand

#' Computes the M step for mu.
#' 
#' @param resp Responsibilities of each particle.
#' @param ylist 
#' @param lambda
#' @param l
#' @param sigma
#' @param sigma_eig_by_clust
#' @param Dl
#' @param Dlp1
#' @param TT
#' @param N
#' @param dimdat
#' @param first_iter
#' @param mus
#' @param Zs
#' @param Ws
#' @param uws
#' @param uzs
#' @param maxdev
#' @param x
#' @param niter
#' @param err_rel
#' @param err_abs
#' @param zerothresh
#' @param local_adapt
#' @param local_adapt_niter
#' @param rho_init
#' @param iter
#'
#' @return
#' @export
#'
#' @examples
Mstep_mu <- function(resp,
                     ylist,  
                     lambda = 0.5,
                     l = 3,
                     sigma,
                     sigma_eig_by_clust = NULL,
                     Dlsqrd,
                     Dl, tDl, Dlp1,  TT, N, dimdat,
                     first_iter = TRUE,
                     e_mat,

                     ## Warm startable variables
                     mus = NULL,
                     Zs = NULL,
                     Ws = NULL,
                     uws = NULL,
                     uzs = NULL,
                     ## End of warm startable variables

                     maxdev = NULL,
                     x = NULL,
                     niter = (if(local_adapt) 1e2 else 1e3),
                     err_rel = 1E-3,
                     err_abs = 0,
                     zerothresh = 1E-6,
                     local_adapt = FALSE,
                     local_adapt_niter = 10,
                     rho_init = 0.01,
                     iter = NULL){

  ####################
  ## Preliminaries ###
  ####################
  TT = length(ylist)
  numclust = ncol(resp[[1]])
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)
  resp.sum = lapply(resp, colSums) %>% do.call(rbind, .)
  N = sum(unlist(resp.sum))

  ## Other preliminaries
  schur_syl_A_by_clust = schur_syl_B_by_clust = term3list = list()
  ybarlist = list()
  ycentered_list = Xcentered_list = yXcentered_list = list()
  Qlist = list()
  sigmainv_list = list()
  for(iclust in 1:numclust){

    ## Retrieve sigma inverse from pre-computed SVD, if necessary
    if(is.null(sigma_eig_by_clust)){
      sigmainv = solve(sigma[iclust,,])
    } else {
      sigmainv = sigma_eig_by_clust[[iclust]]$sigma_inv
    }

    resp.iclust <- lapply(resp, FUN = function(r) matrix(r[,iclust]))

    AB <- get_AB_mats(y = y, resp = resp.iclust, Sigma_inv = sigmainv,
                      e_mat = e_mat, N = N, Dlp1 = Dlp1, Dl = Dl,
                      Dlsqrd = Dlsqrd, rho = rho_init, z = NULL, w = NULL,
                      uz = NULL, uw = NULL)

    ## Store the Schur decomposition
    schur_syl_A_by_clust[[iclust]] = myschur(AB$A)
    schur_syl_B_by_clust[[iclust]] = myschur(AB$B)

    ycentered <- NULL
    ycentered_list[[iclust]] = ycentered
    sigmainv_list[[iclust]] = sigmainv
  }

  ##########################################
  ## Run ADMM separately on each cluster ##
  #########################################
  admm_niters = admm_inner_iters = vector(length = numclust, mode = "list")
  if(first_iter){ mus = vector(length = numclust, mode = "list") }
  Zs <- lapply(1:numclust, function(x) matrix(0, nrow = TT, ncol = dimdat))
  Ws <- lapply(1:numclust, function(x) matrix(0, nrow = TT - l, ncol = dimdat))
  uzs <- lapply(1:numclust, function(x) matrix(0, nrow = TT, ncol = dimdat))
  uws <- lapply(1:numclust, function(x) matrix(0, nrow = TT - l, ncol = dimdat))

  ## For every cluster, run LA-ADMM
  start.time = Sys.time()
  for(iclust in 1:numclust){

    resp.iclust <- lapply(resp, FUN = function(r) matrix(r[,iclust]))

    ## Possibly locally adaptive ADMM, for now just running with rho == lambda
    res = la_admm_oneclust(K = (if(local_adapt) local_adapt_niter else 1),
                           local_adapt = local_adapt,
                           iclust = iclust,
                           niter = niter,
                           TT = TT, N = N, dimdat = dimdat, maxdev = maxdev,
                           schurA = schur_syl_A_by_clust[[iclust]],
                           schurB = schur_syl_B_by_clust[[iclust]],
                           sigmainv = sigmainv_list[[iclust]],
                           rho = rho_init, 
                           rhoinit = rho_init, 
                           ## rho = (if(iclust == 1) rho_init else res$rho/2),
                           ## rhoinit = (if(iclust == 1) rho_init else res$rho/2),
                           sigma = sigma,
                           lambda = lambda,
                           resp = resp.iclust,
                           resp_sum = resp.sum[,iclust],
                           l = l,
                           Dlp1 = Dlp1,
                           Dl = Dl,
                           tDl = tDl,
                           y = ylist,
                           err_rel = err_rel,
                           err_abs = err_abs,
                           zerothresh = zerothresh,
                           sigma_eig_by_clust = sigma_eig_by_clust,
                           em_iter = iter,

                           ## Warm starts from previous *EM* iteration
                           first_iter = first_iter,
                           uw = uws[[iclust]],
                           uz = uzs[[iclust]],
                           z = Zs[[iclust]],
                           w = Ws[[iclust]])


    ## Store the results
    mus[[iclust]] = res$mu
    admm_niters[[iclust]] = res$kk
    admm_inner_iters[[iclust]] = res$inner.iter

    ## Store other things for for warmstart
    Zs[[iclust]] = res$Z
    uzs[[iclust]] = res$uz
    uws[[iclust]] = res$uw
    Ws[[iclust]] = res$W
  }

  ## Aggregate the yhats into one array
  mu_array = array(NA, dim = c(TT, dimdat, numclust))
  for(iclust in 1:numclust){ mu_array[,,iclust] = mus[[iclust]] }

  ## Each are lists of length |numclust|.
  return(list(mns = mu_array,
              admm_niters = admm_niters, 
              admm_inner_iters = admm_inner_iters,

              ## For warmstarts
              Zs = Zs,
              Ws = Ws,
              uws = uws,
              uzs = uzs,
              N = N,

              ## Return the column
              rho = res$rho,

              ## For using in the Sigma M step
              ycentered_list = ycentered_list,
              Xcentered_list = Xcentered_list,
              yXcentered_list = yXcentered_list,
              Qlist = Qlist
  ))
}
