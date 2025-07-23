# Generated from _main.Rmd: do not edit by hand

#' Pool the entire series of cytograms, fit a GMM model, and re-aggregate parameters.
#' 
#' @param ylist Data.
#' @param numclust Number of clusters (only works for 2).
#' @return 
#' @export
overfit_gmm <- function(ylist, numclust = 2, reorder = TRUE){ 

  ## Basic checks
  stopifnot(numclust == 2)

  ## Estimate individual GMM models
  gmm_list <- lapply(ylist, gmm_each, numclust = numclust)

  ## Reorder the cluster memberships sequentially
  if(reorder){
    gmm_list = match_clusters_gmm(gmm_list, numclust = 2)
  }

  ## Obtain the resulting data in long format
  tab_list <- lapply(gmm_list, function(a) a$tab)
  TT = length(ylist)
  names(tab_list) = 1:TT
  tab_long = tab_list %>% bind_rows(.id = "time") %>%
    pivot_longer(-c("time", "y"), values_to = "cluster",
                                   names_to = "type")
  memlist = tab_long %>% subset(type == "soft_cluster") %>% group_by(time) %>%
    group_split() %>% lapply(function(a)pull(a, cluster))

  ## Get the Gaussian distribution parameters
  param_list <- lapply(gmm_list, function(a) a$param)
  param_mat = param_list %>% bind_rows(.id = "time")
  mu = param_mat[,c("mn1", "mn2")] %>% as.matrix()
  prob = param_mat[,c("prob1", "prob2")] %>% as.matrix()
  sigma = param_mat[,c("sd1", "sd2")] %>% .^2 %>% as.matrix()

  ## Responsibilities
  resp_list = gmm_list %>% lapply(function(a) a$resp)

  ## Make a list of each time points' sigmas
  TT = nrow(sigma)
  sigma_list = lapply(1:TT, function(tt){
    one_row = sigma[tt,] %>% as.numeric()
    one_sigma = array(NA, dim = c(2,1,1))
    one_sigma[,1,1] = one_row
    return(one_sigma)
  })

  return(list(tab_long = tab_long,
              memlist = memlist,
              param_mat = param_mat,
              mu = mu,
              prob = prob,
              sigma = sigma,
              sigma_list = sigma_list,
              resp_list = resp_list,
              numclust = numclust))
}
