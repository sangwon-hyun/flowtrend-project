# Generated from _main.Rmd: do not edit by hand

#' Sequentially permute a time series of GMMs.
#'
#' @param gmm_list Output of an lapply of gmm_each over ylist.
#' @param numclust Number of clusters.
#'
#' @return Same format as |gmm_list|
match_clusters_gmm <- function(gmm_list, numclust = 2){

  ## Initialize some objects
  TT <- length(gmm_list)
  
  ## Helper to permute one GMM
  my_reorder <- function(one_gmm, new_order){
    one_gmm_reordered = one_gmm
    one_gmm_reordered$tab = one_gmm_reordered$tab %>%
      mutate(hard_cluster = hard_cluster %>%
               plyr::revalue(replace=c("1"=new_order[1],"2"=new_order[2]))) %>%
      mutate(soft_cluster = soft_cluster %>%
               plyr::revalue(replace=c("1"=new_order[1],"2"=new_order[2])))

    ## reorder the parameters
    param = one_gmm$param
    new_mns = param[,c("mn1", "mn2")][new_order] %>% as.numeric()
    new_sds = param[,c("sd1", "sd2")][new_order] %>% as.numeric()
    new_probs = param[,c("prob1", "prob2")][new_order] %>% as.numeric()
    one_gmm_reordered$param <- data.frame(mn1 = new_mns[1], mn2 = new_mns[2],
                                          sd1 = new_sds[1], sd2 = new_sds[2],
                                          prob1 = new_probs[1], prob2 = new_probs[2])
    ## if(all(new_order == 2:1)) browser()

    ## reorder the responsibilities
    one_gmm_reordered$resp = one_gmm$resp[,new_order]
    return(one_gmm_reordered)
  }
  
  ## Fill in the rest of the rows by sequentially matching clusters from t-1 to
  ## t, using the Hungarian Algorithm
  new_orders = matrix(NA, ncol = 2, nrow = TT)
  new_orders[1,] = c(1, 2)
  gmm_list_copy = gmm_list
  for(tt in 2:TT){

    ## Find order
    param1 = gmm_list_copy[[tt-1]] %>% .$param
    param2 = gmm_list[[tt]] %>% .$param
    klmat <- symmetric_kl_between_gaussians(param1 = param1, param2 = param2)
    new_order <- RcppHungarian::HungarianSolver(klmat) %>%
      .$pairs %>% data.frame() %>% arrange(.[,2]) %>% .[,1]
    new_orders[tt,] = new_order

    ## Perform the reordering
    gmm_list_copy[[tt]] <- my_reorder(gmm_list[[tt]], new_order)
  }
  return(gmm_list_copy)
}

#' Between two sets of Gaussian mean/sd parameters, what is the symmetric KL
#' distance.
#' @param param1 list of mn1, mn2, sd1, sd2
#' @param param2 another such list.
#'
#' @return (numclust x numclust) distance matrix.
symmetric_kl_between_gaussians <- function(param1, param2, numclust = 2){

  ## First model's parameters
  mn_model1 = c(param1$mn1, param1$mn2)
  sd_model1 = c(param1$sd1, param1$sd2)

  ## Second model's parameters
  mn_model2 = c(param2$mn1, param2$mn2)
  sd_model2 = c(param2$sd1, param2$sd2)

  ## Fill out distance matrix
  dist_cs <- matrix(0, ncol = numclust, nrow = numclust)
  for(iclust_1 in 1:numclust){
    for(iclust_2 in 1:numclust){
      dist_cs[iclust_2, iclust_1] <-
        one_symmetric_kl(mu1 = mn_model1[iclust_1], mu2 = mn_model2[iclust_2],
                         sd1 = sd_model1[iclust_1], sd2 = sd_model2[iclust_2])
    }
  }
  rownames(dist_cs) <- rep("C1", numclust)
  colnames(dist_cs) <- rep("C2", numclust)
  return(dist_cs)
}

#' Symmetric KL divergence between two Gaussian distributions.
#' @param mu1 Mean for distribution 1.
#' @param mu2 Mean for distribution 2.
#' @param sd1 Standard deviation for distribution 1.
#' @param sd2 Standard deviation for distribution 2.
one_symmetric_kl <- function(mu1, mu2, sd1, sd2){
  stopifnot(length(mu1)==1 & length(mu2) == 1 &
            length(sd1) == 1 & length(sd2) == 1)
  kl12 <- log(sd2/sd1) + (sd1^2 + (mu1 - mu2)^2)/(2*sd2^2) - 1/2
  kl21 <- log(sd1/sd2) + (sd2^2 + (mu2 - mu1)^2)/(2*sd1^2) - 1/2
  kl_sym <- 0.5*(kl12 + kl21)
  return(kl_sym)
}


#' Apply a gmm in each 1-dimensional cytogram.
#' 
#' @param one_y One-column matrix.
#' @param numclust Number of clusters.
#' @return 
#' @export
gmm_each <- function(one_y, numclust){
  
  ## Basic check
  stopifnot(ncol(one_y) == 1)

  ## Do Gaussian Mixture model
  ## obj <- mclust::Mclust(data = one_y, G = numclust,
##                             modelNames = "V", verbose = FALSE)
  obj <- my_mclust(one_y, numclust, FALSE)

  ## Memberships (soft- and hard-clustered)
  hard_mem = obj$z  %>% apply(1, which.max)
  soft_mem = obj$z  %>%
    apply(1, function(myrow){
      sample(1:2, size = 1, replace = FALSE, prob=myrow)
    })
  tab <- tibble(y = as.numeric(one_y),
                soft_cluster = factor(soft_mem, levels = c(1,2)),
                hard_cluster = factor(hard_mem, levels = c(1,2)))

  ## parameter table
  mn1 = obj$parameters$mean[[1]]
  mn2 = obj$parameters$mean[[2]]
  sds = obj %>% .$parameters %>% .$variance %>% .$sigmasq %>% sqrt()
  sd1 = sds[1]
  sd2 = sds[2]
  prob1 = obj$parameters$pro[1]
  prob2 = obj$parameters$pro[2]

  ## ggplot(data.frame(x=as.numeric(one_y))) +
  ## geom_histogram(aes(x=x)) +
  ## geom_vline(xintercept = mn1, col = 'blue') +
  ## geom_vline(xintercept = mn1 + 2*sd1, col = 'skyblue') +
  ## geom_vline(xintercept = mn1 - 2*sd1, col = 'skyblue') +
  ## geom_vline(xintercept = mn2, col = 'red') +
  ## geom_vline(xintercept = mn2 - 2*sd2, col = 'orange') +
  ## geom_vline(xintercept = mn2 + 2*sd2, col = 'orange')
  param = tibble(mn1 = mn1,
                 mn2 = mn2,
                 sd1 = sd1,
                 sd2 = sd2,
                 prob1 = prob1,
                 prob2 = prob2)
  resp = obj$z

  return(list(tab = tab, param = param, resp = resp))##, sigma_list = sigma_list))
}
