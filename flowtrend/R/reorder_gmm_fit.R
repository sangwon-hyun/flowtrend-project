# Generated from _main.Rmd: do not edit by hand

#' Reordering function for gmm results.
#'
#' @param obj result of \code{gmm_overfit()} or \code{gmm_underfit()}.
#' @param new_order new ordering of clusters to use.
#'
#' @return reordered object.
#' @export
reorder_gmm_fit <- function(obj, new_order){

  ## Setup
  numclust = max(new_order)

  ## Reorder everything
  obj$prob = obj$prob[,new_order]
  obj$mu = obj$mu[,new_order]
  obj$sigma = obj$sigma[,new_order]
  obj$memlist =  lapply(obj$memlist, function(a){
    a_copy = a
    for(iclust in 1:numclust){
      a_copy[which(a == iclust)] = new_order[iclust]
    }
    return(a_copy)
  })

  ## Reorder the responsibilities
  obj$resp_list =  lapply(obj$resp_list, function(oneresp){
    return(oneresp[,new_order])
  })

  return(obj)
}
