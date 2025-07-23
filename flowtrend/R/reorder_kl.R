# Generated from _main.Rmd: do not edit by hand

#' Reorder the cluster numbers for a new flowtrend object \code{newres}; the best
#' permutation (reordering) is to match the original flowtrend object
#' \code{origres}. Also works for flowmix objects.
#'
#' @param newres New flowtrend object to reorder.
#' @param origres Original flowtrend object.
#' @param ylist_particle The particle-level data.
#' @param fac Defaults to 100, to take 1/100'th of the particles from each time point.
#' @param verbose Loud or not?
#'
#' @return Reordered flowtrend object.
#'
#' @export
reorder_kl <- function(newres, origres, ylist_particle, fac = 100, verbose = FALSE){

  ## Randomly sample  1/100 of the original particles (mainly for memory reasons)
  TT = length(ylist_particle)
  N = sapply(ylist_particle, nrow) %>% sum()
  ntlist = sapply(ylist_particle, nrow)
  indlist = lapply(1:TT, function(tt){
    nt = ntlist[[tt]]
    ind = sample(1:nt, round(nt / fac), replace=FALSE)
  })

  ## Sample responsibilities
  ylist_particle_small = Map(function(ind, y){ y[ind,,drop = FALSE]  }, indlist, ylist_particle)

  ## Calculate new responsibilities
  resp_orig_small <- Estep(origres$mn, origres$sigma, origres$prob,
                           ylist = ylist_particle_small,
                           numclust = origres$numclust, first_iter = TRUE)
  resp_new_small <- Estep(newres$mn, newres$sigma, newres$prob,
                          ylist = ylist_particle_small,
                          numclust = newres$numclust, first_iter = TRUE)
  assertthat::assert_that(all(sapply(resp_orig_small, dim) == sapply(resp_new_small, dim)))

  ## Get best ordering (using symm. KL divergence and Hungarian algorithm for
  ## matching)
  best_ord <- get_best_match_from_kl(resp_new_small, resp_orig_small)

  if(verbose) cat("New order is", best_ord, fill=TRUE)
  newres_reordered_kl = newres %>% reorder_clust(ord = best_ord)

  ## Return the reordered object
  return(newres_reordered_kl)
}
