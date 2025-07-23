# Generated from _main.Rmd: do not edit by hand

#' Helper to aggregate CV results and obtain the |res| object, all saved in
#' |destin|.
#'
#' @inheritParams cv_aggregate
#'
#' @return List containing, for every (iprob, imu), the "best" estimated model
#'   out of the |nrestart| replicates (best in the sense that it had the best
#'   likelihood value out of the |nrestart| replicates.)
cv_aggregate_res <- function(destin){

  load(file.path(destin, "meta.Rdata"))

  res.list = list()
  for(iprob in 1:cv_gridsize){
    for(imu in 1:cv_gridsize){ 

      ## If the "best" object has already been created, use it.
      best_refit_filename = make_best_refit_filename(iprob, imu)
      if(file.exists(file.path(destin, best_refit_filename))){
        load(file.path(destin, best_refit_filename), verbose = FALSE)
        bestres = res

      ## Otherwise, aggregate directly from the individual files
      } else {
        obj = rep(NA, nrestart)
        res.list.inner = list()
        for(irestart in 1:nrestart){
          filename = make_refit_filename(iprob, imu, irestart)
          tryCatch({
            load(file.path(destin, filename))
            res.list.inner[[irestart]] = res
            obj[irestart] = res$objectives[length(res$objectives)]
          }, error = function(e){ NULL })
        }
        if(!all(is.na(obj))){
          bestres = res.list.inner[[which.min(obj)]] ## which.min?
        }
      }

      ## Add the "best" object to a list.
      res.list[[paste0(iprob, "-", imu)]] = bestres
    }
  }
  return(res.list)
}
