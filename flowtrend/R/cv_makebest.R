# Generated from _main.Rmd: do not edit by hand

#' From the results saved in \code{destin}, aggregate all |nrestart| files to retain only the "best" restart, and delete the rest.
#' 
#' All meta information (|nfold|, |cv_gridsize|, |nrestart|, |lambda_means|, |lambda_probs|) comes from \code{meta.Rdata}.
#'
#' @param destin Directory with cross-validation output.
#'
#' @export
cv_makebest <- function(destin){  

  ## Read the meta data (for |nfold|, |cv_gridsize|, |nrestart|, |lambda_means|, |lambda_probs|)
  load(file = file.path(destin, 'meta.Rdata'), verbose = FALSE)
  
  ## This loads all the necessary things; just double-checking.
  stopifnot(exists("nrestart"))
  stopifnot(exists("nfold"))
  stopifnot(exists("cv_gridsize"))
  stopifnot(exists(c("lambda_probs")))
  stopifnot(exists(c("lambda_means")))

  ## Aggregate the results
  cvscore.array = array(NA, dim = c(cv_gridsize, cv_gridsize, nfold, nrestart))
  cvscore.mat = matrix(NA, nrow = cv_gridsize, ncol = cv_gridsize)
  for(iprob in 1:cv_gridsize){
    for(imu in 1:cv_gridsize){
      for(ifold in 1:nfold){
        print(c(iprob, imu, ifold))

        ## If the "best" flowtrend object has already been created, do nothing.
        best_filename = make_best_cvscore_filename(iprob, imu, ifold)
        if(file.exists(file.path(destin, best_filename))){
          next

        ## Otherwise, attempt to load from all |nrestart| replicates
        } else {
          objectives = load_all_objectives(destin, iprob, imu, ifold, nrestart)

          ## If all |nrestart| files exist, delete all files but the best model.
          if(all(!is.na(objectives))){
            best_irestart = which(objectives == min(objectives)) %>% .[1] ## If there is a tie, leave it.
            keep_only_best(destin, iprob, imu, ifold, nrestart, best_irestart)
          }
          else {
            print(paste0("iprob=", iprob, " imu=", imu, " ifold=", ifold, " had objectives:  ", objectives))
          }
        }
      }
    }
  }

  ## Also go over the "refit" files
  for(iprob in 1:cv_gridsize){
    for(imu in 1:cv_gridsize){

      ## If the "best" flowtrend object has already been created, do nothing.
      best_filename = make_best_refit_filename(iprob, imu)
      if(file.exists(file.path(destin, best_filename))){

      ##   ## Check if any more jobs have been done since before
      ##   objectives = load_all_refit_objectives(destin, iprob, imu, nrestart)
      ##   if(any(!is.na(objectives))){
      ##     nonmissing_irestart = which(!is.na(objectives))## == min(objectives))
      ##     ## keep_only_best_refit(destin, iprob, imu, nrestart, best_irestart)
      ##     load(file.path(destin, best_filename))
      ##   } else {
      ##     print(paste0("iprob=", iprob, " imu=", imu, " had /refit/ objectives:  ", objectives))
      ##   }
      ## }

        next
        
      ## Otherwise, attempt to load from all |nrestart| replicates
      } else {
        objectives = load_all_refit_objectives(destin, iprob, imu, nrestart)
        if(all(!is.na(objectives))){
          best_irestart = which(objectives == min(objectives)) %>% .[1] ## If there is a tie, leave it.
          keep_only_best_refit(destin, iprob, imu, nrestart, best_irestart)
        } else {
          print(paste0("iprob=", iprob, " imu=", imu, " had /refit/ objectives:  ", objectives))
        }
      }
    }
  }
}


#' Loading all objectives, with NA's for missing files
load_all_objectives <- function(destin, iprob, imu, ifold, nrestart){
  objectives = sapply(1:nrestart, function(irestart){
    filename = make_cvscore_filename(iprob, imu, ifold, irestart)
    tryCatch({
      load(file.path(destin, filename), verbose = FALSE)
      return(objectives[length(objectives)])
    }, error = function(e){ NA })
  })
  return(objectives)
}

#' Loading all objectives, with NA's for missing files
load_all_refit_objectives <- function(destin, iprob, imu, nrestart){
  objectives = sapply(1:nrestart, function(irestart){ 
    filename = make_refit_filename(iprob, imu, irestart)
    ## filename = make_best_cvscore_filename(iprob, imu, ifold)
    tryCatch({
      load(file.path(destin, filename), verbose = FALSE)
      return(res$objectives[length(res$objectives)])
    }, error = function(e){ NA })
  })
  return(objectives)
}

#' Keeping only the output files for the "best" restart, and deleting the rest.
keep_only_best <- function(destin, iprob, imu, ifold, nrestart, best_irestart){
  for(irestart in 1:nrestart){
    filename = make_cvscore_filename(iprob, imu, ifold, irestart)
    if(irestart == best_irestart){
      best_filename = make_best_cvscore_filename(iprob, imu, ifold)
      file.rename(from = file.path(destin, filename),
                    to = file.path(destin, best_filename))
    } else {
      file.remove(file.path(destin, filename))
    }
  }
}

#' Keeping only the output files (among refit models) for the "best" restart,
#' and deleting the rest.
keep_only_best_refit <- function(destin, iprob, imu, nrestart, best_irestart){ 
  for(irestart in 1:nrestart){
    filename = make_refit_filename(iprob, imu, irestart)
    if(irestart == best_irestart){
      best_filename = make_best_refit_filename(iprob, imu)
      file.rename(from = file.path(destin, filename),
                    to = file.path(destin, best_filename))
    } else {
      file.remove(file.path(destin, filename))
    }
  }
}
