# Generated from _main.Rmd: do not edit by hand

#' Initialize the cluster centers.
#'
#' @param ylist  A T-length list of (nt  by 3) datasets.  There should  be T of
#'   such datasets. 3 is actually \code{mulen}.
#' @param numclust Number of clusters (M).
#' @param TT total number of (training) time points.
#'
#' @return An array of dimension (T x dimdat x M).
#' @export
init_mn <- function(ylist, numclust, TT, dimdat, countslist = NULL, seed=NULL){

  if(!is.null(seed)){
    assertthat::assert_that(all((seed %>% sapply(., class)) == "integer"))
    assertthat::assert_that(length(seed) == 7)
    RNGkind("L'Ecuyer-CMRG") 
    .Random.seed <<- seed
  }

  if(is.null(countslist)){
    ntlist = sapply(ylist, nrow)
    countslist = lapply(ntlist, FUN = function(nt) rep(1, nt))
  }

  ## Initialize the means by (1) collapsing to one cytogram (2) random
  ## sampling from this distribution, after truncation,
  TT = length(ylist)
  ylist_downsampled <- lapply(1:TT, function(tt){

    y = ylist[[tt]]
    counts = countslist[[tt]]

    ## Sample so that, in total, we get mean(nt) * 30 sized sample. In the case
    ## of binned data, nt is the number of bins.
    if(nrow(y) > 500) nsize = pmin(nrow(y) / TT * 30, nrow(y)) else nsize = nrow(y)
    some_rows = sample(1:nrow(y),
                       size = nsize,
                       prob = counts/sum(counts))
    y[some_rows,, drop=FALSE]
  })

  ## Jitter the means a bit
  yy = do.call(rbind, ylist_downsampled)
  new_means = yy[sample(1:nrow(yy), numclust),, drop=FALSE]
  jitter_sd = apply(yy, 2, sd) / 100
  jitter_means = MASS::mvrnorm(n = nrow(new_means),
                               mu = rep(0, dimdat),
                               Sigma = diag(jitter_sd, ncol = dimdat))
  new_means = new_means + jitter_means

  ## Repeat TT times == flat/constant initial means across time.
  mulist = lapply(1:TT, function(tt){ new_means })

  ## New (T x dimdat x numclust) array is created.
  muarray = array(NA, dim=c(TT, dimdat, numclust))
  for(tt in 1:TT){
    muarray[tt,,] = as.matrix(mulist[[tt]])
  }
  return(muarray)
}


#' Initialize the covariances.
#'
#' @param data The (nt by 3) datasets. There should be T of them.
#' @param numclust Number of clusters.
#' @param fac Value to use for the diagonal of the (dimdat x dimdat) covariance
#'   matrix.
#'
#' @return An (K x dimdat x dimdat) array containing the (dimdat by dimdat)
#'   covariances.
#' @export
init_sigma <- function(data, numclust, fac = 1){

  ndat = nrow(data[[1]])
  pdat = ncol(data[[1]])
  sigmas = lapply(1:numclust, function(iclust){
    onesigma = diag(fac * rep(1, pdat))
    if(pdat==1) onesigma = as.matrix(fac)
    colnames(onesigma) = paste0("datcol", 1:pdat)
    rownames(onesigma) = paste0("datcol", 1:pdat)
    return(onesigma)
  })
  sigmas = abind::abind(sigmas, along=0)
  return(sigmas)
}
