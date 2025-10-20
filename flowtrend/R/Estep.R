# Generated from _main.Rmd: do not edit by hand

#' E step, which updates the "responsibilities", which are posterior membership probabilities of each particle.
#'
#' @param mn
#' @param sigma covariance
#' @param prob
#' @param ylist
#' @param numclust
#' @param denslist_by_clust
#' @param first_iter
#' @param countslist
#' @param padding A small amount of padding to add to the weighted
#'   densities. Note, a very small value (like 1E-20) should be used, otherwise
#'   the probabilitistic (soft) gating can be affected quite noticably.
#'
#' @return
#' @export
#'
Estep <- function(mn, sigma, prob, ylist = NULL, numclust, denslist_by_clust = NULL,
                  first_iter = FALSE, countslist = NULL, padding = 1E-20){
  ## Basic setup
  TT = length(ylist)
  ntlist = sapply(ylist, nrow)
  resp = list()
  dimdat = dim(mn)[2]
  assertthat::assert_that(dim(mn)[1] == length(ylist))

  ## Helper to calculate Gaussian density for each \code{N(y_{t,k},mu_{t,k} and
  ## Sigma_k)}.
  calculate_dens <- function(iclust, tt, y, mn, sigma, denslist_by_clust,
                             first_iter) {
    mu <- mn[tt, , iclust]
    if (dimdat == 1) {
      dens = dnorm(y, mu, sd = sqrt(sigma[iclust, , ]))
    } else {
       dens = dmvnorm_arma_fast(y, mu, sigma[iclust,,], FALSE)
    }
    return(dens)
  }
  ## Calculate posterior probability of membership of $y_{it}$.
  ncol.prob = ncol(prob)
  empty_row = rbind(1:numclust)[-1,]
  for (tt in 1:TT) {
    ylist_tt = ylist[[tt]]
    if(nrow(ylist_tt) == 0){
      resp[[tt]] <- empty_row
    } else {
      densmat <- sapply(1:numclust, calculate_dens, tt, ylist_tt,
                        mn, sigma, denslist_by_clust, first_iter) 
      wt.densmat <- matrix(prob[tt, ], nrow = ntlist[tt],
                           ncol = ncol.prob, byrow = TRUE) * densmat
      wt.densmat = wt.densmat + padding##1e-20
      wt.densmat <- wt.densmat/rowSums(wt.densmat)
      resp[[tt]] <- wt.densmat
    }
  }

  ## Weight the responsibilities by $C_{it}$.
  if (!is.null(countslist)) {
    resp <- Map(function(myresp, mycount) {
      myresp * mycount
    }, resp, countslist)
  }
  return(resp)
}
