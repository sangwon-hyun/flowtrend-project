# Generated from _main.Rmd: do not edit by hand

#' Evaluating the penalized negative log-likelihood on all data \code{ylist} given parameters \code{mu}, \code{prob}, and \code{sigma}.
#'
#' @param mu
#' @param prob
#' @param prob_link
#' @param sigma
#' @param ylist
#' @param Dl
#' @param l
#' @param lambda
#' @param l_prob
#' @param Dl_prob
#' @param lambda_prob
#' @param alpha
#' @param beta
#' @param countslist
#' @param unpenalized if TRUE, return the unpenalized out-of-sample fit.
#'
#' @return
#' @export
#'
#' @examples
objective <- function(mu, prob, prob_link = NULL, sigma,
                      ylist,
                      Dlp1 = NULL, l = NULL,
                      lambda = 0,
                      l_prob = NULL,
                      Dlp1_prob = NULL,
                      lambda_prob = 0,
                      alpha = NULL, beta = NULL,
                      countslist = NULL,
                      unpenalized = FALSE){

  ## Set some important variables
  TT = dim(mu)[1]
  numclust = dim(mu)[3]
  if(is.null(countslist)){
    ntlist = sapply(ylist, nrow)
  } else {
    ntlist = sapply(countslist, sum)
  }
  N = sum(ntlist)
  dimdat = ncol(ylist[[1]])

  ## Calculate the log likelihood
  loglik = sapply(1:TT, function(tt){
    return(loglik_tt(ylist, tt, mu, sigma, prob, countslist, numclust = numclust, dimdat = dimdat))
  })

  if(unpenalized){
    obj =  -1/N * sum(unlist(loglik)) 
    return(obj)
  } else {

    ## Return penalized likelihood
    mu.splt <- asplit(mu, MARGIN = 2) ## This was 3, which I think produces the same result.
    diff_mu <- sum(unlist(lapply(mu.splt, FUN = function(m) sum(abs(Dlp1 %*% m)))))
    diff_prob <- sum(abs(Dlp1_prob %*% prob_link))
    obj =  -1/N * sum(unlist(loglik)) + lambda * diff_mu + lambda_prob * diff_prob
    return(obj)
  }
}

