# Generated from _main.Rmd: do not edit by hand

#' The M step for the cluster probabilities
#'
#' @param resp Responsibilities.
#' @param H_tf Trend filtering matrix.
#' @param countslist Particle multiplicities.
#' @param lambda_prob Regularization. 
#' @param l_prob Trend filtering degree.
#'
#' @return (T x k) matrix containing the alphas, for \code{prob = exp(alpha)/
#'         rowSums(exp(alpha))}. 
#' @export
#'
Mstep_prob <- function(resp, H_tf, countslist = NULL,
                       lambda_prob = NULL, l_prob = NULL, x = NULL){

  ## Basic setup
  TT <- length(resp)

  ## Basic checks
  stopifnot(is.null(l_prob) == is.null(lambda_prob))

  ## If glmnet isn't actually needed, don't use it.
  if(is.null(l_prob) & is.null(lambda_prob)){

    ## Calculate the average responsibilities, per time point.
    if(is.null(countslist)){
      resp.avg <- lapply(resp, colMeans) %>% do.call(rbind, .)
    } else {
      resp.avg <- lapply(1:TT, FUN = function(ii){
        colSums(resp[[ii]])/sum(countslist[[ii]])
      }) %>% do.call(rbind, .)
    }
    return(resp.avg) 

  ## If glmnet is needed, use it.
  } else {

    lambda_range <- function(lam, nlam = 50, lam.max = max(1, 5*lam)){
      return(exp(seq(log(lam.max), log(lam), length.out = nlam)))
    }

    penalty.facs <- c(rep(0, l_prob+1), rep(1, nrow(H_tf) - l_prob - 1))
    resp.predict <- do.call(rbind, lapply(resp, colSums))
    fac = resp.predict %>% apply(2, function(mycol){max(mycol)/min(mycol)}) %>% max()

    if(FALSE){
    ## Run GLMNET but if it throws an error, reduce the dynamic range |fac| of
    ## each column |resp.predict| until it's okay.
    error_means_null = NULL
    fac = 1024
    while(is.null(error_means_null)){
      print("fac is")
      print(fac)
      error_means_null <- tryCatch({
        resp.predict = resp.predict0 %>% apply(2, function(mycol){mycol = pmax(mycol, max(mycol)/fac)})
        glmnet_obj <- glmnet::glmnet(x = H_tf, y = resp.predict, family = "multinomial",
                                     penalty.factor = penalty.facs, maxit = 1e7,
                                     lambda =  mean(penalty.facs)*lambda_range(lambda_prob),
                                     standardize = F, intercept = FALSE) 
        pred_link <- predict(glmnet_obj, newx = H_tf, type = "link", s = mean(penalty.facs) * lambda_prob)[,,1]
      }, error = function(e){
        return(NULL)
      })
      if(is.null(error_means_null)){
        fac = fac/2
        resp.predict = resp.predict0 %>% apply(2, function(mycol){mycol = pmax(mycol, max(mycol)/fac)})
      }
    }
    } else {
        glmnet_obj <- glmnet::glmnet(x = H_tf, y = resp.predict, family = "multinomial",
                                     penalty.factor = penalty.facs, maxit = 1e7,
                                     lambda =  mean(penalty.facs)*lambda_range(lambda_prob),
                                     standardize = F, intercept = FALSE) 
    }
    pred_link <- predict(glmnet_obj, newx = H_tf, type = "link", s = mean(penalty.facs) * lambda_prob)[,,1]
    return(pred_link)
  } 
}
