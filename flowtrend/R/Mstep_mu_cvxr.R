# Generated from _main.Rmd: do not edit by hand

#' Used as a fallback or to test against \code{Mstep_mu()}.
#' 
#' @param ylist
#' @param resp
#' @param lambda
#' @param l
#' @param Sigma_inv inverse of Sigma
#' @param x covariates
Mstep_mu_cvxr <- function(ylist,
                          resp,
                          lambda,
                          l,
                          Sigma_inv, 
                          x = NULL,
                          thresh = 1E-8,
                          maxdev = NULL,
                          dimdat,
                          N,
                          ecos_thresh = 1E-8,
                          scs_eps = 1E-5){
  
  ## Define dimensions
  TT = length(ylist)
  
  ## Responsibility Weighted Data
  ytildes <- lapply(1:TT, FUN = function(tt){
    yy <- ylist[[tt]]
    g <- resp[[tt]]
    yy <- apply(yy, MARGIN = 2, FUN = function(x) x * g)
    yrow = matrix(c(NA, NA), nrow=1, ncol=dimdat)
    yrow[1,] = colSums(yy)
    yrow
  }) %>% do.call(rbind, .)
  
  ## Auxiliary term, needed to make the objective interpretable
  aux.y <- Reduce("+", lapply(1:TT, FUN = function(tt){
    yy <- ylist[[tt]]
    g <- sqrt(resp[[tt]])
    yy <- apply(yy, MARGIN = 2, FUN = function(x) x * g)
    sum(diag(yy %*% Sigma_inv %*% t(yy)))
  }))
  
  ## Mu, d x T matrix
  mumat <- CVXR::Variable(cols=dimdat, rows=TT)
  
  ## Summed sqrt responsibilities - needed in the objective.
  resp.sum.sqrt <- lapply(resp, FUN = function(x) sqrt(sum(x)))
  
  ## Differencing Matrix, (TT-(l+1)) x TT 
  Dlp1 <- gen_diff_mat(n = TT, l = l+1, x = x)
  # l = 2 is quadratic trend filtering
  # l = 1 is linear trend filtering
  # l = 0 is fused lasso
  
  ## Forming the objective
  obj = 1/(2*N) *( Reduce("+", lapply(1:TT, FUN = function(tt) CVXR::quad_form(t(resp.sum.sqrt[[tt]]*mumat[tt,]), Sigma_inv))) -2 * Reduce("+", lapply(1:TT, FUN = function(tt) t(ytildes[tt,]) %*% Sigma_inv %*% t(mumat[tt,]))) + aux.y) + lambda * sum(CVXR::sum_entries(abs(Dlp1 %*% mumat), axis = 1))

  ## Putting together the ball constraint
  rowmns <- matrix(rep(1, TT^2), nrow = TT)/TT
  mu_dotdot <- rowmns %*% mumat
  constraints = list()
  if(!is.null(maxdev)){
    constraints = list(CVXR::sum_entries(CVXR::square(mumat - mu_dotdot), axis = 2) <= rep(maxdev^2, TT) )
  }
  
  ## Try all two CVXR solvers.
  prob <- CVXR::Problem(CVXR::Minimize(obj), constraints)
  result = NULL
  result <- tryCatch({
    CVXR::solve(prob, solver="ECOS",
          FEASTOL = ecos_thresh, RELTOL = ecos_thresh, ABSTOL = ecos_thresh)
  }, error=function(err){
    err$message = paste(err$message, 
                        "\n", "Lasso solver using ECOS has failed." ,sep="")
    cat(err$message, fill=TRUE)
    return(NULL)
  })
  
  ## If anything is wrong, flag to use SCS solver
  scs = FALSE
  if(is.null(result)){
    scs = TRUE
  } else {
    if(result$status != "optimal") scs = TRUE
  }
  
  ## Use the SCS solver
  if(scs){
    result = CVXR::solve(prob, solver="SCS", eps = scs_eps)
    if(any(is.na(result$getValue(mumat)))){ ## A clumsy way to check.
      stop("Lasso solver using both ECOS and SCS has failed.", sep="")
    }
  }
  
  ## Record Interesting Parameters
  num_iters <- result$num_iters
  status <- result$status
  mumat <- result$getValue(mumat)
  val <- result$value
  
  return(list(mu = mumat, value = val, status = status, num_iters = num_iters))
}
