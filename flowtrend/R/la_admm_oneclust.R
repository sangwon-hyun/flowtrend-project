# Generated from _main.Rmd: do not edit by hand

#' Locally adaptive ADMM; wrapper to \code{admm_oneclust()}.
#'
#' @param K
#' @param ...
#'
#' @return
#'
#' @examples
la_admm_oneclust <- function(K, ...){

  ## Initialize arguments for ADMM.
  args <- list(...)
  p = args$p
  l = args$l
  TT = args$TT
  dimdat = args$dimdat
  rhoinit = args$rhoinit

  ## This initialization can come from the previous *EM* iteration.
  if(args$first_iter){
    mu = matrix(0, nrow = TT, ncol = dimdat)
    Z <- matrix(0, nrow = TT, ncol = dimdat)
    W <- matrix(0, nrow = TT-l, ncol = dimdat)
    uz <- matrix(0, nrow = TT, ncol = dimdat)
    uw <- matrix(0, nrow = TT - l, ncol = dimdat)

    args[['mu']] <- mu
    args[['z']] <- Z
    args[['w']] <- W
    args[['uz']] <- uz
    args[['uw']] <- uw
  }

  cols = c()
  some_admm_objectives = c()


  ## Run ADMM repeatedly with (1) double rho, and (2) previous b
  for(kk in 1:K){
    if(kk > 1){

      ## These ensure warm starts are true
      args[['mu']] <- mu
      args[['z']] <- Z
      args[['w']] <- W
      args[['uz']] <- uz
      args[['uw']] <- uw
      args[['rho']] <- rho
    }

    ## Run ADMM
    args[['outer_iter']] <- kk

    ## Call main function
    argn <- lapply(names(args), as.name)
    names(argn) <- names(args)
    call <- as.call(c(list(as.name("admm_oneclust")), argn))
    res = eval(call, args)

    if(any(abs(res$mu)>1E2)){
      stop("mu is blowing up! Probably because the initial ADMM step size (rho) is too large (and possibly the ball constraint on the means is large.") 
    }

    some_admm_objectives = c(some_admm_objectives, res$single_admm_objective) 

    ## Handling the scenario where the objectives are all zero
    padding = 1E-12
    some_admm_objectives = some_admm_objectives + padding

    ## See if outer iterations should terminate
    if(res$converge){
      res$converge <- T
      break
    }

    ## Update some parameters; double the rho value, and update the B matrix
    rho = 2 * args$rho
    mu = res$mu
    Z = res$Z
    W = res$W
    uz = res$uz
    uw = res$uw
  }
  if(!res$converge) {
    ## warning("ADMM didn't converge for one cluster.") ## Disabling this because it's not really something to fix. The only thing change a user would make is to make |rho_init| smaller 
  }

  ## Record how long the admm took; in terms of # iterations.
  res$kk = kk

  ## Record the final rho
  res$rho = args$rho

  return(res)
}
