# Generated from _main.Rmd: do not edit by hand

#' Soft-gates a responsibility matrix by a bernoulli (or multinoulli) draw.
#' 
#' @param oneresp A 2-column responsibility matrix
#' 
#' @return
#' @export
soft_gate_one_responsibility_matrix <- function(oneresp){

  ## Setup
  numclust = ncol(oneresp)
  vec = rep(0,numclust)

  ## Draw the 0-1 memberships
  zero_one_mat = apply(oneresp, 1, function(myrow){
    draw = sample(1:numclust, size=1, prob=myrow)
    vec[draw]= 1
    vec
  }) %>% t()

  ## Check dimensions and return
  stopifnot(all(dim(zero_one_mat) == dim(oneresp)))
  return(zero_one_mat)
}
