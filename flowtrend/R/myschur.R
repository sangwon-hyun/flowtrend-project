# Generated from _main.Rmd: do not edit by hand

#' @param mat Matrix to Schur-decompose.
myschur <- function(mat){
  stopifnot(nrow(mat) == ncol(mat))
  if(is.numeric(mat) & length(mat)==1) mat = mat %>% as.matrix()
  obj = Matrix::Schur(mat)
  obj$tQ = t(obj$Q)
  obj$orig = mat
  return(obj)
}
