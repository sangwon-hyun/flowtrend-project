# Generated from _main.Rmd: do not edit by hand

#' Makes cluster probability plot (lines over time).
#'
#' @param obj Estimated model (from e.g. \code{flowtrend()})
#'
#' @export
plot_prob <- function(obj, x = NULL){

  ## Basic checks
  if(!is.null(obj)) stopifnot(class(obj) %in% c("flowmix", "flowtrend"))
  if(!is.null(x)){
    times = x
  } else {
    ##stop("must provide x")
    times = 1:(obj$TT)
  }

  numclust = obj$numclust
  probmat = obj$prob %>% as_tibble() %>% setNames(1:numclust) %>% add_column(time = times)
  prob_long = probmat %>% pivot_longer(-time, names_to = "cluster", values_to = "prob")
  prob_long %>% ggplot() +
    geom_line(aes(x=time, y = prob, group = cluster, col = cluster), size = rel(1)) +
    ggtitle("Estimated cluster probability")
}
