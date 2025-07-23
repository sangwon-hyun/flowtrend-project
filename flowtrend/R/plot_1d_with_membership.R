# Generated from _main.Rmd: do not edit by hand

#' Makes 1d plot of data when
#'
#' @param ylist Particlel-level data. A list of (|nt| by |dimdat|) matrices.
#' @param memlist List of memberships for the particles.
#' @param x Time points. Defaults to NULL.
#' @param alpha Between 0 and 1, how transparent to plot the data
#'   points. Defaults to 0.1.
#'
#' @return ggplot object with data, and optionally, a flowtrend model overlaid.
#' @export
plot_1d_with_membership <- function(ylist, memlist, countslist = NULL, x = NULL, alpha = .01){

  ## Basic checks
  if(!is.null(x)){
    stopifnot(length(x) == length(ylist))
    times = x
  } else {
    times = 1:length(ylist)
  }
  dimdat = ncol(ylist[[1]])

  ## If countslist is not provided, make a dummy
  if(is.null(countslist)) countslist = lapply(ylist, function(y) rep(1, nrow(y)))

  ## Make data into long matrix
  ymat <- lapply(1:length(ylist), FUN = function(tt){
    data.frame(time = times[tt], Y = ylist[[tt]],
               counts = countslist[[tt]], cluster = memlist[[tt]])
  }) %>% bind_rows() %>% as_tibble()
  ymat = ymat %>% mutate(cluster = as.factor(cluster))
  colnames(ymat) = c("time", "Y", "counts", "cluster") ## when ylist[[tt]] already has a column name, this is needed.

  gg =
    ymat %>% ggplot() +
    geom_point(aes(x = time, y = Y, group=cluster, col = cluster),##, shape = cluster),
               alpha = alpha) + ##facet_wrap(~cluster)  +
    theme_bw() + ylab("Data") + xlab("Time") 
  return(gg)
}
