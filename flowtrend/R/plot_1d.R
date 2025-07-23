# Generated from _main.Rmd: do not edit by hand

#' Makes 1d plot of data and model
#'
#' @param ylist Data. A list of (|nt| by |dimdat|) matrices
#' @param obj A flowtrend (or flowmix) object. Defaults to NULL.
#' @param x Time points. Defaults to NULL.
#' @param alpha Between 0 and 1, how transparent to plot the data
#'   points. Defaults to 0.1.
#' @param bin If TRUE, the data is binned.
#'
#' @return ggplot object with data, and optionally, a flowtrend model overlaid.
#' @export
plot_1d <- function(ylist, countslist=NULL, obj=NULL, x = NULL, alpha = .1, bin = FALSE,
                    plot_band = TRUE){

  ## Basic checks
  if(!is.null(obj))stopifnot(class(obj) %in% c("flowmix", "flowtrend"))
  if(!is.null(x)){
    stopifnot(length(x) == length(ylist))
    times = x
  } else {
    times = 1:length(ylist)
  }
  dimdat = ncol(ylist[[1]])
  assertthat::assert_that(dimdat == 1)

  ## If countslist is not provided, make a dummy
  if(is.null(countslist)) countslist = lapply(ylist, function(y) rep(1, nrow(y)))


  ## Make data into long matrix
  ymat <- lapply(1:length(ylist), FUN = function(tt){
    data.frame(time = times[tt], Y = ylist[[tt]], counts = countslist[[tt]])
  }) %>% bind_rows() %>% as_tibble()
  colnames(ymat) = c("time", "Y", "counts") ## when ylist[[tt]] already has a column name, this is needed.
  
  if(bin){
    gg = ymat %>% ggplot() +
      geom_raster(aes(x = time, y = Y, fill = counts)) +
      theme_bw() + ylab("Data") + xlab("Time")  +
      scale_fill_gradientn(colours = c("white", "black"))
  }
  if(!bin){
    gg = ymat %>% ggplot() +
      geom_point(aes(x = time, y = Y), alpha = alpha) +
      theme_bw() + ylab("Data") + xlab("Time") 
  }

  ## If there is a flowtrend object to plot, do it.
  if(is.null(obj)){
    return(gg)
  } else {
  
    ## Add the model
    numclust = obj$numclust
    mnmat = obj$mn %>% .[,1,] %>%  `colnames<-`(1:numclust) %>% as_tibble() %>% 
      add_column(time = times)
    probmat = obj$prob %>% as_tibble() %>% setNames(1:numclust) %>% add_column(time = times)
    mn_long = mnmat %>% pivot_longer(-time, names_to = "cluster", values_to = "mean")
    prob_long = probmat %>% pivot_longer(-time, names_to = "cluster", values_to = "prob")
    est_long = full_join(mn_long, prob_long, by = c("time","cluster"))
    gg = gg + geom_path(aes(x = time, y = mean, linewidth = prob, group = cluster, color = cluster),
                        data = est_long,
                        lineend = "round", linejoin="mitre") +
      scale_linewidth(range = c(0.05, 5), limits = c(0, 1))

    if(plot_band){
      ## Add the estimated 95% probability regions for data.
      stdev = obj$sigma %>% .[,,1] %>% sqrt()
      mn_long_by_clust = mn_long %>% group_by(cluster) %>% group_split()
      band_long_by_clust = lapply(1:numclust, function(iclust){
        mn_long_by_clust[[iclust]] %>% mutate(upper = mean + 1.96 * stdev[iclust]) %>% mutate(lower = mean - 1.96 * stdev[iclust])
      })
      band_long = band_long_by_clust %>% bind_rows()
      gg = gg + geom_line(aes(x = time, y = upper, group = cluster, color = cluster),
                     data = band_long, size = rel(.7), alpha = .5) +
        geom_line(aes(x = time, y = lower, group = cluster, color = cluster),
                  data = band_long, size = rel(.7), alpha = .5) +
      guides(size = "none") # To turn off line size from legend
    }
    return(gg)
  }
}
