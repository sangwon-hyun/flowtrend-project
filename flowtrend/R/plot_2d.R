# Generated from _main.Rmd: do not edit by hand

#' Simple plotter for 2d particle data.
#' 
#' @param ylist Data. A list of (|nt| by |dimdat|) matrices
#' @param countslist Count data.
#' @param obj flowtrend (or flowmix) object.
#' @param time Out of 1 through \code{lengthy(list)}, which time point to plot.
#' @param zero_one_list_censored A list of zeros and 1s marking the censored
#'   points.
#' 
#' @export
#' @return ggplot object.
plot_2d <- function(ylist, countslist = NULL, obj = NULL, tt, bin = TRUE,
                    point_color = "blue", raster_colours = c("white", "blue"),
                    zero_one_list_censored = NULL){

  ## Basic checks
  if(!is.null(obj)) stopifnot(class(obj) %in% c("flowmix", "flowtrend"))
  stopifnot(ncol(ylist[[1]]) == 2)
  if(!is.null(obj)) stopifnot(obj$dimdat == 2)
  ##  if(!bin) stop("2d plotting for particle level data isn't supported.")
  

  ## Take data from one time point
  y = ylist %>% .[[tt]]
  if(is.null(colnames(y))){ colnames(y) = paste0("dim", c(1,2)) }
  y = y %>% as_tibble()


  ## Handle counts
  if(is.null(countslist)){
    counts = rep(1, nrow(ylist[[1]]))
  }
  if(!is.null(countslist)){
    counts = countslist[[tt]]
  }

  ## Get variable names
  varnames = y %>% colnames()
  varname1 = varnames[1]
  varname2 = varnames[2]

  ## Get data from one timepoint
  y = y %>% add_column(counts = counts)

  ## If a list of 0's and 1's (marking censored points) is provided, plot it
  if(!is.null(zero_one_list_censored)){
    zero_one = zero_one_list_censored[[tt]]
    y = y %>% add_column(censored = factor(zero_one, levels = c(FALSE, TRUE)))##c(0,1)))
  }

  if(!bin){
    p = y %>% ggplot() +
      geom_point(aes(x = !!sym(varname1), y=!!sym(varname2)), size = rel(counts),
                 alpha = .2, col = point_color) +
      theme_minimal() +
      theme(legend.position = "none")  + 
      ## scale_size()
      scale_size_area()
    if(!is.null(zero_one_list_censored)){
      p =
        y %>% ggplot() + geom_point(aes(x = !!sym(varname1), 
                                        y = !!sym(varname2), col = zero_one), size = rel(counts), alpha = 0.2) +
        theme_minimal() + theme(legend.position = "none") + 
        ## scale_size() +
        scale_size_area() +
        scale_color_manual(values=c("blue", "greenyellow"))
    }

  }
  if(bin){
    p =
      y %>% ggplot() +
      geom_raster(aes(x = !!sym(varname1), y=!!sym(varname2), fill = counts))  +
      scale_fill_gradientn(guide="none", colours = raster_colours) 
  }

  p = p + ggtitle(paste0("Time=", tt))

  ## Adding visualizations of the model |obj|
  if(is.null(obj)){
    return(p)
  } else {
    mnlist = lapply(1:obj$numclust, function(iclust){
      one_mnmat = obj$mn[,,iclust]
      colnames(one_mnmat) = paste0("dim", 1:2)
      one_mnmat %>% as_tibble() %>% add_column(cluster = iclust)
    })
    mnmat = do.call(rbind, mnlist)

    mn_colours = rep("red", obj$numclust)
    for(iclust in 1:obj$numclust){

      ## Add ellipse
      el = ellipse::ellipse(x = obj$sigma[iclust,,],
                            centre = obj$mn[tt,,iclust]) %>% as_tibble()
      p = p + geom_path(aes(x = x, y = y), data = el,
                        colour = mn_colours[iclust], lty = 2,
                        linewidth = pmin(obj$prob[tt,iclust] * 8, 0.8))

      ## Add mean
      p = p + geom_point(aes(x = dim1, y = dim2),
                         data = mnmat %>% subset(cluster == iclust) %>% .[tt,],
                         colour = mn_colours[iclust],
                         ## size = rel(3))
                         size = obj$prob[tt,iclust] * 10)
    }
  }
  return(p)
}




## This is from the many-cruises 02-helpers.R file; come back to this.
if(FALSE){
  p = datobj_2d %>%
    ggplot() +
    theme_minimal() +
    geom_raster(aes(x = !!sym(varname1), y=!!sym(varname2), fill = counts)) +
    scale_fill_gradientn(colours = colours, guide="colorbar")+
    xlim(c(0,8)) + ylim(c(0, 8)) +
    theme(legend.position = "none")
}
