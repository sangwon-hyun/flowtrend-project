# Generated from _main.Rmd: do not edit by hand

#' Makes three 2-dimensional plots of the 3d data
#' 
#' @param ylist Data. A list of (|nt| by |dimdat|) matrices
#' @param countslist Count data.
#' @param obj flowmix or flowtrend object
#' @param tt Time point
#' @param return_list_of_plots If TRUE, return the list of three plots instead
#'   of the combined plot
#' @param zero_one_list_censored A list of zeros and 1s marking the censored
#'   points.
#' 
#' @export
#' @return
#'   
plot_3d <- function(ylist, obj = NULL, tt, countslist = NULL,
                    mn_colours = NULL,
                    labels = NULL,
                    bin = TRUE,
                    plot_title = NULL,
                    return_list_of_plots = FALSE,
                    zero_one_list_censored = NULL){

  ## Basic checks
  if(!is.null(obj)) stopifnot(class(obj) %in% c("flowmix", "flowtrend"))
  stopifnot(ncol(ylist[[1]]) == 3)
  if(!is.null(labels)) assertthat::assert_that(length(labels) == obj$numclust)
##  if(!bin) stop("This function is only for binned 3d data!")

  ## Extract data
  ## y = ylist[[tt]][,dims]
  labs = colnames(ylist[[1]])
  if(is.null(labs))  labs = paste0("dim", 1:3)
  if(is.null(countslist)){
    counts = rep(1, nrow(ylist[[tt]]))
  }
  if(!is.null(countslist)){
    counts = countslist[[tt]]
  }


  ## Aggregate counts into the two dimensions
  y2d_list = list()
  counts2d_list = list()
  for(ii in 1:3){
    dims = list(c(1:2), c(2:3), c(3,1))[[ii]]
    if(bin){
      yy = flowmix::collapse_3d_to_2d(ylist[[tt]], counts, dims)
      y2d = yy[,1:2]
      colnames(y2d) = labs[dims]
      one_counts = yy[,3]
    } else {
      y2d = ylist[[tt]][,dims]
      colnames(y2d) = labs[dims]
      one_counts = counts
    }
    y2d_list[[ii]] = y2d
    counts2d_list[[ii]] = one_counts
  }
  total_range = sapply(counts2d_list, range) %>% range()


  if(bin) zero_one = NULL
  if(!bin & !is.null(zero_one_list_censored)) zero_one = list(zero_one_list_censored[[tt]])
  if(!bin & is.null(zero_one_list_censored)) zero_one = NULL 

  ## Create three 2d plots
  plotlist = list()
  for(ii in 1:3){
    dims = list(c(1:2), c(2:3), c(3,1))[[ii]]

    ## Make data plot 
    ##one_countslist = (if(!is.null(countslist)) list(counts2d_list[[ii]]) else NULL)
    p = plot_2d(list(y2d_list[[ii]]),
                           list(counts2d_list[[ii]]), obj = NULL, tt = 1, bin = bin,
                           zero_one_list_censored = zero_one)
    if(bin){
      p$scales$scales <- list()
      p = p + scale_fill_gradientn(
                  colors = c("white", "blue", "yellow"), limits = total_range,
                  guide = "none")
    }

    ## Adding visualizations of the model |obj| at time tt
    if(!is.null(obj)){
      if(is.null(mn_colours)){
        mn_colours = rep("red", obj$numclust)
      }
      
      ## Make data matrix
      mnlist = lapply(1:obj$numclust, function(iclust){
        one_mnmat = obj$mn[tt,dims,iclust] %>% t()
        colnames(one_mnmat) = paste0("dim", 1:2)
        one_mnmat %>% as_tibble() %>% add_column(cluster = iclust)
      })
      mnmat = do.call(rbind, mnlist)

      for(iclust in 1:obj$numclust){
        
        ## Add ellipse
        el = ellipse::ellipse(x = obj$sigma[iclust,dims,dims],
                              centre = obj$mn[tt,dims,iclust]) %>% as_tibble()
        p = p + geom_path(aes(x = x, y = y), data = el,
                          colour = mn_colours[iclust], lty = 2,
                          linewidth = pmin(obj$prob[tt,iclust] * 8, 0.8))
        
        ## Add mean
        p = p + geom_point(aes(x = dim1, y = dim2),
                           data = mnmat %>% subset(cluster == iclust),
                           colour = mn_colours[iclust],
                           size = obj$prob[tt,iclust] * 10)

        ## Add cluster number as a label
        cex = rel(3)
        fac = 10
        if(is.null(labels)){
          labels = 1:(obj$numclust)
        }
        dt = data.frame(dim1 = mnmat[,1], dim2 = mnmat[,2], prob = obj$prob[tt,] * fac)
        p = p + ggrepel::geom_text_repel(aes(x = dim1, y = dim2, label = labels,
                                             point.size = sqrt(prob)),
                                         col = mn_colours,
                                         cex = cex,
                                         bg.color = "white",
                                         bg.r = 0.1,
                                         fontface = "bold",
                                         force_pull   = 5, # do not pull toward data points
                                         ## data = mnmat,
                                         data = dt,
                                         seed = 1)

      }
    }

    ## Format a bit more and save
    if(ii == 1){
      if(is.null(plot_title)) plot_title = paste0("Time=", tt)
      p = p + ggtitle(plot_title)
    } else {
      p = p + ggtitle("")
    }
    p = p + theme(legend.position = "none") 
    plotlist[[ii]] = p
  }

  p_combined = my_mfrow(plotlist)

  ## Return the plots
  if(return_list_of_plots) return(plotlist)
  return(p_combined)
}
