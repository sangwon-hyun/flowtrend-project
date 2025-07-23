# Generated from _main.Rmd: do not edit by hand

#' A "soft" version of a rand index between two sets of responsibility
#' (membership probability) matrices. Measures the /disagreement/ between the two clusterings.
#' 
#' @param resp_list1 One list of responsibility matrices.
#' @param resp_list2 Another list of responsibility matrices.
#' @param times Optional; if you would like to isolate your attention to some specific times.
#'
#' @return A single soft rand index number
#' @export
rand_old <- function(resp_list1, resp_list2, times = NULL, prop = .25){

  if(!is.null(times)) resp_list1 = resp_list1[times]
  if(!is.null(times)) resp_list2 = resp_list2[times]

  rand_onetime <- function(resp1, resp2){
    stopifnot(nrow(resp1) == nrow(resp2))
    stopifnot(ncol(resp1) == ncol(resp2))
    nt = nrow(resp1)
    
    ## Form the score
    mat11 = resp1 %*% t(resp1)
    mat22 = resp2 %*% t(resp2)
    mat11not = (1-mat11)
    mat22not = (1-mat22)

    ## Make the diagonals not matter anywhere
    diag(mat11) = 0
    diag(mat22) = 0
    diag(mat11not) = 0
    diag(mat22not) = 0

    a = mat11 * mat22
    b = mat11not * mat22not
    c = mat11 * mat22not
    d = mat11not * mat22

    return(c(a = sum(a), b = sum(b), c = sum(c), d = sum(d)))
  }
  abcd_over_time = mapply(rand_onetime, resp_list1, resp_list2)
  abcdmat = abcd_over_time %>% t()
  abcd = abcdmat %>% colSums()
  score = (abcd["a"] + abcd["b"])/  (sum(abcd))
  return(score)
}

#' Rand index between responsibilities.
#'
#' @param resp_list1 One list.
#' @param resp_list2 Another list.
#' @param times A subset of the time points (out of 1:length(resp_list1)) to
#'   examine.
#' @param smaller If TRUE, use only a small (sampled) subset of the particles'
#'   responsibilities for calculating RAND.
#' @param prop How much to downsample; defaults to 0.1.
#' @export
rand <- function(resp_list1, resp_list2, times = NULL, smaller = TRUE, prop = 0.1){
 
  ## Basic checks
  stopifnot(all(sapply(resp_list1, nrow) == sapply(resp_list2, nrow)))

  ## Subset the times if needed
  if(!is.null(times)) resp_list1 = resp_list1[times]
  if(!is.null(times)) resp_list2 = resp_list2[times]

  ## names(everything)
  ## list2env(everything,envir = environment())
  ## resp_list1 = resp_oracle
  ## resp_list2 = true_resp
  if(smaller){
    indslist = lapply(resp_list1, function(oneresp){
      inds = sample(x=1:nrow(oneresp), size=ceiling(nrow(oneresp)*prop), replace=FALSE) %>% sort()
    })
    resp_list1 = mapply(function(a,b)a[b,,drop=FALSE], resp_list1, indslist, SIMPLIFY=FALSE)
    resp_list2 = mapply(function(a,b)a[b,,drop=FALSE], resp_list2, indslist, SIMPLIFY=FALSE)
  }

  ## Make each list into one long matrix
  resp1 = Reduce(rbind, resp_list1)
  resp2 = Reduce(rbind, resp_list2)
  ## resp1 = do.call(rbind, resp_list1)## %>% bind_rows()
  ## resp2 = do.call(rbind, resp_list2)## %>% bind_rows()
  
  stopifnot(nrow(resp1) == nrow(resp2))
  stopifnot(ncol(resp1) == ncol(resp2))
  nt = nrow(resp1)
  
  ## Form the score
  mat11 = resp1 %*% t(resp1)
  mat22 = resp2 %*% t(resp2)
  mat11not = (1-mat11)
  mat22not = (1-mat22)

  ## Make the diagonals not matter anywhere
  diag(mat11) = 0
  diag(mat22) = 0
  diag(mat11not) = 0
  diag(mat22not) = 0

  a = mat11 * mat22
  b = mat11not * mat22not
  c = mat11 * mat22not
  d = mat11not * mat22

  abcd = c(a = sum(a), b = sum(b), c = sum(c), d = sum(d))
  score = (abcd["a"] + abcd["b"])/  (sum(abcd))
  return(score)
}



#' 2 x 2 contigency table from membership vectors.
#' 
#' @param mem1 Membership vector.
#' @param mem2 Another membership vector.
#'
#' @return A 3x3 matrix containing (1) a 2 x 2 table in the first [1:2,1:2]
#'   entries, and (2) row sums and column sums and total sums in the [,3], [3,]
#'   entries.
#' @export
make_contingency_table<-function(mem1, mem2){
  tab = table(mem1, mem2)
  tab = cbind(tab, rowSums(tab))
  tab = rbind(tab, colSums(tab))
  return(tab)
}


#' Calculate RAND index from a contingency table.
#' 
#' @param tab  A matrix containing  a 2 x 2 table in the first [1:2,1:2]
#'   entries; e.g., from make_contingency_table().
#'
#' @return Rand index.
#' @export
get_rand_from_table <- function(tab){
  numer = sum(sapply(as.numeric(tab), function(nij) choose(nij, 2))) +
    tab[1,1] * tab[2,2] + tab[1,2] * tab[2,1]
  denom = (choose(sum(tab), 2))
  ri = numer/denom
  return(ri)
}


#' Faster rand index calculation from membership vectors.
#'
#' @param mem1
#' @param mem2 
#' 
#' @return RAND index.
#' @export
rand_from_mems <- function(mem1, mem2){
  make_contingency_table(mem1, mem2) %>%
    .[1:2,1:2] %>%
    get_rand_from_table()
}
