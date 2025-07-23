# Generated from _main.Rmd: do not edit by hand  
testthat::test_that("Trend filtering regression matrix is created correctly on equally spaced data.",
{
  ## Check that equally spaced data creates same trendfilter regression matrix
  ## Degree 1
  H1 <- gen_tf_mat(10, 1)
  H1_other <- gen_tf_mat(10, 1, x=(1:10)/10)
  testthat::expect_equal(H1, H1_other)

  ## Degree 1
  H2 <- gen_tf_mat(10, 2)
  H2_other <- gen_tf_mat(10, 2, x=(1:10)/10)
  testthat::expect_equal(H2, H2_other)
  
  ## Check the dimension
  testthat::expect_equal(dim(H1), c(10,10))
  
  ## Check against an alternative function.
  for(ord in c(0,1,2,3,4)){
    H <- gen_tf_mat(10, ord)
    H_eq = gen_tf_mat_equalspace(10, ord)
    testthat::expect_true(max(abs(H_eq- H)) < 1E-10)
  }
})

testthat::test_that("Uneven spaced D matrix is formed correctly", {

  x = (1:6)[-(2)]
  ## k = 0 is piecewise constant, k=1 is piecewise linear, etc.
  for(k in 0:3){

    ## Temp
    k=1
    ## End of temp

    l = k+1

    ## Form Dl using our function
    Dl <- flowtrend::gen_diff_mat(n=5, l = l, x=x)
    
    ## Compare it to rows of H matrix, per section 6 of Tibshirani et al. 2014
    gen_tf_mat(n = length(x), k = k, x = x) %>%
      solve() %>%
      `*`(factorial(k)) %>% ## This part is missing in Tibshirani et al. 2014
      tail(length(x)-(k+1)) -> Hx 
    ratiomat = Hx/Dl

    ## Process the ratios of each entry, and check that they're equal to 1
    ratios = ratiomat[!is.nan(ratiomat)]
    ratios = ratios[is.finite(ratios)]
    testthat::expect_true(all.equal(ratios, rep(1, length(ratios)))) 
  }
})

testthat::test_that("Test for softmax",{
  link = runif(100, min = -10, max = 10) %>% matrix(nrow = 10, ncol = 10)
  testthat::expect_true(all(abs(rowSums(softmax(link)) - 1) < 1E-13))
})

testthat::test_that("E step returns appropriately sized responsibilities.",{
  
  ## Generate some fake data
  TT = 10
  dimdat = 1
  ylist = lapply(1:TT, function(tt){ runif(30*dimdat) %>% matrix(ncol = dimdat, nrow = 30)})
  numclust = 3

  ## Initialize a few parameters, not carefully
  sigma = init_sigma(ylist, numclust) ## (T x numclust x (dimdat x dimdat))
  mn = init_mn(ylist, numclust, TT, dimdat)##, countslist = countslist)
  prob = matrix(1/numclust, nrow = TT, ncol = numclust) ## Initialize to all 1/K.

  ## Calculate responsibility
  resp = Estep(mn = mn, sigma = sigma, prob = prob, ylist = ylist, numclust = numclust)

  ## Check these things
  testthat::expect_equal(length(resp), length(ylist))
  testthat::expect_equal(length(resp), length(ylist))
  testthat::expect_equal(sapply(resp, nrow), sapply(ylist, nrow))
  testthat::expect_true(all(sapply(resp, ncol) == numclust))
})

testthat::test_that("Mstep of pi returns a (T x K) matrix.", {

  ## Generate some fake responsibilities and trend filtering matrix
  TT = 100
  numclust = 3
  nt = 10
  resp = lapply(1:TT, function(tt){
    oneresp = runif(nt*numclust) %>% matrix(ncol=numclust)
    oneresp = oneresp/rowSums(oneresp)
  })
  H_tf <- gen_tf_mat(n = TT, k = 0)

  ## Check the size
  pred_link = Mstep_prob(resp, H_tf, l_prob = 0, lambda_prob = 1E-3)
  testthat::expect_equal(dim(pred_link), c(TT, numclust))
  pred_link = Mstep_prob(resp, H_tf)
  testthat::expect_equal(dim(pred_link), c(TT, numclust))

  ## Check the correctness
  pred_link = Mstep_prob(resp, H_tf)
})

testthat::test_that("Test the M step of \pi against CVXR", {})

testthat::test_that("Test the ball projection", {
set.seed(100)
mat = matrix(rnorm(100), ncol=2)
projected_mat = flowtrend:::projCmat(mat, 1)
testthat::expect_true(all(projected_mat %>% apply(1, function(myrow)sum(myrow*myrow)) < 1+1E-8))
})

testthat::test_that("The prediction function returns the right things", {
  ## Generate data
  set.seed(100)
  dt       <- gendat_1d(100, rep(100, 100))
  ylist = dt %>% dt2ylist()
  x = dt %>% pull(time) %>% unique()
  obj <- flowtrend(ylist = ylist,
                    x = x,
                    maxdev = 5,
                    numclust = 3,
                    lambda = 0.02,
                    l = 1,
                    l_prob = 2,
                    lambda_prob = .005, ## 
                    nrestart = 1,
                    niter = 3)
  predobj = predict_flowtrend(obj)
  testthat::expect_named(predobj, c("mn", "prob", "sigma", "x"))
})

## Generate data
set.seed(100)
dt       <- gendat_1d(100, rep(100, 100))
dt_model       <- gendat_1d(100, rep(100, 100), return_model = TRUE)
held_out = 25:35
dt_subset = dt %>% subset(time %ni% held_out)
ylist = dt_subset %>% dt2ylist()
x = dt_subset %>% pull(time) %>% unique()
set.seed(686)
obj <- flowtrend(ylist = ylist, 
                  x = x,
                  maxdev = 5,
                  numclust = 3,
                  l = 2,
                  l_prob = 2,
                  lambda = 0.02,
                  lambda_prob = .1, ## 
                 nrestart = 5,
                 verbose = TRUE)

obj$all_objectives %>% mutate(irestart = as.factor(irestart)) %>%
  ggplot() + geom_line(aes(x=iter, y=objective, group = irestart, col = irestart))

## Also reorder the cluster labels of the truth, to match the fitted model.
ord = obj$mn[,1,] %>% colSums() %>% order(decreasing=TRUE)
lookup <- setNames(c(1:obj$numclust), ord)
dt_model$cluster = lookup[as.numeric(dt_model$cluster)] %>% as.factor()

## Reorder the cluster labels of the fitted model.
obj = reorder_clust(obj)

testthat::test_that("prediction function returns the right things", {
   
  predobj = predict_flowtrend(obj, newtimes = held_out)

  ## Check a few things
  testthat::expect_equal(predobj$x,  held_out)
  testthat::expect_equal(rowSums(predobj$prob),  rep(1, length(held_out)))
  testthat::expect_equal(dim(predobj$mn),  c(length(held_out), 1, 3))
})

testthat::test_that("Objective value decreases over EM iterations.",{
  
  glist = list()
  for(iseed in 1:5){

    ## Generate synthetic data
    set.seed(iseed*100)
    dt       <- gendat_1d(100, rep(10, 100))
    dt_model       <- gendat_1d(100, rep(10, 100), return_model = TRUE)
    ylist = dt %>% dt2ylist()
    x = dt %>% pull(time) %>% unique()
    
    ## Fit model
    obj <- flowtrend_once(ylist = ylist,
                          x = x,
                          maxdev = 5,
                          numclust = 3,
                          lambda = 0.02,
                          l = 1,
                          l_prob = 2,
                          lambda_prob = 0.05)

    ## Test objective monotonicity
    niter_end = length(obj$objective)
    testthat::expect_true(all(diff(obj$objective) < 1E-3))
  

    ## Make a plot
    g = ggplot(tibble(iter=1:niter_end, objective=obj$objectives)) +
      geom_point(aes(x=iter, y=objective)) +
      geom_line(aes(x=iter, y=objective)) +
      ggtitle(paste0("Seed=", iseed*100)) +
      xlab("EM iteration")
    glist[[iseed]] = g
  }
  title = cowplot::ggdraw() + cowplot::draw_label("Objective values over EM iterations", fontface='bold')
  main_plot = cowplot::plot_grid(plotlist = glist, ncol=5, nrow=1) 
  cowplot::plot_grid(title, main_plot, ncol=1, rel_heights=c(0.1, 1)) %>% print()
})

