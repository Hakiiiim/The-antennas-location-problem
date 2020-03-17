### Normal search Algorithm followed by LBFGS #######
# param is just the concatenation of those for normal_search and lbfgs
NS_lbfgs <- function(test_fun, param) {
  res_NS <- normal_search(ofwrapper, param)
  param_lbfgs <- param
  param_lbfgs$xinit <- res_NS$x_best # start lbfgs from best of normal search
  res_lbfgs <- lbfgs(ofwrapper, param_lbfgs)
  # concatenate both results
  res <- list() 
  res$xhist <- rbind(res_NS$xhist,res_lbfgs$xhist)
  res$fhist <- rbind(res_NS$fhist,res_lbfgs$fhist)
  res$x_best <- res_lbfgs$x_best
  res$f_best <- res_lbfgs$f_best
  # res <- list(xhist=x_hist, fhist=y_hist, x_best=best_par, f_best=best_value)
    # param <- list(LB=-5,UB = 5,maxit = 12,dim=zdim, xinit=rep(-4,zdim), trace=0, lmm =1, factr=100,pgtol=1.e-2) # param for lbfgs
  return(res)
}
