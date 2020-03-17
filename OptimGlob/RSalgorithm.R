### Random search Algorithm  #######

random_search <- function(test_fun, param) {
  # recover params
  LB <- param$LB 
  UB <- param$UB  
  budget <- param$budget 
  dim <- param$dim
  #
  best_par <- NA
  best_value <- Inf
  x_hist <- matrix(, budget, dim)
  y_hist <- matrix(, budget, 1)
  for (i in 1:budget) {
    par <- runif(dim, min=LB, max=UB)
    value <- test_fun(par)
    x_hist[i, ] <- par
    y_hist[i, ] <- value
    if (value < best_value) {
      best_par <- par
      best_value <- value
    }      
  }
  res <- list(xhist=x_hist, fhist=y_hist, x_best=best_par, f_best=best_value)
  return(res)
}
