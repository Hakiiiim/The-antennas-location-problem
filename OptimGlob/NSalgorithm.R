### Normal search Algorithm  #######
### i.e. ES-(1+1) with constant sigma

normal_search <- function(test_fun, param) {
  # recover parameters
  budget <- param$budget
  xinit <- param$xinit
  sigma <- param$sigma
  dim <- length(xinit)  
  if ( dim != param$dim) stop("dim and xinit do not have same dimension")
  LB <- param$LB
  UB <- param$UB
  # expand LB and UB to the number of dimensions
  if (length(LB)==1 & dim>1) LB <- rep(LB[1],dim)
  if (length(UB)==1 & dim>1) UB <- rep(UB[1],dim)

  maxnbsamp <- 10000 # max nb of samples for bounds satisfaction
  
  
  current_par <- xinit
  current_value <- test_fun(xinit)
  best_par <- xinit
  best_value <- test_fun(xinit)
  x_hist <- matrix(, budget, dim)
  y_hist <- matrix(, budget, 1)
  x_hist[1, ] <- current_par
  y_hist[1] <- current_value

  for (i in 2:budget) {
    nb_samples <- 0
    repeat { # sampling in-bounds
      par <- current_par + rnorm(n = dim,mean = 0,sd = sigma)
      nb_samples <- nb_samples + 1
      if (all((par < UB) & (par > LB))) break()
      if (nb_samples > maxnbsamp) {
        # too many attempts at satisfying bounds, issue a warning and project on bounds
        warning("too many attempts at satisfying bounds, project on bounds", immediate. = TRUE, call. = TRUE)
        iviol <- par > UB
        par[iviol] <- UB[iviol]
        iviol <- par < LB
        par[iviol] <- LB[iviol]
        break()
      }
    }  
    value <- test_fun(par)
    x_hist[i, ] <- par
    y_hist[i] <- value
    # current point move (change here to implement simulated annealing)
    if (value > best_value) {
      current_par <- par
      current_value <- value
    }   
    # for archiving results
    if (value > best_value) {
      best_par <- par
      best_value <- value
    }   
  }
  res <- list(xhist=x_hist, fhist=y_hist, x_best=best_par, f_best=best_value)
  return(res)
}
