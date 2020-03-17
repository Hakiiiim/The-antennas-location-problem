cmaes <- function (test_fun, param){
  # User defined input parameters 
  xmean <- param$xinit     # initial point
  N <- length(xmean)       # problem dimension
  if (length(param$LB) == 1) 
    lower <- rep(param$LB, N)
  if (length(param$UB) == 1) 
    upper <- rep(param$UB, N)
  sigma <- param$sigma     # step size
  # Strategy parameter setting
  if (is.null(param$lambda)) {
    lambda <- 4 + floor(3 * log(N))        # default population size, offspring number 
  }
  if (is.null(param$mu)) {
    mu <- floor(lambda/2)                  # number of parents, points for recombination
  }
  maxiter <- ceiling( param$budget / lambda)
  weights <- log(mu + 1) - log(1:mu)     # recombination weights
  weights <- weights/sum(weights)        # normalize recombination weights
  mueff <- sum(weights)^2/sum(weights^2) # variance effective
  cc <- (4 + mueff/N) / (N + 4 + 2*mueff/N)  # cumulation for C control, approximation is 4/N   
  cs <- (mueff + 2)/(N + mueff + 5)      # cumulation for sigma control
  c1 <- 2 / ((N+1.3)^2+mueff)      # approximation is 2/N^2, learning rate for rank-one update of C
  cmu <- 2 * (mueff - 2 + 1/mueff) / ((N + 2)^2 + 2*mueff/2)  # approximation is 0.3*lambda/N^2    # learning rate for rank-mu update of C
  damps <-  1 + 2*max(0, sqrt((mueff - 1)/(N + 1)) - 1) + cs # damping for sigma 
  # Checking the algorithm initial parameters
  stopifnot(length(upper) == N)
  stopifnot(length(lower) == N)
  stopifnot(all(lower < upper))
  stopifnot(length(sigma) == 1)
  best.fit <- Inf
  best.par <- NULL
  sigma.log <- numeric(maxiter)
  eigen.log <- matrix(, maxiter, N)
  xmeanhist <- matrix(, maxiter,N)
  ymeanhist <- numeric(maxiter)
  
  if ( (is.null(param$autolog)) || (param$autolog==TRUE) ) {
    store_hist <<- TRUE # start history recording
    nbcalls <<- 0
    glob_xhist<<- matrix(, 1, ncol = N)
    glob_fhist<<- matrix(, 1, 1)
  }

  pc <- rep(0, N) ; ps <- rep(0, N)  # evolution paths for C and sigma
  B <- diag(N)       # matrix of eigenvectors of C
  D <- diag(N)       # diagonal matrix of eigenvalues of C
  BD <- B %*% D      
  C <- BD %*% t(BD)  # covariance matrix
  chiN <- sqrt(N) * (1 - 1/(4 * N) + 1/(21 * N^2))  # expectation of ||N(0, I)||
  iter <- 0L
  counteval <- 0L
  cviol <- 0L
  msg <- NULL
  arx <- matrix(0, nrow = N, ncol = lambda)
  arfitness <- numeric(lambda)
  while (iter < maxiter) {
    iter <- iter + 1L
    sigma.log[iter] <- sigma
    # Generate lambda offspring 
    arz <- matrix(rnorm(N * lambda), ncol = lambda) # matrix of size N*(lambda) contains standard normally distributed vectors
    arx <- xmean + sigma * (BD %*% arz)
    vx <- ifelse(arx > lower, ifelse(arx < upper, arx, upper), 
                 lower)               # enforce the sample points to be inside the bounds 
    pen <- 1 + colSums((arx - vx)^2)  # penalty of the offspring inside the bounds is 1, otherwise >1  
    pen[!is.finite(pen)] <- .Machine$double.xmax/2
    y <- apply(vx, 2, test_fun)
    y<- -y
   
    counteval <- counteval + lambda
    arfitness <- y * pen

    arindex <- order(arfitness)
    arfitness <- arfitness[arindex]
    aripop <- arindex[1:mu]
    selx <- arx[, aripop]
    xmean <- drop(selx %*% weights)
    selz <- arz[, aripop]
    zmean <- drop(selz %*% weights)
    ymean <- drop(y[aripop] %*% weights)

    xmeanhist[iter,] <- xmean
    ymeanhist[iter] <- ymean
    # save the bests, 2 cases, noisy or not
    if (glob_noisy) {
      if (ymean < best.fit) {
        best.fit <- ymean
        best.par <- xmean
      }
    }
    else {
      valid <- pen <= 1
      if (any(valid)) {
        wb <- which.min(y[valid])
        if (y[valid][wb] < best.fit) {
          best.fit <- y[valid][wb]
          best.par <- arx[, valid, drop = FALSE][, wb]
        }
      }
    }


    norm <- function(x) drop(sqrt(crossprod(x)))
    ps <- (1 - cs) * ps + sqrt(cs * (2 - cs) * mueff) * (B %*% zmean)
    hsig <- drop((norm(ps)/sqrt(1 - (1 - cs)^(2 * counteval/lambda))/chiN) < (1.4 + 2/(N + 1)))
    pc <- (1 - cc) * pc + hsig * sqrt(cc * (2 - cc) * mueff) * drop(BD %*% zmean)
    BDz <- BD %*% selz
    # update of C
    C <- (1 - c1 -cmu) * C + c1 * (pc %o% pc + (1 - hsig) * cc * (2 - cc) * C) + cmu * BDz %*% diag(weights) %*% t(BDz)
    # Adapt step size sigma
    sigma <- sigma * exp((norm(ps)/chiN - 1) * cs/damps)
    e <- eigen(C, symmetric = TRUE)
    eigen.log[iter, ] <- rev(sort(e$values))
    if (!all(e$values >= sqrt(.Machine$double.eps) * abs(e$values[1]))) {
      msg <- "Covariance matrix 'C' is numerically not positive definite."
      break
    }
    B <- e$vectors
    D <- diag(sqrt(e$values), length(e$values))
    BD <- B %*% D
   }
  # add record and plot of sigma and eigenvalues
  res <- list(xhist=glob_xhist, fhist=glob_fhist, x_best=best.par, f_best= -best.fit, sigmahist=sigma.log, xmeanhist=xmeanhist, ymeanhist=-ymeanhist, eigenhist=eigen.log)

  if ( (is.null(param$autolog)) || (param$autolog==TRUE) ) {
    store_hist <<- FALSE  # end history recording
  }

  return(res)
}
