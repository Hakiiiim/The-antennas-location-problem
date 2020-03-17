## set of test functions
antennes2d<-load("antennes_2d_train.Rdata")
X.train<-C
(X.train)

# Ackley function. Global minimum at glob_xstar
ackley <- function(xx, a=20, b=0.2, c=2*pi){
  if (glob_noisy) {
    noise <-  rnorm(1, 0, glob_tau) 
  } else noise <- 0
  xx <- xx - glob_xstar
  aa <- 6.4*xx
  d <- length(aa)
  sum <- -a*exp(-b*sqrt((1/d)*sum(aa*aa))) - exp((1/d)*sum(cos(c*aa)))+ a + exp(1)
  y <- sum + noise
  return(y)
}

## Rastrigin function. Global minimum at glob_xstar
rastrigin <- function(xx){
  if (glob_noisy) {
    noise <-  rnorm(1, 0, glob_tau) 
  } else noise <- 0
  scaling_fact <- 1.024
  xx <- xx - glob_xstar
  aa <- scaling_fact*xx
  d <- length(aa)
  sum <- sum(aa*aa - 10*cos(2*pi*aa))
  y <- sum + 10*d
  return(y + noise)
}

## Schwefel function.

schwefel <- function(xx){
  if (glob_noisy) {
    noise <-  rnorm(1, 0, glob_tau) 
  } else noise <- 0
  xx <- xx - glob_xstar + 1
  aa <- 100*xx
  d <- length(aa)
  sum <- 418.9829*d - sum(aa*sin(sqrt(abs(aa))))
  y <- sum
  return(y + noise)
}

## Sphere function. Global Minimum at glob_xstar

sphere <- function(xx){
  if (glob_noisy) {
    noise <-  rnorm(1, 0, glob_tau) 
  } else noise <- 0
  xx <- xx - glob_xstar
  #aa <- 1.024*xx
  aa <- xx
  d <- length(aa)
  sum <- sum(aa*aa)
  y <- sum
  return(y + noise)
}

##### michalewicz function #########
#glob min in 2D: -1.83, 5D: -4.71 , 10D: -9.64
michalewicz <- function(xx, m=10){
  if (glob_noisy) {
    noise <-  rnorm(1, 0, glob_tau) 
  } else noise <- 0
  aa <- 0.1*pi*(xx + 5)
#   aa <- xx
  d <- length(aa)
  i <- 1:d
  sum <- sum(sin(aa)*sin((i*aa*aa)/pi)^(2*m))
  y <- -sum
  return(y + 2 + noise)
}

##### quadratic function #########"
quadratic <- function(xx){
  if (glob_noisy) {
    noise <-  rnorm(1, 0, glob_tau) 
  } else noise <- 0
  cond.no <- 1.e3
  aa <- 1.024*xx
  d <- length(aa)
  #  xstar <- 2.5 # default value, change it if needed
  xstar <- glob_xstar
  lambdas <- diag(seq(1, cond.no, ,d))
  # matrix with arbitrary orientation. The seed number decides the orientation.
  if (!exists("glob_umat")) {
    set.seed(1) # change this seed to change the orientation of the quadratic function
    glob_umat <<- qr.Q(qr(matrix(runif(d*d),nrow=d,ncol=d)))
    glob_umat <<- glob_umat[,sample(seq(1,d))]
    #  glob_umat <<- diag(1, d) # to generate a quadratic function whose principal axes are aligned with coordinates
  }
  H <- glob_umat%*%lambdas%*%t(glob_umat)  
  y <- 0.5*t(aa - xstar)%*%H%*%(aa - xstar)
  return(y + noise)
}

#### Tunnel function #######"
tunnel <- function(xx, b=0.5){
  if (glob_noisy) {
    noise <-  rnorm(1, 0, glob_tau) 
  } else noise <- 0
  d <- length(xx)
  xx <- xx - glob_xstar
  in_tunnel <- TRUE
  if (d > 1) {
    for(i in 2:d){
      if ((xx[i] < -b) | (xx[i] > b)){
        in_tunnel <- FALSE
      }
    }
    if (in_tunnel) {
      y <- -exp(-norm(as.matrix(xx), type="F")/10)
    } else{
      y <- t(xx)%*%xx
    }
  } else{
    y <- xx^2
  }
  return(y+noise)
}


### 
##### function with hierarchical sensitivities and local optima #########
quad_wave <- function(xx){
  print("grid")
  print(xx)
  print("train")
  print(dim(X.train))
  cond.no <- 1.e1
  d <- length(xx)
  lambdas <- diag(seq(1, cond.no, ,d))
  glob_umat <<- diag(1, d) # to generate a quadratic function whose 
      # principal axes are aligned with coordinates. Could be made into a rotating matrix.
  H <- glob_umat%*%lambdas%*%t(glob_umat)  
  # y <- 10+0.5*t(xx)%*%H%*%(xx)-10*cos(2*pi*sqrt(t(xx)%*%xx))
  # y <- 10+0.5*t(xx)%*%H%*%(xx)-10*cos(2*pi*norm(x = as.matrix(xx),type = "I"))
   y <- 1+0.5*t(xx)%*%H%*%(xx)-0.2*sum(cos(4*pi*xx[c(1,2)]))
  return(y)
}

krigemean <- function(x){
  x <- matrix(x,ncol=2)
  print("grid")
  print(x)
  print("train")
  print(dim(X.train))
  antennes2d<-load("antennes_2d_train.Rdata")
  X.train<-C
  S.train<-S
  m1 <- km(design=X.train, response=S.train ,covtype="matern3_2")
  p1<- predict(m1, newdata=x , type="SK")
  y <- p1[["mean"]]
  return(mean)
}


antennes2d<-load("antennes_2d_train.Rdata")
X.train<-C
S.train<-S
GPmodel <- km(design=X.train, response=S.train ,covtype="matern3_2")
meanofGP <- function(x){
  x <- matrix(x,ncol=2)
  z <- predict(object = GPmodel,newdata=data.frame(x),type="SK")
  return(z$mean)
}

antennes6d<-load("antennes_6d_train.Rdata")
X.train<-C
S.train<-S
GPmodel <- km(design=X.train, response=S.train ,covtype="matern3_2")
meanofGP6D <- function(x){
  x <- matrix(x,ncol=6)
  z <- predict(object = GPmodel,newdata=data.frame(x),type="SK")
  return(z$mean)
}
