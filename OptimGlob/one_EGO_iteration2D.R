# template for one EGO iteration
#
rm(list=ls())

library(DiceDesign)
library(DiceKriging)
library("rgl") # library for plots
source("./test_functions.R")
source("./RSalgorithm.R")
source("./NSalgorithm.R")
source("./lbfgs.R")
source("./cmaes.R")
source("./of_wrapper.R")

##### user data #####
## 2d problem 
load('antennes_2d_train.Rdata')  # or antennes_6d_train.Rdata
nb_antennas <- 1  # make it 3 for the 6d antennas pb
C_init <- C
S_init <- S
fmin <- min(S_init)
LB = -6
UB = 8
# calculate a GP model
GPmodel <- km(~1, design=C_init, response=S_init, covtype = "matern5_2")
zdim <- GPmodel@d
# x = runif(zdim, min=LB,max=UB)
# x <- matrix(x,ncol=GPmodel@d)
x1=seq(-6,5,length=10)
x2=seq(-6,8.7,length=10)
grid <- expand.grid(x1=x1, x2=x2)
grid<-cbind(grid[,1],grid[,2])
y <- apply(grid, 1, meanofGP)
p <- predict(object = GPmodel,newdata=data.frame(grid),type="UK")

##### define -EI function
#   GPmodel and fmin are passed as global variables... a bit ugly
mEI <- function(x) {
  x <- matrix(x,ncol=GPmodel@d)
  p <- predict(object = GPmodel,newdata=data.frame(x),type="UK")
  m <- p$mean
  s <- p$sd
  # !!!  here do the calculations for EI
  imp = -fmin+m
  w = imp / s
  if (s==0 ){
    EI=0
  }
  EI <- imp*pnorm(w)+s*dnorm(w)
  
  # !!!  end here do the calculations for EI
  return(EI)
  # return(as.numeric(-1*EI))
}
budget=3000

# controls for noisy functions and other dirty global variables
glob_noisy <- FALSE # is the function noisy
glob_tau <- 1 # noise std deviation
# glob_estim_noise <- FALSE # this should go in KNF parameters
glob_xstar <- rep(2.5,zdim)
store_hist <<- FALSE # TRUE only inside lbfgs.R, see file.

# optimize with CMA
paramCMA <- list(LB=LB,UB = UB,budget = budget, dim=zdim, xinit=rep(-4, zdim),sigma=2.) # param for cmaes  
optresCMA <- normal_search(mEI, paramCMA)
fhist<- apply(optresCMA$xhist, 1, meanofGP)
# print out results
cat("xbest=")
optresCMA$x_best
cat("fbest=")
optresCMA$f_best
meanofGP(optresCMA$x_best)
open3d()
surface3d(x1, x2, y, col= "lightblue")
points3d(optresCMA$xhist[,1], optresCMA$xhist[,2],fhist , pch=19, col="red", size=10)
title3d("mean of kriging using ES-(1+1)", col="blue", font=4)
decorate3d()
aspect3d(1, 1, 1)

# plot(optresCMA$ymeanhist,type="l",
#      xlab="no. of iterations",ylab="f of xmean")
# title("CMA-ES")

# save the results, put your names in the output file
fname <- "Adref_Benechehab_Gueddari_ego_2d.Rdata"
x_solution_ego_2d <- matrix(optresCMA$x_best,ncol=zdim)
y_solution_ego_2d <- predict(object = GPmodel,newdata=data.frame(x_solution_ego_2d),type="UK")
# y_solution_cma_2d$mean for GPmodel prediction 
# and y_solution_cma_2d$sd for GPmodel standard deviation at x_solution
# just saving the 2d solution, you'll have to add the 6d solution as well
save(x_solution_ego_2d,y_solution_ego_2d,file=fname)
# at the end for example, load("temporary_files.Rdata") to get back 
# other solutions (6d, with EGO) and save again in a final file 
# save(... all objects ..., file="Name1_Name2_Name3.Rdata")
