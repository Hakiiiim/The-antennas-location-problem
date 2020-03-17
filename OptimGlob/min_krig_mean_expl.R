# Correction: maximization of a kriging mean
#
# For the antenna project: replace Xinit and Yinit 
# (which come here from one of the functions used in the TP) 
# by data coming from the antenna files. zdim, LB, UB will change.
# Also, pay attention to the fact that the surface is maximized 
# while optimizers by default minimize.
# Then this answers pretty much the question 7.1.
# <== there is an example to cut and paste in one_EGO_iteration.R
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
zdim <- 2
budget <- 3000
LB = -5
UB = 5

# controls for noisy functions and other dirty global variables
glob_noisy <- FALSE # is the function noisy
glob_tau <- 1 # noise std deviation
# glob_estim_noise <- FALSE # this should go in KNF parameters
glob_xstar <- rep(2.5,zdim)
store_hist <<- FALSE # TRUE only inside lbfgs.R, see file.


# create a LHS DoE
Xinit <- data.frame(LB + (UB-LB)*lhsDesign(n = 4*zdim,dimension = zdim)$design)
# calculate associated objective function
Yinit <- apply(X = Xinit,MARGIN = 1,FUN = rastrigin)
# make a kriging model. 
# I use a fairly high lower bound on the thetas here to prevent degenrated models
GPmodel <-km(design = Xinit,response = Yinit,covtype="matern3_2",lower = rep(0.8,zdim),multistart = 20)

meanofGP <- function(x){
  x <- matrix(x,ncol=zdim)
  z <- predict(object = GPmodel,newdata=data.frame(x),type="UK")
  return(z$mean)
}

# One evaluation of the function, for example the point (1,1)
meanofGP(matrix(c(1,1),nrow=1))

# One MINIMIZATION of the kriging mean
# set the global fun variable to the meanofGP
fun <- meanofGP

# optimize with normal search
paramNS <- list(LB=LB,UB = UB,budget = budget,dim=zdim, xinit=rep(1,zdim),sigma=0.3) # param for normal_search
optresNS <- normal_search(fun, paramNS)
# print out results
cat("xbest=")
optresNS$x_best
cat("fbest=")
optresNS$f_best
plot(optresNS$fhist,type="l",
     xlab="no. calls to f",ylab="f")
title("ES-(1+1)")

# optimize with CMA-ES
paramCMA <- list(LB=LB,UB = UB,budget = budget, dim=zdim, xinit=rep(-4, zdim),sigma=2.) # param for cmaes  
optresCMA <- cmaes(fun, paramCMA)
# print out results
cat("xbest=")
optresCMA$x_best
cat("fbest=")
optresCMA$f_best
plot(optresCMA$ymeanhist,type="l",
     xlab="no. of iterations",ylab="f of xmean")
title("CMA-ES")

# save the results, put your names in the output file
fname <- "Name1_Name2_Name3_cma_2d.Rdata"
x_solution_cma_2d <- matrix(optresCMA$x_best,ncol=zdim)
y_solution_cma_2d <- predict(object = GPmodel,newdata=data.frame(x_solution_cma_2d),type="UK")
# y_solution_cma_2d$mean for GPmodel prediction 
# and y_solution_cma_2d$sd for GPmodel standard deviation at x_solution
# just saving the 2d solution, you'll have to add the 6d solution as well
save(x_solution_cma_2d,y_solution_cma_2d,file=fname)
# at the end for example, load("temporary_files.Rdata") to get back 
# other solutions (6d, with EGO) and save again in a final file 
# save(... all objects ..., file="Name1_Name2_Name3.Rdata")

# interactions with matlab
# you can also save in a matlab format, 
# which can then be loaded in matlab with load("filename.mat") 
library('R.matlab')
writeMat("Name1_Name2_Name3_cma_2d.mat",x_solution_cma_2d=x_solution_cma_2d,y_solution_cma_2d=y_solution_cma_2d)
