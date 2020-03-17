source("./cmaes.R")
source("./test_functions.R")
source("./RSalgorithm.R")
source("./NSalgorithm.R")
source("./lbfgs.R")
source("./of_wrapper.R")


antennes6d<-load("antennes_6d_train.Rdata")
C.train<-C
S.train<-S
GPmodel <- km(design=C.train, response=S.train ,covtype="matern3_2")

library(DiceKriging)
library("rgl") # library for plots

zdim <- 6
budget <- 3000
nb_antennas<- 3
LB = -6.5
UB = 8


#Creating the grid of 4^6 elements
x1=seq(-6,5,length=4)
x2=seq(-6.5,8.7,length=4)
grid <- expand.grid(x1=x1, x2=x2,x3=x1,x4=x2,x5=x1,x6=x2)
grid<-cbind(grid[,1],grid[,2],grid[,3],grid[,4],grid[,5],grid[,6])

#Applying the function meanofGP6D
y <- apply(grid, 1, meanofGP6D)

fun <- meanofGP6D
# optimize with normal search
paramNS <- list(LB=LB,UB = UB,budget = budget,dim=6, xinit=rep(-3.2,6),sigma=0.3) # param for normal_search
optresNS <- normal_search(fun, paramNS)
print(optresNS$x_best)
print(optresNS$f_best)


# controls for noisy functions and other dirty global variables
glob_noisy <- FALSE # is the function noisy
glob_tau <- 1 # noise std deviation
# glob_estim_noise <- FALSE # this should go in KNF parameters
glob_xstar <- rep(2.5,zdim)
store_hist <<- FALSE # TRUE only inside lbfgs.R, see file.

paramCMA <- list(LB=LB,UB = UB,budget = budget, dim=6, xinit=rep(-3.2,6),sigma=1.) # param for cmaes  
optresCMA <- cmaes(fun, paramCMA)
print(optresCMA$f_best)

par(mfrow=c(1,1))
# print out results
cat("xbest=")
optresNS$x_best
cat("fbest=")
optresNS$f_best
plot(optresNS$fhist,type="l",
     xlab="no. calls to f",ylab="f")
title("ES-(1+1)")

# print out results
cat("xbest=")
optresCMA$x_best
cat("fbest=")
optresCMA$f_best
plot(optresCMA$ymeanhist,type="l",
     xlab="no. of iterations",ylab="f of xmean")
title("CMA-ES")


# save the results, put your names in the output file
fname <- "Adref_Benechehab_Gueddari_cma_6d.Rdata"
x_solution_cma_6d <- matrix(optresNS$x_best,ncol=zdim)
y_solution_cma_6d <- predict(object = GPmodel,newdata=data.frame(x_solution_cma_6d),type="UK")
# y_solution_cma_2d$mean for GPmodel prediction 
# and y_solution_cma_2d$sd for GPmodel standard deviation at x_solution
# just saving the 2d solution, you'll have to add the 6d solution as well
save(x_solution_cma_6d,y_solution_cma_6d,file=fname)
# at the end for example, load("temporary_files.Rdata") to get back 
# other solutions (6d, with EGO) and save again in a final file 
# save(... all objects ..., file="Name1_Name2_Name3.Rdata")