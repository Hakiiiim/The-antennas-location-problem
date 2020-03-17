source("./cmaes.R")
source("./test_functions.R")
source("./RSalgorithm.R")
source("./NSalgorithm.R")
source("./lbfgs.R")
source("./of_wrapper.R")

antennes2d<-load("antennes_2d_train.Rdata")
X.train<-C
S.train<-S
summary(antennes2d)
library(DiceKriging)
library("rgl") # library for plots
nb_antennas=1
zdim <- 2
budget <- 3000
LB = -6
UB = 8
x1=seq(-6,5,length=10)
x2=seq(-6,8.7,length=10)
grid <- expand.grid(x1=x1, x2=x2)
grid<-cbind(grid[,1],grid[,2])
# grid<-as.matrix(grid)
# y<- meanofGP(grid)
# y <- apply(grid, 1, krigemean)
GPmodel <- km(design=X.train, response=S.train ,covtype="matern3_2")
y <- apply(grid, 1, meanofGP)

# p <- predict(m1, X.test , type="SK")
fun <- meanofGP
# optimize with normal search
paramNS <- list(LB=LB,UB = UB,budget = budget,dim=zdim, xinit=rep(-3.2,2),sigma=0.3) # param for normal_search
optresNS <- normal_search(fun, paramNS)

open3d()
surface3d(x1, x2, y, col= "lightblue")
points3d(optresNS$xhist[,1], optresNS$xhist[,2], optresNS$fhist, pch=19, col="red", size=10)
title3d("mean of kriging using ES-(1+1)", col="blue", font=4)
decorate3d()
aspect3d(1, 1, 1)

# controls for noisy functions and other dirty global variables
glob_noisy <- FALSE # is the function noisy
glob_tau <- 1 # noise std deviation
# glob_estim_noise <- FALSE # this should go in KNF parameters
glob_xstar <- rep(2.5,zdim)
store_hist <<- FALSE # TRUE only inside lbfgs.R, see file.



paramCMA <- list(LB=LB,UB = UB,budget = budget, dim=zdim, xinit=rep(-3.3,2),sigma=1.) # param for cmaes  
optresCMA <- cmaes(fun, paramCMA)

open3d()
surface3d(x1, x2, y, col= "lightblue")
points3d(optresCMA$xmeanhist[,1], optresCMA$xmeanhist[,2], optresCMA$ymeanhist, pch=19, col="red", size=10)
title3d("meanofGP using CMAES", col="blue", font=4)
decorate3d()
aspect3d(1, 1, 1)


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
fname <- "Adref_Benechehab_Gueddari_cma_2d.Rdata"
x_solution_cma_2d <- matrix(optresNS$x_best,ncol=zdim)
y_solution_cma_2d <- predict(object = GPmodel,newdata=data.frame(x_solution_cma_2d),type="UK")
# y_solution_cma_2d$mean for GPmodel prediction 
# and y_solution_cma_2d$sd for GPmodel standard deviation at x_solution
# just saving the 2d solution, you'll have to add the 6d solution as well
save(x_solution_cma_2d,y_solution_cma_2d,file=fname)
# at the end for example, load("temporary_files.Rdata") to get back 
# other solutions (6d, with EGO) and save again in a final file 
# save(... all objects ..., file="Name1_Name2_Name3.Rdata")


