# Example file to show 2D functions in 3D plots
# Assume here that dim = 2 (dim is the number of variables)
# R. Le Riche
rm(list=ls())

library("rgl") # library for plots
source("./test_functions.R")
source("./RSalgorithm.R")
source("./NSalgorithm.R")
source("./lbfgs.R")
source("./cmaes.R")

##### user data #####
nameoffun <- "meanofGP" # "quadratic", "ackley", "michalewicz", "sphere", "rastrigin", "schwefel", "tunnel" ,"meanofGP" ...
nameofoptim <- "cmaes" # "random_search" , "normal_search", "lbfgs", "cmaes"

###parameters of the optimizer as a list. 
zdim <- 2
budget <- 30000
# parameters of the optimizer as a list
if (nameofoptim == "random_search") {
  param <- list(LB=-5,UB = 5,budget = budget,dim=zdim) # param for random_search
} else if (nameofoptim == "normal_search") {
  param <- list(LB=-5,UB = 5,budget = budget,dim=zdim, xinit=c(-3.2,3.2),sigma=0.3) # param for normal_search
} else if (nameofoptim == "cmaes") {
  param <- list(LB=-6,UB = 8,budget = budget,dim=zdim, xinit=c(0,0),sigma=1.) # cmaes
} else if (nameofoptim == "lbfgs") {
  param <- list(LB=-5,UB = 5,maxit = 100,dim=zdim, xinit=c(-3.2,3.2), trace=5, lmm =1, factr=100,pgtol=1.e-2) # param for lbfgs
} else {
  stop("unknown nameofoptim: ",nameofoptim,sep="")
}

# controls for noisy functions
glob_noisy <- FALSE # is the function noisy
glob_tau <- 1 # noise std deviation
# glob_estim_noise <- FALSE # this should go in KNF parameters
glob_xstar <- rep(2.5,zdim)

##### end user data #####

store_hist <<- FALSE # TRUE only inside lbfgs.R, see file.

eval(parse(text=paste("fun <- ",nameoffun)))
eval(parse(text=paste("opt <- ",nameofoptim)))

source("./of_wrapper.R")
no.grid <- 10
x <- seq(param$LB, param$UB, ,no.grid)
x.grid <- expand.grid(x, x)
z <- apply(x.grid, 1, fun)
z.grid <- matrix(z, no.grid)

#
# start an optimization
cat("Running optimizer ", nameofoptim,"\n")
optres <- opt(ofwrapper, param)
cat("Done with ", nameofoptim,"\n")

#### 3D rgl plot
open3d()
surface3d(x, x, z.grid, col= "lightblue")
points3d(optres$xhist[,1], optres$xhist[,2], optres$fhist, pch=19, col="red", size=10)
title3d("mean of kriging using CMA-ES", col="blue", font=4)
decorate3d()
aspect3d(1, 1, 1)
# text3d(optres$xhist[,1],optres$xhist[,2],optres$fhist, texts=1:length(optres$xhist[,1]), adj=c(1,1), col="green")
rgl.snapshot("./fileofplot.png", fmt="png", top=T)
# 


par(mfrow=c(1,1))
# print out results
cat("xbest=")
optres$x_best
cat("fbest=")
optres$f_best
plot(optres$fhist,type="l",
     xlab="no. calls to f",ylab="f")
title("dakchi")

### 2D contour plot
# save the contour in the current directory
png(filename="./contour.png") 
 contour(x, x, z.grid, nlevels=20, xlab="x", ylab="y")
 points(optres$xhist[,1], optres$xhist[,2], pch=20, col="blue")
 points(param$xinit[1], param$xinit[2], pch=19, col="red")
dev.off()
# control the number of plots on the same window with "par"
par(mfrow=c(1,3))  
contour(x, x, z.grid, nlevels=20, xlab="x1", ylab="x2", col="green")
points(optres$xhist[,1], optres$xhist[,2], pch=20, col="blue")
text(optres$xhist, labels=1:length(optres$xhist[,1]), pos=3, cex=1.0)     # (un)comment for labeling (or not) nb of calls to f when points created
points(optres$xhist[1,1], optres$xhist[1,2], pch=19, col="red") # initial point drawn in red

# plot of fhist and fbest vs nbcalls
matplot(optres$fhist, type="o", pch=19, ylab = "f_hist",xlab="nbcalls")
matplot(apply(optres$fhist, 2, cummin), type="o", pch=19, ylab = "f_best",xlab="nbcalls")

