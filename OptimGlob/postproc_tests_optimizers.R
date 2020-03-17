# Loads, plots and compares the results of several tests of optimizers
# Author : Rodolphe Le Riche
# Input : the names of the *.RData files containing each test results (matrix "fbhmat") 
# Output : comparison plots
# Assumption : the data files contain the list hist, which has the element hist[[i]]$fhist where i ranges from 1 to the no of experiments performed
# 

### user data
# filenames <- c("./cmaes_quadratic_5.RData","./lbfgs_quadratic_5.RData")
filenames <- c("./lbfgs_ackley_5.RData")
# filenames <- c("./random_search_sphere_5.RData","./normal_search_sphere_5_sigma_0.1.RData","./normal_search_sphere_5_sigma_1.RData","./normal_search_sphere_5_sigma_10.RData")
# filenames <- c("./normal_search_michalewicz_5_sigma_2.RData","./cmaes_michalewicz_5.RData","./lbfgs_michalewicz_5.RData")
# filenames <- c("./random_search_sphere_5.RData","./normal_search_sphere_5_sigma_1.RData","./cmaes_sphere_5.RData")
# filenames <- c("./random_search_sphere_5.RData","./normal_search_sphere_5_sigma_1.RData","./EGO_sphere_5.RData","./cmaes_sphere_5.RData")
# filenames <- c("./EGO_ackley_5.RData","./normal_search_ackley_5_sigma_5.RData")
# filenames <- c("./normal_search_ackley_5_sigma_5.RData","./lbfgs_restart_ackley_5.RData")
# filenames <- c("./cmaes_sphere_5.RData","./lbfgs_restart_sphere_5.RData")
thelegends <- filenames # could be changed for more explicit legends
### end user data

source(file="./utility_optim_functions.R")
no_files <- length(filenames)
fbhmat <- list()
no_test <- c()
for (i in 1:no_files) {
  load(filenames[i])
  fbhmat[[i]] <- t(apply(hist2hmat(hist),1,cummin))
  no_test <- c(no_test,dim(fbhmat[[i]])[1])
}

### plot on log scale of median and trajectories
theylim <- range(sapply(fbhmat,range))+c(-0.5,0.5)
if (theylim[1] <= 1.e-4) {theylim[1]<- 1.e-4}
theylim <- log(theylim)
thexlim <- c(0,max(sapply(fbhmat,dim)[2,]))

for (i in 1:no_files) {
  if (i==1){
    plot(log(apply(fbhmat[[i]], 2, median)),type="l", lwd=2.5, xlim=thexlim, ylim=theylim,col=i, ylab="log fbest",xlab="no. calls to f", cex.axis=0.9, cex.lab=0.9)    
  } else {
    lines(log(apply(fbhmat[[i]], 2, median)),type="l", lwd=2.5, ylim=theylim,col=i, ylab="log fbest",xlab="no. calls to f", cex.axis=0.9, cex.lab=0.9)
  }
  for (j in 1:no_test[i]) {
    lines(log(fbhmat[[i]][j, ]), type="l", lwd=0.5, col=i)
  }
}

txtwdth <- floor(dim(fbhmat[[1]])[2]/1.5)

legend(x="topright", legend=thelegends, col=seq(1,no_files),lty=rep(1, no_files) , cex=0.9, text.width=txtwdth)

### plot of mean +/- std deviation

for (i in 1:no_files) {
  if (i==1) {
    plotfmat(fbhmat[[i]],xlab="no. function evaluations",ylab=" fbest +/- std dev",xlim=thexlim, icol=i)
  } else {
    plotfmat(fbhmat[[i]],add=TRUE,icol=i)
  }
}
legend(x="topright", legend=thelegends, col=seq(1,no_files),lty=rep(1,no_files) , cex=1.2, text.width=txtwdth)

