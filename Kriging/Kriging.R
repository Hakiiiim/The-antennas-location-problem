#application
rm(list=ls())
load("antennes_6d_test.Rdata")
C_6d_test <- C
load ("antennes_6d_train.Rdata")
C_6d <- C
S_6d <- S
load("antennes_2d_test.Rdata")
C_test <- C
load ("antennes_2d_train.Rdata")


library(ModelMetrics)

######################################################################"
#Noyaux

#Noyau gaussian
kGauss <- function(x,y,Theta,Sigma,param=NULL){
  if(is.null(param)) param <-1
  Gaussian <- function(a,b){
    (Sigma^2)*exp(-(norm(as.matrix(a-b))^2)/(2*Theta^2))}
  
  n = dim(x)[1]
  m = dim(y)[1]
  A = matrix(nrow = n,ncol=m)
  for (i in 1:n){
    for (j in 1:m) {
      A[i,j] = Gaussian(x[i,],y[j,])
    }
  }
  return(A)
}

# Noyau exponentiel
kExp <- function(x,y,theta=NULL,sigma=NULL){
  if(is.null(theta)) {theta <-1}
  if(is.null(sigma)) {sigma <-1}
  exponentiel <- function(x,y){
    return((sigma^2)*exp(-norm(as.matrix(x-y))/theta))
  }
  n = dim(x)[1]
  m = dim(y)[1]
  A = matrix(nrow = n,ncol=m)
  for (i in 1:n){
    for (j in 1:m) {
      A[i,j] = exponentiel(x[i,],y[j,])
    }
  }
  return(A)
}

# #Noyau sinus (fonctionne pas)
# kSin <- function(x,y,Theta,Sigma,param=NULL){
#   if(is.null(param)) param <-1
#   Sin <- function(x,y){
#     (Sigma^2)*(Theta/norm(as.matrix(x-y)))*sin(norm(as.matrix(x-y))/(Theta))}
#   
#   n = dim(x)[1]
#   m = dim(y)[1]
#   A = matrix(nrow = n,ncol=m)
#   for (i in 1:n){
#     for (j in 1:m) {
#       A[i,j] = Sin(x[i,],y[j,])
#     }
#   }
#   return(A)
# }

# Noyau matern 5/2
kMatern5 <- function(x,y,theta=NULL,sigma=NULL){
  if(is.null(theta)) {theta <-1}
  if(is.null(sigma)) {sigma <-1}
  matern5 <- function(x,y){
    return((sigma^2)*(1+(sqrt(5)/theta)*norm(as.matrix(x-y)) + (5/3*theta^2)*(norm(as.matrix(x-y)))^2)*exp(-sqrt(5)*norm(as.matrix(x-y))/theta))
  }
  n = dim(x)[1]
  m = dim(y)[1]
  A = matrix(nrow = n,ncol=m)
  for (i in 1:n){
    for (j in 1:m) {
      A[i,j] = matern5(x[i,],y[j,])
    }
  }
  return(A)
}


# Noyau matern 3/2
kMatern3 <- function(x,y,theta=NULL,sigma=NULL){
  if(is.null(theta)) {theta <-1}
  if(is.null(sigma)) {sigma <-1}
  matern3 <- function(x,y){
    return((sigma^2)*(1+(sqrt(3)/theta)*norm(as.matrix(x-y)))*exp(-sqrt(3)*norm(as.matrix(x-y))/theta))
  }
  n = dim(x)[1]
  m = dim(y)[1]
  A = matrix(nrow = n,ncol=m)
  for (i in 1:n){
    for (j in 1:m) {
      A[i,j] = matern3(x[i,],y[j,])
    }
  }
  return(A)
}

###################################################################
# leave-one-out quality measure

#fonction qui rÃ©alise une validation croisÃ©e en renvoyant la somme des erreurs en enlevant 20 un seul point
test <- function(kernel,theta,sigma,C,S) {
  ngrid <- n.expl
  xgrid <- C
  
  pred = rep(0,ngrid)
  
  #On enleve à chaque fois un element
  for (i in 1:ngrid) {
    #donnees d'apprentissage
    X <- xgrid[-i,]
    y <- S[-i]
    if (kernel == "gauss"){
      kXX <- kGauss(X,X,theta,sigma)
      invkXX <- solve(kXX)
      m <- function(x){
        kxx <- kGauss(x,x,theta,sigma)
        kxX <- kGauss(x,X,theta,sigma)
        return(kxX%*%invkXX%*%y)
      }
      mu <- m(xgrid)
    } else if (kernel == "exponentiel") {
      kXX <- kExp(X,X,theta,sigma)
      invkXX <- solve(kXX)
      
      m <- function(x){
        kxx <- kExp(x,x,theta,sigma)
        kxX <- kExp(x,X,theta,sigma)
        
        return(kxX%*%invkXX%*%y)
      }
      
      mu <- m(xgrid)
    } else if (kernel == "matern5") {
      kXX <- kMatern5(X,X,theta,sigma)
      invkXX <- solve(kXX)
      
      m <- function(x){
        kxx <- kMatern5(x,x,theta,sigma)
        kxX <- kMatern5(x,X,theta,sigma)
        
        return(kxX%*%invkXX%*%y)
      }
      
      mu <- m(xgrid)
    } else if (kernel == "matern3") {
      kXX <- kMatern3(X,X,theta,sigma)
      invkXX <- solve(kXX)
      
      m <- function(x){
        kxx <- kMatern3(x,x,theta,sigma)
        kxX <- kMatern3(x,X,theta,sigma)
        
        return(kxX%*%invkXX%*%y)
      }
      
      mu <- m(xgrid)
    } 
    pred[i] <- mu[i]
  }
  return(rmse(S,pred))
} 


theta = matrix(seq(from=0.01, to=10, length=100),ncol=1)
sigma = matrix(seq(from=0.01, to=10, length=10),ncol=1)
kernel = c("gauss","exponentiel","matern3","matern5")

Rmse=100000

#reste à trouver l'Ã©quivalent des autres noyaux dans le cas vectoriel 
#erreur pour le noyau sinus cardinal : renvoi de valeurs infinies


#le temps de calcul est important
for (k in 1:length(kernel)) { #juste gauss et exp parce que les autres ne sont pas encore codÃ©s
  for (j in 1:length(sigma)) {
    for (i in 1:length(theta)) {
      result = test(kernel = kernel[k], theta[i], sigma[j],C,S)
      if (result <= Rmse){
        Rmse = result
        thetabest=theta[i]
        sigmabest=sigma[j]
        kernelbest=kernel[k]
      }
    }
  }
}


print(c("meilleure valeure de theta :",thetabest))
print(c("meilleure valeure de sigma :",sigmabest))
print(c("meilleur noyau :",kernelbest))
print(c("l'erreur quadratique moyenne :",Rmse))

theta <- 4.5509
sigma <- 3.34
#kernel : matern3/2
#rmse : 1.8459


################################################
# Pédiction manuelle test

kXX <- kMatern3(C,C,theta,sigma)
invkXX <- solve(kXX)

m <- function(x){
  kxx <- kMatern3(x,x,theta,sigma)
  kxX <- kMatern3(x,C,theta,sigma)
  
  return(kxX%*%invkXX%*%S)
}

covcond <- function(x1,x2){
  kxxp <- kMatern3(x1,x2,theta,sigma)
  kxX1 <- kMatern3(x1,C,theta,sigma)
  kxX2 <- kMatern3(x2,C,theta,sigma)
  
  return(kxxp - kxX1%*%invkXX%*%t(kxX2))
}

m_test <- m(C_test)

fileConn<-file("output.txt")
writeLines(toString(t(m_test)), fileConn)
close(fileConn)

text <- read.table("output.txt")
print(text)

###############################################
# Projet

#1
#load data

#2
install.packages("DiceKriging")
library(DiceKriging)

Model_2d <- km(formula=~1,design = C, response = S, covtype="matern3_2")

Result <- predict.km(object = Model_2d, newdata = C_test, type = "SK")

Couv_2d_mean_krg <- Result$mean
Couv_2d_sdv_krg <- Result$sd

Model_6d <- km(formula=~1,design = C_6d, response = S_6d, covtype="matern3_2")

Result <- predict.km(object = Model_6d, newdata = C_6d_test, type = "SK")

Couv_6d_mean_krg <- Result$mean
Couv_6d_sdv_krg <- Result$sd

fname <- "Benechehab_Gueddari_Adref_cma_2d.Rdata"
save(Couv_2d_mean_krg,Couv_2d_sdv_krg,Couv_6d_mean_krg,Couv_6d_sdv_krg,file=fname)
#3
show(Model_2d)

plot(Model_2d)

#4

#########################"
#plot surface

#regular

#set.seed(55)
#regx = (5-(-6))*runif(25)-6
#regy = (8-(-6))*runif(25)-6
regx = seq(from=-6, to=5, length=25)
regy = seq(from=-6, to=8, length=25)
reg = matrix(nrow = 625,ncol = 2)
k=0
for (j in 1:25) {
  for (i in 1:25){
    k=k+1
    reg[k,1] = regx[i]
    reg[k,2] = regy[j]
  }
}

#uniform
min1 <- -6
max1 <- 5
min2 <- -6.5
max2 <- 8.7
n <- 100
X1 <- sort((max1-min1)*runif(n)+min1)
X2 <- sort((max2-min2)*runif(n)+min2)
X = matrix(nrow = n^2,ncol = 2)
k=0
for (j in 1:n) {
  for (i in 1:n){
    k=k+1
    X[k,1] = X1[i]
    X[k,2] = X2[j]
  }
}
Result <- predict.km(object = Model_2d, newdata = X, type = "SK")
#Prediction and confidence interval
mu <- Result$mean
up <- Result$upper95
low <- Result$lower95
library(rgl)
open3d()
surface3d(X1,X2,mu,color = "red")
surface3d(X1,X2,up,color = "black", alpha = 0.5, lit = FALSE, front = "lines", back = "lines")
surface3d(X1,X2,low,color = "black", alpha = 0.5, lit = FALSE, front = "lines", back = "lines")
axes3d()







