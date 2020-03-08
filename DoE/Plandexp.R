rm(list=ls())
load("antennes_2d_test.Rdata")
C_test <- C
load ("antennes_2d_train.Rdata")

plot(C[,1],C[,2], col="blue")
grid(nx = 20, ny = 20, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
plot(C_test[,1],C_test[,2], col="red")
grid(nx = 100, ny = 100, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

regx = seq(from=-6, to=5, length=5)
regy = seq(from=-6.5, to=8.7, length=5)
reg = matrix(nrow = 25,ncol = 2)
k=0
for (j in 1:5) {
  for (i in 1:5){
    k=k+1
    reg[k,1] = regx[i]
    reg[k,2] = regy[j]
  }
}

min1 <- -6
max1 <- 5
min2 <- -6.5
max2 <- 8.7

n <- 20
unif1 <- (max1-min1)*runif(n)+min1
unif2 <- (max2-min2)*runif(n)+min2

plot(reg[,1],reg[,2])
grid(nx = 20, ny = 20, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)

plot(unif1,unif2)
grid(nx = 20, ny = 20, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)

maximin <- function(X) {
  return(min(dist(X)))
}
print(c("maximin pour grille reguliere", maximin(reg)))
print(c("maximin pour grille donnée", maximin(C)))
print(c("maximin pour grille données test", maximin(C_test)))
print(c("maximin pour grille uniforme", maximin(c(unif1,unif2))))

install.packages("DiceDesign")
library("DiceDesign")
print(c("discrepancy pour grille reguliere", discrepancyCriteria(reg,type='L2')))
print(c("discrepancy pour grille donnée", discrepancyCriteria(C,type='L2')))
print(c("discrepancy pour grille uniforme", discrepancyCriteria(c(unif1,unif2),type='L2')))

object <- rss2d(C,c(-6,-6.5),c(5,8.7))

#Maximin
NewCall1 <- maximinSA_LHS(C, T0=1, c=0.95, it=1000, p=50, profile="GEOM", Imax=100)
NewC1 <- NewCall1$design
plot(NewC1[,1],NewC1[,2], col="blue",legend= "New Design")
points(C[,1],C[,2], col="red",legend= "First Design")
grid(nx = 20, ny = 20, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
legend(-6, 4, legend=c("New Design","First Design"),
       col=c("blue", "red"), lty=1:2, cex=0.8)
print(c("maximin pour grille améliorée avec maximin critère", maximin(NewC1)))

#Discrepancy
Discrepcrit <- function(nb_points,nb_candidates){
  X_final = matrix(nrow = nb_points,ncol = 2)
  unif <- cbind((max1-min1)*runif(2)+min1,(max2-min2)*runif(2)+min2)
  d <- discrepancyCriteria(unif,type='L2')
  d <- d$DisL2
  X_final[1:2,] <- unif
  for(i in 3:nb_points){
    set.seed(i)
    unif <- cbind((max1-min1)*runif(nb_candidates)+min1,(max2-min2)*runif(nb_candidates)+min2)
    k=1
    X_final[i,] <- unif[1,]
    d <- discrepancyCriteria(X_final[1:i,],type='L2')
    d <- d$DisL2
    mincrit <- d
    for (j in 2:nb_candidates) {
      X_final[i,] <- unif[j,]
      d <- discrepancyCriteria(X_final[1:i,],type='L2')
      d <- d$DisL2
      if (d < mincrit) {
        mincrit <- d
        k=j
      }
    }
    X_final[i,] <- unif[k,]
  }
  return(X_final)
}
discrepdesign <- Discrepcrit(20,5000)
points(discrepdesign[,1],discrepdesign[,2], col="blue")
plot(C[,1],C[,2], col="red")
grid(nx = 20, ny = 20, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
legend(-6, 5.5, legend=c("New Design (discrepancy)","First Design"),
       col=c("blue", "red"), lty=1:2, cex=0.8)

print(c("maximin pour grille améliorée avec discrepance critère", maximin(discrepdesign)))
print(c("discrepancy pour grille améliorée avec discrepance critère", discrepancyCriteria(discrepdesign,type='L2')))

#IMSE
library(DiceKriging)
#install.packages("KrigInv")
library("KrigInv")
Model <- km(design = C, response = S)
integcontrol <- list(distrib="imse",n.points=20,n.candidates=5000,init.distrib="MC")
integ.param <- integration_design(integcontrol=integcontrol,lower=c(-6,-6.5),upper=c(5,8.7), model=Model)
imsedesign <- integ.param$integration.points
plot(imsedesign[,1],imsedesign[,2], col="blue")
grid(nx = 20, ny = 20, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
print(c("discrepancy pour grille génerée avec le critère IMSE", discrepancyCriteria(imsedesign,type='L2')))
print(c("maximin pour grille génerée avec le critère IMSE", maximin(imsedesign)))

#MaxVar
MaxVar <- function(nb_points,nb_candidates,object_km){
  X_final = matrix(nrow = nb_points,ncol = 2)
  unif <- cbind((max1-min1)*runif(1)+min1,(max2-min2)*runif(1)+min2)
  Result <- predict.km(object = object_km, newdata = unif, type = "SK")
  sd <- Result$sd
  mincrit <- sd[1]
  X_final[1,] <- unif
  for(i in 2:nb_points){
    set.seed(i)
    unif <- cbind((max1-min1)*runif(nb_candidates)+min1,(max2-min2)*runif(nb_candidates)+min2)
    Result <- predict.km(object = object_km, newdata = unif, type = "SK")
    sd <- Result$sd
    k=0
    for (j in 1:nb_candidates) {
      if (sd[j] <= mincrit) {
        mincrit <- sd[j]
        k=1
        break
      }
    }
    if (k==1) {
      X_final[i,] <- unif[j,]
    } else {
      X_final[i,] <- unif[1,]
    }
  }
  return(X_final)
}
maxvardesign <- MaxVar(20,500,Model)
plot(maxvardesign[,1],maxvardesign[,2], col="blue")
grid(nx = 20, ny = 20, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
print(c("discrepancy pour grille génerée avec le critère MaxVar", discrepancyCriteria(maxvardesign,type='L2')))
print(c("maximin pour grille génerée avec le critère MaxVar", maximin(maxvardesign)))


#######################################
#6D

rm(list=ls())
load("antennes_6d_test.Rdata")
C_test <- C
load ("antennes_6d_train.Rdata")

# plot(C[,1],C[,2], col="blue")
# grid(nx = 300, ny = 300, col = "lightgray", lty = "dotted",
#      lwd = par("lwd"), equilogs = TRUE)
# plot(C_test[,1],C_test[,2], col="red")
# grid(nx = 300, ny = 300, col = "lightgray", lty = "dotted",
#      lwd = par("lwd"), equilogs = TRUE)
# plot(C[,3],C[,4], col="blue")
# grid(nx = 300, ny = 300, col = "lightgray", lty = "dotted",
#      lwd = par("lwd"), equilogs = TRUE)
# plot(C_test[,3],C_test[,4], col="red")
# grid(nx = 300, ny = 300, col = "lightgray", lty = "dotted",
#      lwd = par("lwd"), equilogs = TRUE)
# plot(C[,5],C[,6], col="blue")
# grid(nx = 300, ny = 300, col = "lightgray", lty = "dotted",
#      lwd = par("lwd"), equilogs = TRUE)
# plot(C_test[,5],C_test[,6], col="red")
# grid(nx = 300, ny = 300, col = "lightgray", lty = "dotted",
#      lwd = par("lwd"), equilogs = TRUE)

regx1 = seq(from=-6, to=5, length=3)
regx2 = seq(from=-6, to=5, length=3)
regx3 = seq(from=-6, to=5, length=3)
regx4 = seq(from=-6.5, to=8.7, length=3)
regx5 = seq(from=-6.5, to=8.7, length=3)
regx6 = seq(from=-6.5, to=8.7, length=3)
reg6d = expand.grid(regx1,regx2,regx3,regx4,regx5,regx6)

min1 <- -6
max1 <- 5
min2 <- -6.5
max2 <- 8.7

n <- 300
set.seed(55)
unif1 <- (max1-min1)*runif(n)+min1
unif2 <- (max1-min1)*runif(n)+min1
unif3 <- (max1-min1)*runif(n)+min1
unif4 <- (max2-min2)*runif(n)+min2
unif5 <- (max2-min2)*runif(n)+min2
unif6 <- (max2-min2)*runif(n)+min2

unif6d <- cbind(unif1,unif2,unif3,unif4,unif5,unif6)

print(c("maximin pour grille reguliere", maximin(reg6d)))
print(c("maximin pour grille donnée", maximin(C)))
print(c("maximin pour grille uniforme", maximin(unif6d)))

print(c("discrepancy pour grille reguliere", discrepancyCriteria(reg6d,type='L2')))
print(c("discrepancy pour grille donnée", discrepancyCriteria(C,type='L2')))
print(c("discrepancy pour grille uniforme", discrepancyCriteria(unif6d,type='L2')))










