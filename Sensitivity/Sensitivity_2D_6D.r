#---------------------------------#
# Projet fil rouge : SA_antennes. #
# ADREF - BENCHHAB - GUEDARI      #
# Sensitivity part 2D / 6D        #
#---------------------------------#



#------------------------------------------------------------------------------------#
#                          Sensitivity for 2D                                        #
#------------------------------------------------------------------------------------#

# Loading data 
load("antennes_2d_train.Rdata")

library(sensitivity)
library(DiceKriging)

# Create the kriging model using the C, response S and the kernel = matern3_2 that we got at metamodeling part
modelk = km(S~., design= data.frame(C) ,response = S,covtype="matern3_2")

# Define the function that will predict for us the couvred surface Sn using our modelK, and the Xn (the data that we will genrate by later)
S_tild = function(xn , model = modelk)
{
  t= predict(modelk , newdata = data.frame(xn),checkNames = FALSE , type = "UK")
  t$mean
  return (t$mean)
}

# Generate our data that follow a Uniform distribution on [LB, UB]
X1 <- data.frame(matrix(runif(n=400,min =  -6, max = 5), nrow = 400))
X2 <- data.frame(matrix(runif(n=400,min =  -6.5, max = 8.7), nrow = 400))
X <- cbind(X1,X2)

# Geting Sn surface couvred using S_tild() that we created 
Sn = S_tild(X, modelk)

# S0 = E(S(X)) ~ mean(Sn)
mu0 <- mean(Sn)

# ------------------------
# Decomposition de Sobol
# ------------------------

# Pour X1
cat("Projection sur la variable pour X1")
plot(Sn - mu0 ~ X[,1], xlab=paste("X", 1, sep=""), ylab=expression(f(X)-mu[0]))

cat("Esperance conditionnelle simulee") ; 
m.ss1 <- smooth.spline(X[,1], Sn-mu0)
lines(m.ss1, lwd=10, col="blue")


# Pour X2
cat("Projection sur la variable pour X2")
plot(Sn - mu0 ~ X[,2], xlab=paste("X", 2, sep=""), ylab=expression(f(X)-mu[0]))

cat("Esperance conditionnelle simulee") ; 
m.ss2 <- smooth.spline(X[,2], Sn-mu0)
lines(m.ss2, lwd=10, col="red")

cat("variance de S1 est : ",var(m.ss1$y))
cat("variance de S2 est : ",var(m.ss2$y) )

cat("Indices de Sobol") ;

SA.Antenne <- fast99(model = S_tild, 
                     factors = 2, 
                     n = 1000,
                     q = "qunif", 
                     q.arg = list(list(-6,5), list(-6.5,8.7)))

plot(SA.Antenne); box()
print(SA.Antenne)

#---------------------
# methode de morris
#---------------------

m2 <- morris(model=S_tild, factors=2, r=10, 
             design = list(type = "oat", 
                           levels=5, 
                           grid.jump=2), 
             binf=c(-6,-6.5), bsup=c(5,8.7))

plot(m2, xlim = c(-5, 20),ylim=c(12,22))
print(m2)


#------------------------------------------------------------------------------------#
#                          Sensitivity for 6D                                        #
#------------------------------------------------------------------------------------#


# Loading data 
load("antennes_6d_train.Rdata")

library(sensitivity)
library(DiceKriging)

# Create the kriging model using the C, response S and the kernel = matern3_2 that we got at metamodeling part
modelk = km(formula = S~., design= data.frame(C) ,response = S, covtype="matern3_2")

# Define the function that will predict for us the couvred surface Sn using our modelK, and the Xn (the data that we will genrate by later)
S_tild = function(xn, model = modelk) {
  t = predict(modelk, newdata = data.frame(xn), checkNames = FALSE, type = "SK")
  t$mean
  return (t$mean)
}

# Generate our data that follow a Uniform distribution on [LB, UB]

# longitudes  
X1 <- data.frame(matrix(runif(n=10000,min =  -6, max = 5), nrow = 10000))
X3 <- data.frame(matrix(runif(n=10000,min =  -6, max = 5), nrow = 10000))
X5 <- data.frame(matrix(runif(n=10000,min =  -6, max = 5), nrow = 10000))

#latitudes
X2 <- data.frame(matrix(runif(n=10000,min =  -6.5, max = 8.7), nrow = 10000))
X4 <- data.frame(matrix(runif(n=10000,min =  -6.5, max = 8.7), nrow = 10000))
X6 <- data.frame(matrix(runif(n=10000,min =  -6.5, max = 8.7), nrow = 10000))

X <- cbind(X1, X2, X3, X4, X5, X6)

# Geting Sn surface couvred using S_tild() that we created 
Sn = S_tild(X, modelk)

# S0 = E(S(X)) ~ mean(Sn)
S0 <- mean(Sn)

# ------------------------
# Decomposition de Sobol
# ------------------------

# Pour X1
cat("Projection sur la variable pour X1")
plot(Sn - S0 ~ X[,1], xlab=paste("X", 1, sep=""), ylab=expression(f(X)-S0[0]))

cat("Esperance conditionnelle simulee") ; 
m.ss1 <- smooth.spline(X[,1], Sn-S0)
lines(m.ss1, lwd=10, col="blue")


# Pour X2
cat("Projection sur la variable pour X2")
plot(Sn - S0 ~ X[,2], xlab=paste("X", 2, sep=""), ylab=expression(f(X)-S[0]))

cat("Esperance conditionnelle simulee") ; 
m.ss2 <- smooth.spline(X[,2], Sn-S0)
lines(m.ss2, lwd=10, col="red")

# Pour X3
cat("Projection sur la variable pour X3")
plot(Sn - S0 ~ X[,3], xlab=paste("X", 3, sep=""), ylab=expression(f(X)-S[0]))

cat("Esperance conditionnelle simulee") ; 
m.ss3 <- smooth.spline(X[,3], Sn-S0)
lines(m.ss3, lwd=10, col="blue")

# Pour X4
cat("Projection sur la variable pour X4")
plot(Sn - S0 ~ X[,4], xlab=paste("X", 4, sep=""), ylab=expression(f(X)-S[0]))

cat("Esperance conditionnelle simulee") ; 
m.ss4 <- smooth.spline(X[,4], Sn-S0)
lines(m.ss4, lwd=10, col="red")

# Pour X5
cat("Projection sur la variable pour X5")
plot(Sn - S0 ~ X[,5], xlab=paste("X", 5, sep=""), ylab=expression(f(X)-S[0]))

cat("Esperance conditionnelle simulee") ; 
m.ss5 <- smooth.spline(X[,5], Sn-S0)
lines(m.ss5, lwd=10, col="blue")

# Pour X6
cat("Projection sur la variable pour X6")
plot(Sn - S0 ~ X[,6], xlab=paste("X", 6, sep=""), ylab=expression(f(X)-S[0]))

cat("Esperance conditionnelle simulee") ; 
m.ss6 <- smooth.spline(X[,6], Sn-S0)
lines(m.ss6, lwd=10, col="red")

cat("variance de S1 est : ",var(m.ss1$y))
cat("variance de S2 est : ",var(m.ss2$y) )
cat("variance de S3 est : ",var(m.ss3$y) )
cat("variance de S4 est : ",var(m.ss4$y) )
cat("variance de S5 est : ",var(m.ss5$y) )
cat("variance de S6 est : ",var(m.ss6$y) )

cat("Indices de Sobol") ;

SA.Antenne <- fast99(model = S_tild, 
                     factors = 6, 
                     n = 1000,
                     q = "qunif", 
                     q.arg = list(list(-6,5),list(-6,5),list(-6,5), list(-6.5,8.7),list(-6.5,8.7),list(-6.5,8.7)))

plot(SA.Antenne); box()
print(SA.Antenne)


#-------------------------
# methode de morris
#-------------------------

m2 <- morris(model=S_tild, factors=6, r=10, 
             design = list(type = "oat", 
                           levels=6, 
                           grid.jump=3),
             binf=c(-6, -6.5, -6, -6.5, -6, -6.5), bsup=c(5, 8.7, 5, 8.7, 5, 8.7)
                           )

plot(m2, xlim = c(10, 20),ylim=c(10,22))
print(m2)


