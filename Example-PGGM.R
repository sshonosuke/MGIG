###----------------------------------------------------###
###           Oneshot simulation study for             ###
###         Partial Gaussian Graphical models          ###
###----------------------------------------------------###
rm(list=ls())


## load packages 
library(MASS)
library(coda)
source("PGGM-function.R")
source("MGIG-sampler.R")    # used in "PGGM-function.R"



## settings
set.seed(1)
n <- 100
q <- 3      # dimension of Y
p <- 10     # dimension of X
mc <- 12000    # length of MCMC
burn <- 2000   # length of burn-in

  
## data generation
X <- matrix(rnorm(n*p), n, p)
Inv_Om_y <- matrix(NA, q, q)
for(j in 1:q){
  for(k in 1:q){
    Inv_Om_y[j, k] <- (0.5)^(abs(j-k)) / 2
  }
}
Om_y <- solve(Inv_Om_y)
De <- matrix(0, q, p)
for(k in 1:p){
  if(runif(1)<0.5){
    De[,k] <- mvrnorm(1, rep(0, q), Om_y)
  }
}
B <- (-1)*t(De)%*%Inv_Om_y
R <- Inv_Om_y

E <- mvrnorm(n, rep(0, q), Inv_Om_y)
Y <- X%*%B + E


## MCMC (with Gibbs sampler)
tt <- proc.time()[3]
pos1 <- PGGM(Y=Y, X=X, mc=mc, burn=burn, MGIG_sampler="GS")
time1 <- proc.time()[3] - tt

## MCMC (with Riccati approximation)
tt <- proc.time()[3]
pos2 <- PGGM(Y=Y, X=X, mc=mc, burn=burn, MGIG_sampler="RC")
time2 <- proc.time()[3] - tt

## MCMC (with MH sampler)
tt <- proc.time()[3]
pos3 <- PGGM(Y=Y, X=X, mc=mc, burn=burn, MGIG_sampler="MH")
time3 <- proc.time()[3] - tt

## MCMC (with Hit-and-Run sampler)
tt <- proc.time()[3]
pos4 <- PGGM(Y=Y, X=X, mc=mc, burn=burn, MGIG_sampler="HR")
time4 <- proc.time()[3] - tt




## MSE of Delta
hDe1 <- apply(pos1$De, c(1,2), mean)
hDe2 <- apply(pos2$De, c(1,2), mean)
hDe3 <- apply(pos3$De, c(1,2), mean)
hDe4 <- apply(pos4$De, c(1,2), mean)

mean( (hDe1-De)^2 )  # Gibbs sampler
mean( (hDe2-De)^2 )  # Riccati approximation
mean( (hDe3-De)^2 )  # MH sampler
mean( (hDe4-De)^2 )  # Hit-and-Run sampler


# MSE of Omega
hOm1 <- apply(pos1$Om, c(1,2), mean)
hOm2 <- apply(pos2$Om, c(1,2), mean)
hOm3 <- apply(pos3$Om, c(1,2), mean)
hOm4 <- apply(pos4$Om, c(1,2), mean)

mean( (hOm1-Om_y)^2 )  # Gibbs sampler
mean( (hOm2-Om_y)^2 )  # Riccati approximation
mean( (hOm3-Om_y)^2 )  # MH sampler
mean( (hOm4-Om_y)^2 )  # Hit-and-Run sampler


# ESS of Omega
mean( apply(pos1$Om, c(1,2), effectiveSize) )  # Gibbs sampler
mean( apply(pos4$Om, c(1,2), effectiveSize) )  # Hit-and-Run sampler
mean( apply(pos1$Om, c(1,2), effectiveSize) ) / time1   # Gibbs sampler
mean( apply(pos4$Om, c(1,2), effectiveSize) ) / time4   # Hit-and-Run sampler

# ESS of Delta
mean( apply(pos1$De, c(1,2), effectiveSize) )  # Gibbs sampler
mean( apply(pos4$De, c(1,2), effectiveSize) )  # Hit-and-Run sampler
mean( apply(pos1$De, c(1,2), effectiveSize) ) / time1   # Gibbs sampler
mean( apply(pos4$De, c(1,2), effectiveSize) ) / time4   # Hit-and-Run sampler

