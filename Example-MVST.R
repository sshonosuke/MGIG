###----------------------------------------------------###
###           Fitting matrix skew-t distribution       ###
###                to landsat satellite data           ###
###----------------------------------------------------###
rm(list=ls())

## load packages 
library(coda)
source("MVST-function.R")
source("MGIG-sampler.R")    # used in "MVST-function.R"


## settings 
sel <- 2    # class number (2, 3, 4, 5)
nu <- 5    # degrees of freedom
mc <- 2000
bn <- 500


## load dataset
data <- read.table("sat.trn", header=F)
Y_all <- array(unlist(data[,1:36]), dim=c(4435,4,9))
n <- dim(data)[1]
p <- 4
q <- 9
class <- data[,37]
Y <- Y_all[class==sel,,]
n <- dim(Y)[1]

vY <- matrix(NA, n, p*q)   # vector observation 
for(i in 1:n){
  vY[i,] <- as.vector( Y[i,,] )
}



## Fitting of matrix skew-t (Gibbs sampler)
tt <- proc.time()[3]
fit1 <- MVST(Y, nu=nu, mc=mc, burn=bn, MGIG_sampler="GS", print=T)
cpt1 <- proc.time()[3] - tt

## Fitting of matrix skew-t (MH sampler)
tt <- proc.time()[3]
fit2 <- MVST(Y, nu=nu, mc=mc, burn=bn, MGIG_sampler="MH", print=T)
cpt2 <- proc.time()[3] - tt

## Fitting of matrix skew-t (Hit-and-Run sampler)
tt <- proc.time()[3]
fit3 <- MVST(Y, nu=nu, mc=mc, burn=bn, MGIG_sampler="HR", print=T)
cpt3 <- proc.time()[3] - tt

## Fitting (matrix-t)
fit_MVT <- MVT(Y, nu=nu, mc=mc, burn=bn, print=T)



## ESS (W)
mW_pos <- apply(fit1$W, c(1,3,4), mean)
mean( apply(mW_pos, c(2,3), effectiveSize) )   # Gibbs sampler

mW_pos <- apply(fit2$W, c(1,3,4), mean)
mean( apply(mW_pos, c(2,3), effectiveSize) )   # MH sampler

mW_pos <- apply(fit3$W, c(1,3,4), mean)
mean( apply(mW_pos, c(2,3), effectiveSize) )   # Hit-and-Run sampler


## ESS (B)
mean( apply( fit1$B, c(2,3), effectiveSize) )
mean( apply( fit2$B, c(2,3), effectiveSize) )
mean( apply( fit3$B, c(2,3), effectiveSize) )


## ESS (Psi)
mean( apply( fit1$Psi, c(2,3), effectiveSize) )
mean( apply( fit2$Psi, c(2,3), effectiveSize) )
mean( apply( fit3$Psi, c(2,3), effectiveSize) )


# ESS (Omega)
mean( apply( fit1$Om, c(2,3), effectiveSize) )
mean( apply( fit2$Om, c(2,3), effectiveSize) )
mean( apply( fit3$Om, c(2,3), effectiveSize) )




## Posterior predictive loss of matrix skew-t (Gibbs sampling)
PPL <- c()
mc <- dim(fit1$M)[1]
Y_pos <- array(NA, c(mc, n, p*q))
for(k in 1:mc){
  for(i in 1:n){
    V <- kronecker(fit1$Om[k,,], fit1$W[k,i,,])
    m <- as.vector( fit1$M[k,,] + fit1$W[k,i,,]%*%fit1$B[k,,] )
    Y_pos[k,i,] <- mvrnorm(1, m, V)
  }
}
mY <- apply(Y_pos, c(2,3), mean)
V <- apply(Y_pos, c(2,3), var)
(n/(n+1))*mean((vY-mY)^2) + mean(V)    # PPL


## Posterior predictive loss of matrix-t 
Y_pos_MVT <- array(NA, c(mc, n, p*q))
for(k in 1:mc){
  m <- as.vector( fit_MVT$M[k,,] )
  for(i in 1:n){
    V <- kronecker(fit_MVT$Om[k,,], fit_MVT$W[k,i,,])
    Y_pos_MVT[k,i,] <- mvrnorm(1, m, V)
  }
}

mY_MVT <- apply(Y_pos_MVT, c(2,3), mean)
V_MVT <- apply(Y_pos_MVT, c(2,3), var)
(n/(n+1))*mean((vY-mY_MVT)^2) + mean(V_MVT)    # PPL

