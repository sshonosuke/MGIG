library(GIGrvg)
library(control)   
library(MASS)
library(MCMCpack)

###-------------------------------------------------------###
###-------------------------------------------------------###
###            Matrix skew-t distribution                 ###
###-------------------------------------------------------###
###-------------------------------------------------------###

## IMPUT
# Y: (n,p)-data matrix 
# nu: degrees of freedom of skew-t distribution 
# mc: length of MCMC 
# bn: length of burn-in 
# MGIG_sampler: samplers for MGIG distribution (the following three options)
#   - GS: Gibbs sampler
#   - MH: Metropolis-Hastings (naive Wishart proposal)
#   - HR: Hit-and-Run sampler 
#    (Default is 'GS')

## OUTPUT (list object)
# M: MCMC samples of mean matrix 
# B: MCMC samples of skewness matrix 
# Psi: MCMC samples of (p,p)-covariance matrix
# Om: MCMC samples of (q,q)-covariance matrix
# W: MCMC samples of latent (p,p)-matrix

MVST <- function(Y, nu=10, mc=1000, burn=200, MGIG_sampler="GS", print=F, rho=NULL){
  # preparation
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  q <- dim(Y)[3]
  
  # prior for (p,q)-mean matrix (default settings)
  A_0M <- matrix(0, p, q)
  IU_0M <- (0.01)*diag(p)
  IV_0M <- (0.01)*diag(q)
  # prior for (p,q)-skewness matrix (default settings)
  A_0B <- matrix(0, p, q)
  IU_0B <- (0.01)*diag(p)
  IV_0B <- (0.01)*diag(q)
  # prior for (p,p)-covariance matrix (default settings)
  Psi0 <- diag(p)
  Inv_Psi0 <- solve(Psi0)
  eta0 <- 1
  # prior for (q,q)-covariance matrix (default settings)
  Om0 <- diag(q)
  xi0 <- 1
  
  # initial values
  W <- IW <- array(NA, dim=c(n, p, p))
  for(i in 1:n){
    IW[i,,] <- diag(runif(p, 0.8, 1.2))
  }
  M <- apply(Y, c(2,3), mean)
  B <- A_0B    # prior mean
  Psi <- Psi0    # prior mean
  Om <- Om0     # prior mean
  
  # objects to store posterior samples
  M_pos <- array(NA, dim=c(mc, p, q))
  B_pos <- array(NA, dim=c(mc, p, q))
  Psi_pos <- array(NA, dim=c(mc, p, p))
  Om_pos <- array(NA, dim=c(mc, q, q))
  W_pos <- array(NA, dim=c(mc, n, p, p))
  
  # select MGIG sampler 
  if(MGIG_sampler=="GS"){   
    rMGIG <- rMGIG_GS_single 
  }
  if(MGIG_sampler=="MH"){
    rMGIG <- function(Si_latest, la, Psi, Ga){
      N <- dim(Psi)[1]
      rMGIG_MH_single(X_latest=Si_latest, ga=la+(N+1)/2, R=Psi, TT=Ga)
    }
  }
  if(MGIG_sampler=="HR") {
    rMGIG <- function(Si_latest, la, Psi, Ga){
      rMGIG_HR_single(Si_latest=Si_latest, ka=la, Psi=Ga, Phi=Psi)
    }
  }
  
  ## MCMC iterations 
  for(itr in 1:mc) {
    Inv_Om <- solve(Om)
    
    # W (generate inverse W)
    la <- (nu + q - p - 1)/2
    Phi_tilde <- B%*%Inv_Om%*%t(B)
    for(i in 1:n){
      Ga_tilde <- Psi + (Y[i,,]-M)%*%Inv_Om%*%t(Y[i,,]-M)
      IW[i,,] <- rMGIG(Si_latest=IW[i,,], la=la, Psi=Ga_tilde, Ga=Phi_tilde)
      W[i,,] <- solve(IW[i,,]+10^(-10)*diag(p))
    }
    W_pos[itr,,,] <- W
    
    # M
    Mat <- kronecker(IV_0M, IU_0M)
    D_tilde <- solve( kronecker(Inv_Om, apply(IW, c(2,3), sum)) + Mat )
    d_tilde <- as.vector( Mat%*%as.vector(A_0M) )
    for(i in 1:n){
      d_tilde <- d_tilde + as.vector( kronecker(Inv_Om, IW[i,,])%*%as.vector(Y[i,,]-W[i,,]%*%B) )
    }
    M <- matrix(mvrnorm(1, D_tilde%*%d_tilde, D_tilde), p, q)
    M_pos[itr,,] <- M
    
    # B
    Mat <- kronecker(IV_0B, IU_0B)
    D_tilde <- solve( kronecker(Inv_Om, apply(W, c(2,3), sum)) + Mat )
    d_tilde <- as.vector( Mat%*%as.vector(A_0B) )
    for(i in 1:n){
      d_tilde <- d_tilde + as.vector( kronecker(Inv_Om, diag(p))%*%as.vector(Y[i,,]-M) )
    }
    B <- matrix(mvrnorm(1, D_tilde%*%d_tilde, D_tilde), p, q)
    B_pos[itr,,] <- B
    
    # Psi
    S_tilde <- solve( apply(IW, c(2,3), sum) + Inv_Psi0 )
    Psi <- rwish(v=eta0+n*nu, S=S_tilde)
    Psi <- Psi/Psi[1,1]
    Psi_pos[itr,,] <- Psi
    
    # Omega
    S_tilde <- Om0
    for(i in 1:n){
      resid <- Y[i,,]-M-W[i,,]%*%B
      S_tilde <- S_tilde + t(resid)%*%IW[i,,]%*%resid
    }
    Om <- riwish(v=xi0+n*p, S=S_tilde)
    Om_pos[itr,,] <- Om
  
    # print
    if(print & itr%%100==0){ print(itr) }
  }
  
  # summary 
  omit <- 1:burn
  M_pos <- M_pos[-omit,,]
  B_pos <- B_pos[-omit,,]
  Psi_pos <- Psi_pos[-omit,,]
  Om_pos <- Om_pos[-omit,,]
  W_pos <- W_pos[-omit,,,]
  
  Result <- list(M=M_pos, B=B_pos, Psi=Psi_pos, Om=Om_pos, W=W_pos)
  return( Result )
}






###-------------------------------------------------------###
###-------------------------------------------------------###
###                 Matrix t distribution                 ###
###-------------------------------------------------------###
###-------------------------------------------------------###

MVT <- function(Y, nu=10, mc=1000, burn=200, print=F){
  # preparation
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  q <- dim(Y)[3]
  
  # prior for (p,q)-mean matrix (default settings)
  A_0M <- matrix(0, p, q)
  IU_0M <- (0.01)*diag(p)
  IV_0M <- (0.01)*diag(q)
  # prior for (p,p)-covariance matrix (default settings)
  Psi0 <- diag(p)
  Inv_Psi0 <- solve(Psi0)
  eta0 <- 1
  # prior for (q,q)-covariance matrix (default settings)
  Om0 <- diag(q)
  xi0 <- 1
  
  # initial values
  W <- IW <- array(NA, dim=c(n, p, p))
  for(i in 1:n){
    W[i,,] <- IW[i,,] <- diag(p)
  }
  M <- apply(Y, c(2,3), mean)
  Psi <- Psi0    # prior mean
  Om <- Om0     # prior mean
  
  # objects to store posterior samples
  M_pos <- array(NA, dim=c(mc, p, q))
  Psi_pos <- array(NA, dim=c(mc, p, p))
  Om_pos <- array(NA, dim=c(mc, q, q))
  W_pos <- array(NA, dim=c(mc, n, p, p))
  
  ## MCMC iterations 
  for(itr in 1:mc) {
    Inv_Om <- solve(Om)
    
    # W (generate inverse W)
    for(i in 1:n){
      Ga_tilde <- Psi + (Y[i,,]-M)%*%Inv_Om%*%t(Y[i,,]-M) 
      IW[i,,] <- rwish(v=nu+q, S=solve(Ga_tilde))
      W[i,,] <- solve(IW[i,,] + 10^(-5)*diag(p))
    }
    W_pos[itr,,,] <- W
    
    # M
    Mat <- kronecker(IV_0M, IU_0M)
    D_tilde <- solve( kronecker(Inv_Om, apply(IW, c(2,3), sum)) + Mat )
    d_tilde <- as.vector( Mat%*%as.vector(A_0M) )
    for(i in 1:n){
      d_tilde <- d_tilde + as.vector( kronecker(Inv_Om, IW[i,,])%*%as.vector(Y[i,,]) )
    }
    M <- matrix(mvrnorm(1, D_tilde%*%d_tilde, D_tilde), p, q)
    M_pos[itr,,] <- M
    
    # Psi
    S_tilde <- solve( apply(IW, c(2,3), sum) + Inv_Psi0 ) + 10^(-5)*diag(p)
    Psi <- rwish(v=eta0+n*nu, S=S_tilde)
    Psi <- Psi/Psi[1,1]
    Psi_pos[itr,,] <- Psi
    
    # Omega
    S_tilde <- Om0 + 10^(-5)*diag(q)
    for(i in 1:n){
      resid <- Y[i,,]-M
      S_tilde <- S_tilde + t(resid)%*%IW[i,,]%*%resid
    }
    Om <- riwish(v=xi0+n*p, S=S_tilde)
    Om_pos[itr,,] <- Om
    
    # print
    if(print & itr%%100==0){ print(itr) }
  }
  
  # summary 
  omit <- 1:burn
  M_pos <- M_pos[-omit,,]
  Psi_pos <- Psi_pos[-omit,,]
  Om_pos <- Om_pos[-omit,,]
  W_pos <- W_pos[-omit,,,]
  
  Result <- list(M=M_pos, Psi=Psi_pos, Om=Om_pos, W=W_pos)
  return( Result )
}

