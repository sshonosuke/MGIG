###----------------------------------------------------###
###----------------------------------------------------###
###       Gibbs sampler for Matrix Generalized         ###
###       Inverse Gaussian (MGIG) distribution         ###
###----------------------------------------------------###
###----------------------------------------------------###
library(GIGrvg)
library(control)   # used in FB sampler




###---------------------------------------###
###     (proposed)   Gibbs sampler        ###
###---------------------------------------###
## INPUT 
# la: shape parameter (scalar)
# Phi: shape parameter (matrix)
# Ga: shape parameter (matrix)
# mc: length of Gibbs sampler 
# burn: burn-in period 

## OUTPUT
# random samples from MGIG

rMGIG_GS <- function(la, Psi, Ga, mc=1000, burn=mc/2, print=T){
  # function for univariate GIG distribution 
  myrgig <- function(n, Lambda, Chi, Psi){
    LambdaChiPsi <- cbind(Lambda, Chi, Psi)   
    apply(LambdaChiPsi, 1, function(x) rgig(1, lambda=x[1], chi=x[2], psi=x[3]))
  }
  
  # dimension
  p <- dim(Psi)[1]
  
  # initial values
  a <- rep(1, p)
  B <- diag(p)
  
  # arrays to store random samples
  Si_rn <- array(NA, dim=c(p, p, mc))
  
  # Gibbs sampler (when p=1)
  if(p==1){
    for (item in 1:mc){
      # a 
      inv_B <- solve(B)
      Ga_tilde <- inv_B %*% Ga %*% t(inv_B)
      Chi_a <- diag(Ga_tilde)
      Psi_a <- diag(t(B) %*% Psi %*% B)
      Lambda_a <- la + p - (1:p) + 1
      a <- myrgig(n = 1, Lambda = Lambda_a, Chi = Chi_a, Psi = Psi_a)
      # Sigma
      t_sqrt_Si <- sqrt(a) * t(B)
      Si_rn[,,item] <- crossprod(t_sqrt_Si)
    }
  }
  
  # Gibbs sampler (when p>1)
  if(p>1){
    for (item in 1:mc){
      # a 
      inv_B <- solve(B)
      Ga_tilde <- inv_B %*% Ga %*% t(inv_B)
      Chi_a <- diag(Ga_tilde)
      Psi_a <- diag(t(B) %*% Psi %*% B)
      Lambda_a <- la + p - (1:p) + 1
      a <- myrgig(n = 1, Lambda = Lambda_a, Chi = Chi_a, Psi = Psi_a)
      
      # b (i = 1)
      P <- crossprod(inv_B, (1 / a) * inv_B)
      M1 <- Psi
      R1 <- B
      R1[(2:p), 1] <- 0
      R1 <- t(sqrt(a) * t(R1))
      M2 <- Ga
      R2 <- t(inv_B)
      R2[1, ] <- R2[1, , drop=FALSE] + t(B[(2:p), 1]) %*% R2[(2:p), , drop=FALSE]
      R2 <- t((1 / sqrt(a)) * t(R2))
      vn <- (-1) * M1[2:p, 1:p, drop=FALSE] %*% (R1 %*% R1[1, ]) + R2[2:p, , drop=FALSE] %*% t(M2[1, , drop=FALSE] %*% R2)
      N <- a[1] * Psi[2:p, 2:p, drop=FALSE] + M2[1, 1] * P[2:p, 2:p, drop=FALSE]
      chol_N <- chol(N)
      z <- rnorm(n = p - 1, mean = 0, sd = 1) + backsolve(r = t(chol_N), x = vn, upper.tri=FALSE)
      B[2:p, 1] <- backsolve(r = chol_N, x = z, upper.tri=TRUE)
      #  b (i > 1)
      for (i in seq(from = 2, length.out = p - 2)) {
        M1[, i - 1] <- M1[, i - 1] + M1[, i:p, drop=FALSE] %*% B[i:p, i - 1]
        M1[i - 1, ] <- M1[i - 1, , drop=FALSE] + t(B[i:p, i - 1]) %*% M1[i:p, , drop=FALSE]
        R1[(i + 1):p, ] <- R1[(i + 1):p, , drop=FALSE] - outer(B[(i + 1):p, i], R1[i, ])
        M2[, i:p] <- M2[, i:p, drop=FALSE] - outer(M2[, i - 1], B[i:p, i - 1])
        M2[i:p, ] <- M2[i:p, , drop=FALSE] - outer(B[i:p, i - 1], M2[i - 1, ])
        R2[i, ] <- R2[i, , drop=FALSE] + t(B[(i + 1):p, i]) %*% R2[(i + 1):p, , drop=FALSE]
        vn <- (-1) * M1[(i + 1):p, 1:p, drop=FALSE] %*% (R1 %*% R1[i, ]) + R2[(i + 1):p, , drop=FALSE] %*% t(M2[i, , drop=FALSE] %*% R2)
        N <- a[i] * Psi[(i + 1):p, (i + 1):p, drop=FALSE] + M2[i, i] * P[(i + 1):p, (i + 1):p, drop=FALSE]
        chol_N <- chol(N)    
        z <- rnorm(n = p - i, mean = 0, sd = 1) + backsolve(r = t(chol_N), x = vn, upper.tri=FALSE)
        B[(i + 1):p, i] <- backsolve(r = chol_N, x = z, upper.tri=TRUE)
      }
      
      # Si 
      t_sqrt_Si <- sqrt(a) * t(B)
      Si_rn[,,item] <- crossprod(t_sqrt_Si)
      
      # print
      if(print & item%%1000==0){ print(item) }
    }
  }
    
  # output 
  return( Si_rn[,,-(1:burn), drop=FALSE] )
}





###---------------------------------------###
###             MH sampler                ###
###    (mode-adjusted Wishart proposal)   ###
###---------------------------------------###
rMGIG_MH1 <- function(Psi, Phi, nu, rho, mc=1000, burn=mc/2, print=T){
  # functions 
  ARE <- function(al, Phi, Psi){
    N <- dim(Phi)[1]
    A <- al * diag(N)
    B <- t(chol(Phi))
    R <- 1
    eigen_Psi <- eigen(Psi)
    V <- eigen_Psi$vectors
    lam <- eigen_Psi$values
    Q <- V %*% diag(sqrt(pmax(0, lam)), nrow = N) %*% t(V)
    X <- care(A, B, Q, R)
    return(X)
  }
  
  Wishart_R <- function(L, rho){
    N <- dim(L)[1]
    P <- matrix(0, N, N)
    diag(P) <- sqrt( rchisq(N, df=rho-((1:N)-1)) )
    P[upper.tri(P)] <- rnorm( N*(N-1)/2 )
    return(P%*%L)
  }
  
  weight <- function(R, Phi, Psi, inv_Si, nu, rho, log=T){
    N <- dim(R)[1]
    La <- crossprod(R)
    inv_R <- backsolve(R, diag(N), upper.tri=T)
    inv_La <- tcrossprod(inv_R)
    if(log){
      value <- sum((2*nu - rho) * log(diag(R)) - (1/2) * diag(Psi%*%inv_La + (Phi - inv_Si)%*%La))
    }else{
      value <- (prod(diag(R)))^(2*nu - rho) * exp(-(1/2) * sum(diag(Psi%*%inv_La + (Phi - inv_Si)%*%La)))
    }
    return(value)
  }

  # preparation
  N <- dim(Psi)[1]
  al <- nu - (N + 1) / 2
  La_star <- ARE(al, Phi, Psi)$X
  Si <- La_star / (rho - N - 1)
  inv_Si <- solve(Si)
  L <- chol(Si)
  
  # initial values
  R <- diag(N)
  La <- diag(N)
  
  # arrays to store random samples
  La_rn <- array(NA, dim=c(N, N, mc))
  
  # iterations (independent MH algorithm)
  for (item in 1:mc) {
    proposal_R <- Wishart_R(L=L, rho=rho)
    w1 <- weight(R=proposal_R, Phi=Phi, Psi=Psi, inv_Si=inv_Si, nu=nu, rho=rho, log=T)
    w2 <- weight(R=R, Phi=Phi, Psi=Psi, inv_Si=inv_Si, nu=nu, rho=rho, log=T)
    log_MHratio <- w1 - w2
    log_uniform <- log(runif(1))
    if(log_uniform <= log_MHratio){
      R <- proposal_R
      La <- crossprod(R)
    }
    La_rn[,,item] <- La
    
    # print  
    if(print & item%%1000==0){ print(item) }
  }
  
  # output
  return( La_rn[,,-(1:burn), drop=FALSE] )
}



###---------------------------------------###
###             MH sampler                ###
###     (naive Wishart proposal)          ###
###---------------------------------------###
rMGIG_MH2 <- function(ga, R, TT, mc=1000, burn=mc/2, print=T){
  # preparation
  M <- dim(R)[1]
  inv_R <- solve(R)
  df <- 2*ga
  
  # initial values
  X <- diag(M)
  inv_X <- solve(X)
  
  # arrays to store random samples 
  XX <- array(NA, dim=c(M, M, mc))
  
  # iterations (independent MH)
  for (item in 1:mc) {
    proposal <- rWishart(n=1, df=df, Sigma=inv_R)[,,1]
    inv_proposal <- solve(proposal)
    log_MHratio <- (-1) * (1 / 2) * sum(diag(TT %*% (inv_proposal - inv_X)))
    log_uniform <- log(runif(1))
    if(log_uniform <= log_MHratio){
      X <- proposal
      inv_X <- inv_proposal
    }
    XX[,,item] <- X
    
    # print
    if( print & item%%1000==0){ print(item) }
  }
  
  # output
  return( XX[,,-(1:burn), drop=FALSE] )
}
  



###---------------------------------------###
###        Hit-and-Run sampler            ###
###---------------------------------------###
rMGIG_HR <- function(ka, Psi, Phi, Si0=NULL, mc=1000, burn=mc/2, print=T){
  # function
  fun <- function(d){
    p <- length(d)
    d <- sort(d, decreasing=T)
    Mat_1 <- matrix(d, p, p, byrow=F)
    Mat_2 <- matrix(d, p, p, byrow=T)
    output <- sum( log((Mat_1-Mat_2)[upper.tri(Mat_1, diag=F)]) )
    return(output)
  }
  
  # preparation
  p <- dim(Psi)[1]
  L <- matrix(0, p, p)
  Lower <- lower.tri(L, diag = TRUE)
  
  # initial values
  if(is.null(Si0)){
    Si <- diag(seq(0, 2, length.out=p+1)[-1], nrow=p)
  }else{
    Si <- Si0
  }
  eigen_Si <- eigen(Si)
  eigen_Si_vectors <- eigen_Si$vectors
  eigen_Si_values <- eigen_Si$values  
  inv_Si <- eigen_Si_vectors%*%diag(1/eigen_Si_values, nrow=p)%*%t(eigen_Si_vectors)
  eigen_log_Si_values <- log(eigen_Si_values)
  
  # arrays to store random samples
  Si_rn <- array(NA, dim = c(p, p, mc))
  
  # iterations
  for (item in 1:mc) {
    l <- rnorm(p*(p+1)/2, mean=0, sd=1)
    L[Lower] <- l
    s <- sqrt( sum(l^2) )
    tilde_L <- L/s
    D <- tilde_L + t(tilde_L) - diag(diag(tilde_L), nrow=p)
    la <- rnorm(1, mean=0, sd=1)
    log_proposal <- eigen_Si_vectors%*%diag(eigen_log_Si_values, nrow=p)%*%t(eigen_Si_vectors) + la*D
    eigen_log_proposal <- eigen(log_proposal)
    eigen_log_proposal_vectors <- eigen_log_proposal$vectors
    eigen_log_proposal_values <- eigen_log_proposal$values    
    eigen_proposal_values <- exp(eigen_log_proposal_values)
    proposal <- eigen_log_proposal_vectors%*%diag(eigen_proposal_values, nrow=p)%*%t(eigen_log_proposal_vectors)
    inv_proposal <- eigen_log_proposal_vectors%*%diag(1/eigen_proposal_values, nrow=p)%*%t(eigen_log_proposal_vectors)
    log_ratio <- sum(ka*(eigen_log_proposal_values - eigen_log_Si_values) - (0.5)*diag((inv_proposal-inv_Si)%*%Psi + (proposal-Si)%*% Phi)) + sum(eigen_log_proposal_values - eigen_log_Si_values) + (fun(eigen_proposal_values) - fun(eigen_log_proposal_values)) - (fun(eigen_Si_values) - fun(eigen_log_Si_values))
    log_uniform <- log(runif(1, min=0, max=1))
    if (log_uniform <= log_ratio) {
      Si <- proposal
      eigen_Si_vectors <- eigen_log_proposal_vectors
      eigen_Si_values <- eigen_proposal_values
      inv_Si <- inv_proposal
      eigen_log_Si_values <- eigen_log_proposal_values
    }
    Si_rn[,,item] <- Si
    
    # print
    if(print & item%%1000==0){ print(item) }
  }
  
  # output
  return( Si_rn[,,-(1:burn), drop=FALSE] )
}
  




