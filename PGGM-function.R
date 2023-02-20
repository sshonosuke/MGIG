###-------------------------------------------------------###
###-------------------------------------------------------###
###      MCMC for Partial Gaussian Graphical models       ###
###-------------------------------------------------------###
###-------------------------------------------------------###

# packages 
library(GIGrvg)
library(control)   



###------------------------------------------###
###             Main function                ###
###------------------------------------------###
## INPUT 
# Y: (n,q) response matrix 
# X: (n,p) covariate matrix 
# mc: length of MCMC 
# bn: length of burn-in 
# MGIG_sampler: samplers for MGIG distribution (the following four options)
#   - GS: Gibbs sampler
#   - MH: Metropolis-Hastings (naive Wishart proposal)
#   - HR: Hit-and-Run sampler 
#   - RC: Updating by the mode obtained by the Riccati approximation 
#    (Default is 'GS')

## OUTPUT (list object)
# Om: MCMC samples of 'Omega_y'
# De: MCMC samples of 'Delta'

PGGM <- function(Y, X, mc=2000, burn=500, MGIG_sampler="GS", print=F){
  # preparation
  nq <- dim(Y)
  n <- nq[1]
  q <- nq[2]
  p <- dim(X)[2]
  
  # prior (default settings)
  al <- (q + 1) / 2
  u <- q
  V <- diag(q) / u
  a <- b <- 1
  l <- rep(1, p)
  
  # useful quantities
  crossprod_Y <- crossprod(Y)
  crossprod_X <- crossprod(X)
  XXX <- as.list(rep(NA, p))
  for (i in 1:p) {
    X_i <- X[, i]
    `X_{-i}` <- X[, (-i), drop = FALSE]
    XXX[[i]] <- t(`X_{-i}`) %*% X_i
  }
  inv_V <- solve(V)
  colSums_X2 <- colSums(X^2)
  
  # initial values
  De <- matrix(0, nrow = q, ncol = p)
  z <- apply(De, 2, function(De_i){ all(De_i==0) })
  z <- as.integer(z)
  N0 <- sum(z)
  inv_Om_y <- var(Y)
  Om_y <- solve(inv_Om_y)
  la <- rep(1, p)
  pit <- 1 / 2
  
  # matrices/arrays to store posterior samples
  De_pos <- array(NA, dim=c(q, p, mc))
  Om_y_pos <- array(NA, dim=c(q, q, mc))
  la_pos <- matrix(NA, mc, p)
  pi_pos <- rep(NA, mc)
  N0_pos <- rep(NA, mc)
  
  # select MGIG sampler 
  if(MGIG_sampler=="GS"){   
    rMGIG <- rMGIG_GS_single 
  }
  if(MGIG_sampler=="MH"){
    rMGIG <- function(Si_latest, la, Psi, Ga){
      M <- dim(Psi)[1]
      rMGIG_MH_single(X_latest=Si_latest, ga=la+(M+1)/2, R=Psi, TT=Ga)
    }
  }
  if(MGIG_sampler=="HR") {
    rMGIG <- function(Si_latest, la, Psi, Ga){
      rMGIG_HR_single(Si_latest=Si_latest, ka=la, Psi=Ga, Phi=Psi)
    }
  }
  if(MGIG_sampler=="RC") {
    rMGIG <- function(Si_latest, la, Psi, Ga){
      Phi <- Psi
      Psi <- Ga
      N <- dim(Phi)[1]
      nu <- la + (N+1)/2
      alpha <- nu - (N+1)/2
      A <- alpha * diag(N)
      B <- t(chol(Phi))
      R <- 1
      eigen_Psi <- eigen(Psi)
      V <- eigen_Psi$vectors
      lam <- eigen_Psi$values
      Q <- V %*% diag(sqrt(pmax(0, lam)), nrow = N) %*% t(V)
      return( care(A, B, Q, R)$X )
    }
  }
  
  ## MCMC iterations 
  for(item in 1:mc) {
    # Delta 
    for(i in 1:p){
      la_i <- la[i]
      `De_{-i}` <- De[, (-i), drop=F]
      H_i <- Om_y %*% (t(Y) %*% X[, i]) + `De_{-i}` %*% XXX[[i]]
      s_i <- la_i / (1 + la_i * colSums_X2[i])
      p_i <- pit / (pit + ((1 - pit) / (1 + la_i * colSums_X2[i])^(q / 2)) * exp(s_i * (t(H_i) %*% inv_Om_y %*% H_i) / 2))
      z_i <- rbinom(1, size = 1, prob = p_i)
      if(z_i == 1){
        De_i <- rep(0, q)
      }else{
        De_i <- mvrnorm(1, mu=(-1)*s_i*H_i, Sigma=s_i*Om_y)
      }
      De[, i] <- De_i
      z[i] <- z_i
    }
    N0 <- sum(z)
    De_pos[,,item] <- De
    N0_pos[item] <- N0
    
    # Om_y 
    la.temp <- (n - p + N0 + u) / 2
    la.temp <- la.temp - (q + 1) / 2
    Psi.temp <- crossprod_Y + inv_V
    Ga.temp <- De %*% (crossprod_X + diag(1 / la, nrow = p)) %*% t(De)
    Om_y <- rMGIG(Si_latest=Om_y, la=la.temp, Psi=Psi.temp, Ga=Ga.temp)
    inv_Om_y <- solve(Om_y)
    Om_y_pos[,,item] <- Om_y
    
    # la 
    la_yes <- rgamma(p, shape=al, rate=l)
    Chi <- apply(De, 2, function(De_i){ t(De_i)%*%inv_Om_y%*%De_i })
    la_no <- myrgig(1, Lambda=1/2, Chi=Chi, Psi=2*l)
    la <- ifelse(test=(z==1), yes=la_yes, no=la_no)
    la_pos[item, ] <- la
    
    # non-null probability
    pit <- rbeta(1, shape1=N0+a, shape2=p-N0+b)
    pi_pos[item] <- pit
   
    # print
    if(print & item%%1000==0){ print(item) }
  }
  
  # summary 
  Om_y_pos <- Om_y_pos[,,-(1:burn), drop=F]
  De_pos <- De_pos[,,-(1:burn), drop=F]
  return( list(Om=Om_y_pos, De=De_pos) )
}


