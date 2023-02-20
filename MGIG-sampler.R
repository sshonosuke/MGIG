###-------------------------------------------------------###
###-------------------------------------------------------###
###                   MGIG samplers                       ###
###-------------------------------------------------------###
###-------------------------------------------------------###

## This code includes three MGIG samplers that can be implemented within MCMC 
# rMGIG_GS_single: Gibbs sampler
# rMGIG_MH_single: MH sampler
# rMGIG_HR_single: Hit-and-Run sampler


###    Gibbs sampler (GS)    ###
## INPUT
# Si_latest: current value of 'Sigma' whose full conditional is MGIG 
# la, Psi, Ga: parameters of target MGIG distribution

## OUTPUT
# updated 'Sigma'
rMGIG_GS_single <- function(Si_latest, la, Psi, Ga){
  p <- dim(Psi)[1]
  B_tilde <- t(chol(Si_latest))
  a <- (diag(B_tilde))^2
  B <- t(t(B_tilde) / sqrt(a))
  if(p==1){
    # a 
    inv_B <- solve(B)
    Ga_tilde <- inv_B %*% Ga %*% t(inv_B)
    Chi_a <- diag(Ga_tilde)
    Psi_a <- diag(t(B)%*%Psi%*%B)
    Lambda_a <- la + p - (1:p) + 1
    a <- myrgig(n=1, Lambda=Lambda_a, Chi=Chi_a, Psi=Psi_a)
  } 
  if(p>1){
    # a 
    inv_B <- solve(B)
    Ga_tilde <- inv_B %*% Ga %*% t(inv_B)
    Chi_a <- diag(Ga_tilde)
    Psi_a <- diag(t(B) %*% Psi %*% B)
    Lambda_a <- la + p - (1:p) + 1
    a <- myrgig(n = 1, Lambda = Lambda_a, Chi = Chi_a, Psi = Psi_a)
    # b (i=1)
    P <- crossprod(inv_B, (1 / a) * inv_B)
    M1 <- Psi
    R1 <- B
    R1[(2:p), 1] <- 0
    R1 <- t(sqrt(a) * t(R1))
    M2 <- Ga
    R2 <- t(inv_B)
    R2[1, ] <- R2[1, , drop = FALSE] + t(B[(2:p), 1]) %*% R2[(2:p), , drop = FALSE]
    R2 <- t((1 / sqrt(a)) * t(R2))
    vn <- (-1) * M1[2:p, 1:p, drop = FALSE] %*% (R1 %*% R1[1, ]) + R2[2:p, , drop = FALSE] %*% t(M2[1, , drop = FALSE] %*% R2)
    N <- a[1] * Psi[2:p, 2:p, drop = FALSE] + M2[1, 1] * P[2:p, 2:p, drop = FALSE]
    chol_N <- chol(N)
    z <- rnorm(n = p - 1, mean = 0, sd = 1) + backsolve(r = t(chol_N), x = vn, upper.tri = FALSE)
    B[2:p, 1] <- backsolve(r = chol_N, x = z, upper.tri = TRUE)
    # b (i>1)
    for (i in seq(from = 2, length.out = p - 2)) {
      M1[, i - 1] <- M1[, i - 1] + M1[, i:p, drop = FALSE] %*% B[i:p, i - 1]
      M1[i - 1, ] <- M1[i - 1, , drop = FALSE] + t(B[i:p, i - 1]) %*% M1[i:p, , drop = FALSE]
      R1[(i + 1):p, ] <- R1[(i + 1):p, , drop = FALSE] - outer(B[(i + 1):p, i], R1[i, ])
      M2[, i:p] <- M2[, i:p, drop = FALSE] - outer(M2[, i - 1], B[i:p, i - 1])
      M2[i:p, ] <- M2[i:p, , drop = FALSE] - outer(B[i:p, i - 1], M2[i - 1, ])
      R2[i, ] <- R2[i, , drop = FALSE] + t(B[(i + 1):p, i]) %*% R2[(i + 1):p, , drop = FALSE]
      vn <- (-1) * M1[(i + 1):p, 1:p, drop = FALSE] %*% (R1 %*% R1[i, ]) + R2[(i + 1):p, , drop = FALSE] %*% t(M2[i, , drop = FALSE] %*% R2)
      N <- a[i] * Psi[(i + 1):p, (i + 1):p, drop = FALSE] + M2[i, i] * P[(i + 1):p, (i + 1):p, drop = FALSE]
      chol_N <- chol(N)    
      z <- rnorm(n = p - i, mean = 0, sd = 1) + backsolve(r = t(chol_N), x = vn, upper.tri = FALSE)
      B[(i + 1):p, i] <- backsolve(r = chol_N, x = z, upper.tri = TRUE)
    }
  }
  
  # Sigma 
  t_sqrt_Si <- sqrt(a) * t(B)
  Si <- crossprod(t_sqrt_Si)
  return( Si )
}



##  MH sampler                
rMGIG_MH_single <- function(X_latest, ga, R, TT){
  # preparation
  M <- dim(R)[1]
  inv_R <- solve(R)
  df <- 2*ga
  
  # initial values
  X <- X_latest
  inv_X <- solve(X)
  
  # accept & reject 
  proposal <- rWishart(n=1, df=df, Sigma=inv_R)[,,1]
  inv_proposal <- solve(proposal)
  log_MHratio <- (-1) * (1 / 2) * sum(diag(TT %*% (inv_proposal - inv_X)))
  log_uniform <- log(runif(1))
  if(log_uniform <= log_MHratio){
    X <- proposal
  }
  # output
  return( X )
}








##   Hit-and-Run (HR) sampler           
rMGIG_HR_single <- function(Si_latest, ka, Psi, Phi){
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
  Si <- Si_latest
  eigen_Si <- eigen(Si)
  eigen_Si_vectors <- eigen_Si$vectors
  eigen_Si_values <- eigen_Si$values  
  inv_Si <- eigen_Si_vectors%*%diag(1/eigen_Si_values, nrow=p)%*%t(eigen_Si_vectors)
  eigen_log_Si_values <- log(eigen_Si_values)
  
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
  }
  
  # output
  return( Si )
}



##   univariate GIG distribution 
myrgig <- function(n, Lambda, Chi, Psi){
  LambdaChiPsi <- cbind(Lambda, Chi, Psi)   
  apply(LambdaChiPsi, 1, function(x) rgig(1, lambda=x[1], chi=x[2], psi=x[3]))
}




