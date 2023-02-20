###-----------------------------------------------------###
###     Example for generating random samples from      ###
###    MGIG distribution using four sampling methods    ###
###-----------------------------------------------------###
rm(list=ls())


# load packages
library(coda)
source("MGIG-MCMC-function.R")


# settings 
set.seed(1234)
p <- 15         # dimension 
mc <- 6000      # length of MCMC
bn <- 1000      # length of burn-in 
scenario <- 1   # scenarios of matrix parameter of MGIG distribution (from 1 to 3) 
la <- 2         # shape parameter of MGIG distribution 


# parameter values 
Ga <- diag(rep(1,p))
if(scenario==1){ Psi <- diag(rep(1,p)) }
if(scenario==2){ Psi <- diag(c(rep(1,p-2), 10, 50)) }
if(scenario==3){ Psi <- diag(1:p) }


# tuning constants 
N_FB <- dim(Psi)[1]    # tuning constant used in MH1 sampler 
rho_difference <- 5    # tuning constant used in MH1 sampler 
M_YTMG <- dim(Psi)[1]  # tuning constant used in MH2 sampler


# (proposed) Gibbs sampler
tt <- proc.time()[3]
rn1 <- rMGIG_GS(la, Psi, Ga, mc=mc, burn=bn, print=F)
cpt1 <- proc.time()[3] - tt

# MH1 sampler
tt <- proc.time()[3]
rn2 <- rMGIG_MH1(Psi=Ga, Phi=Psi, nu=la+(N_FB+1)/2, rho=rho_difference+N_FB+1, mc=mc, burn=bn, print=F)
cpt2 <- proc.time()[3] - tt

# MH2 sampler
tt <- proc.time()[3]
rn3 <- rMGIG_MH2(ga=la+(M_YTMG+1)/2, R=Psi, TT=Ga, mc=mc, burn=bn, print=F)
cpt3 <- proc.time()[3] - tt

# Hit-and-Run sampler
tt <- proc.time()[3]
rn4 <- rMGIG_HR(ka=la, Psi=Ga, Phi=Psi, mc=mc, burn=bn, print=F)
cpt4 <- proc.time()[3] - tt


# ESS 
mean( apply(rn1, c(1,2), effectiveSize) )
mean( apply(rn2, c(1,2), effectiveSize) )
mean( apply(rn3, c(1,2), effectiveSize) )
mean( apply(rn4, c(1,2), effectiveSize) )


# ESS (per second) 
mean( apply(rn1, c(1,2), effectiveSize) ) / cpt1
mean( apply(rn2, c(1,2), effectiveSize) ) / cpt2
mean( apply(rn3, c(1,2), effectiveSize) ) / cpt3
mean( apply(rn4, c(1,2), effectiveSize) ) / cpt4

