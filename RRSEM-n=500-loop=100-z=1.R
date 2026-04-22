rm(list=ls());#50%
library(mvtnorm)
library(invgamma)
#install.packages("mvtnorm")
#install.packages("invgamma")
##defined constants##
rep_num <- 100
n <- 500
z <- 4
s <- 1
p <- 3
q <- 1
G <- 5 

n.mcmc = 20000 
n.kept = 10000 
n.burn = n.mcmc - n.kept
n.thin = 1
sig.mh = 3
sig.mh.g = 0.5

PI = 0   #for 1 latent mediator case


####### total number of acceptance for each variable in MH algorithm
##########################
# true parameter values #

Delta.true <- matrix(c(
   0, 0.2, 0, 0.5, 0.5, 0.5
), nrow = q, ncol = (1 + z + s), byrow = TRUE)

L.se.true <- cbind(Delta.true)
psd.true <- rep(0.3, q)  

B.true <- matrix(0, nrow = p, ncol = q) 
B.true[,1] <- c(1, 0.9, 0.7)

L.me.true <- B.true
psi.true <- rep(0.2, p)  

# for calculation

gam.z.true.t <- c(0.5, 0.5, 0.3, 0.5)
gam.s.true.t <- 0.5
gam.m.true.t <- 0.5
gam.true.t <- c(gam.z.true.t, gam.s.true.t, gam.m.true.t)

lam0.true.t <- rep(1, G)

gam.z.true.c <- c(0.3, 0.1, 0.1, 0.4)
gam.s.true.c <- 0.7
gam.m.true.c <- 0.1
gam.true.c <- c(gam.z.true.c, gam.s.true.c, gam.m.true.c)

lam0.true.c <- rep(1, G)

osigmaT= 0.5
omegaT = rnorm(n, 0, osigmaT); 
alphaT = 0.5;

########## Model Identification #########
### for the loading matrix of CFA model

Id.B <- matrix(c(
  0,
  1,
  1
), nrow = p, ncol = q, byrow = TRUE )

Id.B <- (Id.B > 0)
## free loadings in each row of the loading matrix B
n.b.row <- rowSums(Id.B)  
n.B <- sum(Id.B)

###for the measurment equation
Id.me <- Id.B
n.me.row <- rowSums(Id.me)  
n.me <- sum(Id.me)
 
Id.psi <- rep(1, p)
Id.psi <- as.logical(Id.psi)
n.psi <- sum(Id.psi)
 
Id.Delta <- matrix(rep(1, q*(1 + s + z)), nrow  = q)
Id.Delta <- (Id.Delta > 0)

## for the structual equation
Id.se <- cbind(Id.Delta)
n.se.row <- rowSums(Id.se)  
n.se <- sum(Id.se)

Id.psd <- rep(1, q)
Id.psd <- as.logical(Id.psd)
n.psd <- sum(Id.psd)

### for the PH model
Id.gam.t <- as.logical(rep(1, s + z + q))
n.gam.t <- sum(Id.gam.t)
Id.gam.c <- as.logical(rep(1, s + z + q))
n.gam.c <- sum(Id.gam.c)

## for the multivariate regression

##########prior setting
rho.scale <- 3

Delta0 <- matrix(rep(0, 1 + z + s), nrow = q, ncol = (1 + z + s), byrow = TRUE)

sig.delta <- rep(0.001,n.se)         # prior precision for   

L.se0 <- cbind(Delta0)
 
alpha.psd <- 9      #   0_epsilon
beta.psd <- 4       #   0_epsilon
 
B0 <- matrix(c(
  0,
  0,
  0
), nrow = p, ncol = q, byrow = TRUE)
sig.b <- 0.001                     

L.me0 <- B0

alpha.psi <- 9      #   0_  
beta.psi <- 4       #   0_  

gam0.t <- rep(0, s + z + q)
sig.gamma.t <- rep(0.001,n.gam.t)          
 
gam0.c <- rep(0, s + z + q)
sig.gamma.c <- rep(0.001,n.gam.c)     

#### For the piecewise constant   0
#### we choose G = 5, thus    = (  1,   2    3    4    5), corresponding to time interval
####                          [d1 = 0,d2],(d2,d3],(d3,d4],(d4,d5],(d5,d6 = max(T)]; 
####                          with dj = (j - 1)/5th quantile of T, vector d = (d1, d2, ..., d6), 
####                          interval length corresponding to   j is d[j+1] - d[j] 
 
alpha.lam0.t <- 1        # prior shape parameter for   j in piecewise constant baseline
beta.lam0.t <- 0.01      # prior rate parameter for   j in piecewise constant baseline

alpha.lam0.c <- 1        # prior shape parameter for   j in piecewise constant baseline
beta.lam0.c <- 0.01      # prior rate parameter for   j in piecewise constant baseline

##### some useful functions
##### calculate log-likelihood of the full conditional dist. of M|X, V, theta, & observed data
Loglike <- function(B, psi, lam0.t, gam.t, lam0.c, gam.c, omage, alphat, osigma, d, u, Delta, PI = matrix(0, ncol = q, nrow = q), psd, V, M, Z, S, X){  
  
    short <- array()
    invP  <- solve(diag(psi))
    gam.z.t <- gam.t[1:z]
    gam.s.t <- gam.t[(z + 1):(z + s)]
    gam.m.t <- gam.t[(z + s + 1):(z + s + q)]

    gam.z.c <- gam.c[1:z]
    gam.s.c <- gam.c[(z + 1):(z + s)]
    gam.m.c <- gam.c[(z + s + 1):(z + s + q)]

    for (i in 1:n) {
      short[i] <- (V[i, ] - M[i, ] %*% t(B)) %*% invP %*% t(V[i,] - M[i, ] %*% t(B))
    }
    log.like1 <- -0.5*(p*log(2*pi) + log(abs(det(diag(psi)))) + short)
    
    short1 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      part11 <- u[,j]*X[,2]*(log(lam0.t[j]) + Z %*% gam.z.t + S %*% gam.s.t  + M %*% gam.m.t + omage)
      part12 <- u[,j]*(1-X[,2])*(log(lam0.c[j]) + Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c + alphat*omage)
      
      tempa <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempa <- tempa + as.numeric(lam0.t[g]*(d[g + 1] - d[g]))
        }
      }
      part21 <- -u[,j]*(lam0.t[j]*(X[,1] - d[j]) + tempa)*exp(Z %*% gam.z.t + S %*% gam.s.t + M %*% gam.m.t + omage)
        
      tempb <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempb <- tempb + as.numeric(lam0.c[g]*(d[g + 1] - d[g]))
        }
      }
       part22 <- -u[,j]*(lam0.c[j]*(X[,1] - d[j]) + tempb)*exp(Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c + alphat*omage)
      short1[, j] <- part12 + part22 + part11 + part21

    }
    log.like2 <- apply(short1, 1, sum)
   
 ###### 
 
    const <- rep(1,n)
    Dat1 <- cbind(const, Z, S)
    Coe1 <- cbind(Delta)
    Mcen <- Dat1 %*% t(Coe1) 
    invPI0 <- solve(diag(1, q) - PI)
    short2 <- array()
    SIG <- invPI0 %*% diag(psd,q) %*% t(invPI0)
    ISIG <- solve(SIG)
    for (i in 1:n) {
      short2[i] <- t(M[i, ] - invPI0 %*% Mcen[i, ]) %*% ISIG %*% (M[i, ] - invPI0 %*% Mcen[i, ])
    }
    log.like3 <- -0.5*(q*log(2*pi) + log(abs(det(SIG))) + short2)
    
    log.like <- log.like1 + log.like2 + log.like3 
    return(log.like)
}
#Loglike(B, psi, lam0.t, gam.t, lam0.c, gam.c, omage, alphat, d, u, Delta,PI, psd, V, M, Z, S, X)

Loglike.omage <- function(B, psi, lam0.t, gam.t, lam0.c, gam.c, omage, alphat, osigma, d, u, Delta, PI = matrix(0, ncol = q, nrow = q), psd, V, M, Z, S, X){
  
    gam.z.t <- gam.t[1:z]
    gam.s.t <- gam.t[(z + 1):(z + s)]
    gam.m.t <- gam.t[(z + s + 1):(z + s + q)]

    gam.z.c <- gam.c[1:z]
    gam.s.c <- gam.c[(z + 1):(z + s)]
    gam.m.c <- gam.c[(z + s + 1):(z + s + q)]

    short1 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      part11 <- u[,j]*X[,2]*(log(lam0.t[j]) + Z %*% gam.z.t + S %*% gam.s.t + M %*% gam.m.t + omage)
      part12 <- u[,j]*(1-X[,2])*(log(lam0.c[j]) + Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c + alphat*omage)
      
      tempa <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempa <- tempa + as.numeric(lam0.t[g]*(d[g + 1] - d[g]))
        }
      }
      part21 <- -u[,j]*(lam0.t[j]*(X[,1] - d[j]) + tempa)*exp(Z %*% gam.z.t + S %*% gam.s.t + M %*% gam.m.t + omage)
        
      tempb <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempb <- tempb + as.numeric(lam0.c[g]*(d[g + 1] - d[g]))
        }
      }
       part22 <- -u[,j]*(lam0.c[j]*(X[,1] - d[j]) + tempb)*exp(Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c + alphat*omage)
      short1[, j] <- part12 + part22 + part11 + part21

    }
    log.like2 <- apply(short1, 1, sum)
###########    
    log.like4 <- sum(-(t(omage)%*%omage)/(2*osigma))        
    log.like <- log.like2 + log.like4 
    return(log.like)
} 
### calculate the loglikelihood of full conditional dist. of gam|X, theta, & observed data

Loglike_GA.t <- function(lam0.t, gam.t, lam0.c, gam.c, omage, alphat,osigma, d, u, M, Z, S, X){

  gam.z.t <- gam.t[1:z]
  gam.s.t <- gam.t[(z + 1):(z + s)]
  gam.m.t <- gam.t[(z + s + 1):(z + s + q)]

  sum1 <- 0
  for (j in 1:G) {
    part11 <- u[,j]*X[,2]*(log(lam0.t[j]) + Z %*% gam.z.t + S %*% gam.s.t + M %*% gam.m.t + omage)
       
    tempe <- 0
    if (j > 1) {
      for (g in 1:(j - 1)) {
        tempe <- tempe + as.numeric(lam0.t[g]*(d[g + 1] - d[g]))
      }
    }

    part21 <- -u[,j]*(lam0.t[j]*(X[,1] - d[j]) + tempe)*exp(Z %*% gam.z.t + S %*% gam.s.t + M %*% gam.m.t+ omage)
   
    sum1 <- sum1 + sum(part11) + sum(part21)  
  } 

  log.like1 <- sum1

#######
 
  invP.t  <- diag(sig.gamma.t, n.gam.t)
  gam.t <- matrix(gam.t, ncol = 1)
  short2.t <- t(gam.t - gam0.t) %*% invP.t %*% (gam.t - gam0.t)
  log.like2.t <- -0.5*((n.gam.t)*log(2*pi) + log(abs(det(chol2inv(chol(invP.t))))) + short2.t)

  log.like <- log.like1 + log.like2.t
  return(log.like)
}

Loglike_GA.c <- function(lam0.t, gam.t, lam0.c, gam.c, omage, alphat, osigma,d, u, M, Z, S, X){

  gam.z.c <- gam.c[1:z]
  gam.s.c <- gam.c[(z + 1):(z + s)]
  gam.m.c <- gam.c[(z + s + 1):(z + s + q)]

  sum1 <- 0
  for (j in 1:G) {
    part12 <- u[,j]*(1-X[,2])*(log(lam0.c[j]) + Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c +alphat*omage)

     tempf <- 0
    if (j > 1) {
      for (g in 1:(j - 1)) {
        tempf <- tempf + as.numeric(lam0.c[g]*(d[g + 1] - d[g]))
      }
    }

    part22 <- -u[,j]*(lam0.c[j]*(X[,1] - d[j]) + tempf)*exp(Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c +alphat*omage)
 
    sum1 <- sum1 + sum(part12) + sum(part22) 
  } 

  log.like1 <- sum1
 
 
  invP.c  <- diag(sig.gamma.c, n.gam.c)
  gam.c <- matrix(gam.c, ncol = 1)
  short2.c <- t(gam.c - gam0.c) %*% invP.c %*% (gam.c - gam0.c)
  log.like2.c <- -0.5*((n.gam.c)*log(2*pi) + log(abs(det(chol2inv(chol(invP.c))))) + short2.c)

  log.like <- log.like1 + log.like2.c
  return(log.like)
}

Loglike_alphat <- function(lam0.t, gam.t, lam0.c, gam.c, omage, alphat, osigma,d, u, M, Z, S, X){
  sum1 <- 0
  for (j in 1:G) {
    part12 <- u[,j]*(1-X[,2])*(log(lam0.c[j]) + Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c +alphat*omage)

     tempf <- 0
    if (j > 1) {
      for (g in 1:(j - 1)) {
        tempf <- tempf + as.numeric(lam0.c[g]*(d[g + 1] - d[g]))
      }
    }

    part22 <- -u[,j]*(lam0.c[j]*(X[,1] - d[j]) + tempf)*exp(Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c +alphat*omage)
 
    sum1 <- sum1 + sum(part12) + sum(part22) 
  } 

  log.like1 <- sum1
   
  log.like2.c <- dnorm(alphat, 0, 0.5, log=TRUE)

  log.like <- log.like1 + log.like2.c
  return(log.like)
}
 
  
####### Store the mean and sd. of each parameter in all replications
mean.B <- sd.B <- array(0, dim = c(rep_num, n.B))
mean.psi <- sd.psi <- array(0, dim = c(rep_num, n.psi))
mean.L.se <- sd.L.se <- array(0, dim = c(rep_num, n.se))
mean.psd <- sd.psd <- array(0, dim = c(rep_num, n.psd))
mean.gam.t <- sd.gam.t <- array(0, dim = c(rep_num, n.gam.t))
mean.lam0.t <- sd.lam0.t <- array(0, dim = c(rep_num, G))
mean.gam.c <- sd.gam.c <- array(0, dim = c(rep_num, n.gam.c))
mean.lam0.c <- sd.lam0.c <- array(0, dim = c(rep_num, G))
mean.alphat <- sd.alphat <- array(0, dim = c(rep_num, 1))
mean.omage <- sd.omage <- matrix(NA, nrow = rep_num, ncol = n)

####### Store the quantile of each parameter in all replications
q.B <- array(0, dim = c(rep_num, n.B*2))
q.psi <- array(0, dim = c(rep_num, n.psi*2))
q.L.se <- array(0, dim = c(rep_num, n.se*2))
q.psd <- array(0, dim = c(rep_num, n.psd*2))
q.gam.t <- array(0, dim = c(rep_num, n.gam.t*2))
q.lam0.t <- array(0, dim = c(rep_num, G*2))
q.gam.c <- array(0, dim = c(rep_num, n.gam.c*2))
q.lam0.c <- array(0, dim = c(rep_num, G*2))
q.alphat <- array(0, dim = c(rep_num, 1*2))
q.omage <- array(NA, dim = c(rep_num, n, 2))

##############################################################################
 
##### for storing results
result <- list('B' = array(NA, dim = c(n.mcmc, n.B)),
               'alphat' = array(NA, dim = c(n.mcmc, 1)),
               'omage' = matrix(NA, nrow = n.mcmc, ncol = n),
               'L.se' = array(NA, dim = c(n.mcmc, n.se)),
               'gam.t' = array(NA, dim = c(n.mcmc, n.gam.t)),
               'gam.c' = array(NA, dim = c(n.mcmc, n.gam.c)),
               'psd' = array(NA, dim = c(n.mcmc, q)),
               'psi' = array(NA, dim = c(n.mcmc, n.psi)), 
               'lam0.t' = array(NA, dim = c(n.mcmc, G)),
               'lam0.c' = array(NA, dim = c(n.mcmc, G)))

t0 <- Sys.time()
cat('Generating', rep_num, 'sets of data...  Please wait \n')
#############
#simulate data with censoring
 
for (crep in 1:rep_num) {
    set.seed(201 + crep *10) 
  ################################
  ######  data generation   ######
S <- matrix(rbinom(n, 1, 0.2), ncol = 1)
  #S <- S[, 2:3]
  Z <- rmvnorm(n, mean = rep(0,z), sigma = diag(1,z))
 
  X <- array(NA, dim = c(n, 2))
 
  Dat <- cbind(rep(1,n), Z, S)
  Coe <- cbind(Delta.true)
  M <- Dat %*% t(Coe) + rmvnorm(n, mean = rep(0, q), sigma = diag(psd.true, q))
  M.true <- M <- M %*% t(solve((diag(1, q) - PI)))
  const <- rep(1,n)
  Dat1 <- cbind(M)
  V <- Dat1 %*% t(L.me.true) + rmvnorm(n, mean = rep(0, p), sigma = diag(psi.true))
  
  # censoring time and survival time
  Dat2 <- cbind(Z, S, M)
  temp.t <- exp(Dat2 %*% gam.true.t + omegaT)
  temp.c <- exp(Dat2 %*% gam.true.c + alphaT*omegaT)
  # with lambda_0(t) = 1
  T.t <- (-1/(1*temp.t))*log(runif(n))
  T.c <- (-1/(1*temp.c))*log(runif(n))
  # with lamda_0(t) = 2t + 1
  #T.t <- sqrt((-1/temp.t)*log(runif(n)) + 1/4) - 1/2
  #T.c <- sqrt((-1/temp.c)*log(runif(n)) + 1/4) - 1/2
  
  X[, 1] <- (T.t <= T.c)*T.t + (T.t > T.c)*T.c
  X[, 2] <- ifelse(T.t <= T.c, 1, 0)
 
   m0<-length(which(X[, 2]==0))
   m0/(n) 
### generate initial values for latent mediator M1
###init value 1 for parameters

B <- matrix(c(
  1,
  0,
  0
), nrow = p, ncol = q, byrow = TRUE)

L.me <- B

psi <- rep(1, p)  
Delta <- matrix(rep(0, q*(1 + s + z)), nrow = q, ncol = (1 + z + s), byrow = TRUE)

L.se <- cbind(Delta)

psd <- rep(1, q)  

gam.z.t <- c(0, 0, 0, 0)
gam.s.t <- 0
gam.m.t <- 0
gam.t <- c(gam.z.t, gam.s.t, gam.m.t)

lam0.t <- rep(1, G)


gam.z.c <- c(0, 0, 0, 0)
gam.s.c <- 0
gam.m.c <- 0
gam.c <- c(gam.z.c, gam.s.c, gam.m.c)

lam0.c <- rep(1, G)

alphat = 0.5
osigma = 0.5
omage = omegaT

# to facilitate computation later 
iv.psi <- 1/psi
iv.sqrt.psi <- sqrt(iv.psi)

iv.psd <- 1/psd
iv.sqrt.psd <- sqrt(iv.psd)

if (q > 0) {
  const <- rep(1,n)
  Dat <- cbind(const, Z, S)
  Coe <- cbind(Delta)
  M <- Dat %*% t(Coe) + rmvnorm(n, mean = rep(0, q), sigma = diag(psd, q))
  M <- M %*% t(solve((diag(1, q) - PI)))   
}
  
  
  d <- c(quantile(X[ , 1], probs = seq(0, 1, 1/G)))
  d[1] <- 0
  u <- array(0,dim = c(n, G))
  for (j in 1:G) {
    al <- as.logical(d[j] < X[,1])*(X[,1] <= d[j + 1])
    u[, j] <- al
  }

n.accept.gam.t = 1
n.accept.M = rep(1,n) 
n.accept.gam.c = 1
n.accept.omage = rep(1,n) 
n.accept.alphat = 1

   ## iteration
   it=1
   while (it<n.mcmc+1)
{
cat("loopp", crep, "Iteration", it, fill=TRUE);  
 

  # step1: update latent mediators M
  if (q > 0) {
    PI0 <- (diag(1, q) - PI)
    ISG <- crossprod(iv.sqrt.psi * B) + crossprod(iv.sqrt.psd * PI0)     #   _  ^-1
    SIG <- chol2inv(chol(ISG))                     
    SIG <- sig.mh*SIG
    cSIG <- chol(SIG)
  }
  #### Random Walk Metropolis, M_t + N[0,   _mh^2*  _  _t]
  M.new <- M + t(crossprod(cSIG, matrix(rnorm(q*n), nrow = q)))
  ll1 <- Loglike(B, psi, lam0.t, gam.t, lam0.c, gam.c, omage, alphat, osigma,d, u, Delta,PI, psd, V, M, Z, S, X)
  ll2 <- Loglike(B, psi, lam0.t, gam.t, lam0.c, gam.c, omage, alphat,osigma, d, u, Delta, PI, psd, V, M.new, Z, S, X)
  ## acceptance ratio
  p.accept <- exp(ll2 - ll1)
  accept <- (runif(n) < p.accept)  
  M[accept, ] <- M.new[accept, ]
  n.accept.M <- n.accept.M + accept
  #n.accept.M/n.mcmc
  # step2: update loading matrix B and residual variance psi
  count.B <- 1
  L.me <- B
  omg.me <- M
  for (k in 1:p) {
    free <- Id.me[k, ]                  # which column on row K is free/fixed
    len <- n.me.row[k]
    Vcen <- t(L.me[k, !free, drop = F] %*% t(omg.me)[!free, , drop = F])
    Mk <- t(omg.me)[free, , drop = F]   # data of the latent variable that the kth indicator loads on
    Psiginv <- rep(sig.b, len)          # H_0yk


    Vk.star <- V[ , k] - Vcen
    alpha.psi.star <- alpha.psi + 0.5*n
    beta.psi.star <- beta.psi + 0.5*sum(Vk.star^2)
    if (len > 0) {
      A_vk <- chol2inv(chol(diag(Psiginv, len) + tcrossprod(Mk)))
      temp <- Psiginv*L.me0[k, free] + Mk %*% Vk.star
      a_vk <- A_vk %*% temp
      beta.psi.star <- beta.psi.star + 0.5*(sum(L.me0[k, free]*Psiginv*L.me0[k, free]) - sum(temp*a_vk))
    }

    iv.psi[k] <- rgamma(1, shape = alpha.psi.star, rate = beta.psi.star)
    psi[k] <- 1/iv.psi[k]                  
    iv.sqrt.psi[k] <- sqrt(iv.psi[k])

    if (len > 0) {
      L.me[k, free] <- rmvnorm(1, a_vk, sigma = psi[k]*A_vk)
      if (n.b.row[k] > 0) {
        B[k, ] <- L.me[k, ]
      }
 
        if (n.b.row[k] > 0) {
          result$B[it, count.B:(count.B + n.b.row[k] - 1)] <- B[k, Id.B[k, ]]
        }
 
      count.B <- count.B + n.b.row[k]
    }  
  }

  # step3: update regression coefficients delta and mean_m    and error variance psd
  ## Note!! Need revise if the number of latent mediator changed!
  count.se <- 1
  L.se <- cbind(Delta)
  omg.se <- cbind(rep(1,n), Z, S)
  for (k in 1:q) {
    free <- Id.se[k, ]  # which column on row K is free/fixed
    len <- n.se.row[k]
    Mcen <- t(L.se[k, !free, drop = F] %*% t(omg.se)[!free, , drop = F])
    Mk.star <- M[ , k] - as.vector(Mcen)
    alpha.psd.star <- alpha.psd + 0.5*n
    beta.psd.star <- beta.psd + 0.5*sum(Mk.star^2)

    if (len > 0) {
      Yk <- omg.se[, free, drop = F]  
      iH0dk <- diag(sig.delta, len)
      L.se0k <- L.se0[k, free]
      A_dk <- chol2inv(chol(iH0dk + crossprod(Yk)))
      temp <- iH0dk %*% L.se0k + t(Yk) %*% Mk.star
      a_dk <- A_dk %*% temp
      beta.psd.star <- beta.psd.star + 0.5*(crossprod(crossprod(iH0dk, L.se0k), L.se0k)
                                    - crossprod(a_dk, temp))
    }

    iv.psd[k] <- rgamma(1, shape = alpha.psd.star, rate = beta.psd.star)
    psd[k] <- 1/iv.psd[k]
    iv.sqrt.psd[k] <- sqrt(iv.psd[k])

    if (len > 0) {
      L.se[k, free] <- rmvnorm(1, a_dk, sigma = psd[k]*A_dk)

        result$L.se[it, count.se:(count.se + len - 1)] <- L.se[k, free]

      count.se <- count.se + len
    }
  }

  Delta <- L.se[, 1:(1 + z + s), drop = F]
 
    result$psd[it, ] <- psd
    result$psi[it, ] <- psi
 

  # step4: update regression coefficients gam, M-H I guess
  #### Random Walk Metropolis, gam_t + N[0,   _mh^2*I}

 gam.new.t <- gam.t + rmvnorm(1, rep(0, n.gam.t), sigma = diag(rep(0.0015, n.gam.t)))
 llgt1 <- Loglike_GA.t(lam0.t, gam.t, lam0.c, gam.c, omage, alphat,osigma, d, u, M, Z, S, X)
 llgt2 <- Loglike_GA.t(lam0.t,  gam.new.t, lam0.c, gam.c, omage, alphat,osigma, d, u, M, Z, S, X)
 ## acceptance ratio
 p.accept <- exp(llgt2 - llgt1)
 accept <- (runif(1) < p.accept)
 if (accept) gam.t <-  gam.new.t
 gam.z.t <- gam.t[1:z]
 gam.s.t <- gam.t[(z + 1):(z + s)]
 gam.m.t <- gam.t[(z + s + 1):(z + s + q)]
 n.accept.gam.t <- n.accept.gam.t + accept
 #n.accept.gam.t/n.mcmc
    result$gam.t[it, ] <- gam.t
 
  # step4: update regression coefficients gam, M-H I guess
  #### Random Walk Metropolis, gam_t + N[0,   _mh^2*I}

 gam.new.c <- gam.c + rmvnorm(1, rep(0, n.gam.c), sigma = diag(rep(0.0025, n.gam.c)))
 llgc1 <- Loglike_GA.c(lam0.t, gam.t, lam0.c, gam.c, omage, alphat, osigma,d, u, M, Z, S, X)
 llgc2 <- Loglike_GA.c(lam0.t,  gam.t, lam0.c, gam.new.c, omage, alphat, osigma,d, u, M, Z, S, X)
 ## acceptance ratio
 p.accept <- exp(llgc2 - llgc1)
 accept <- (runif(1) < p.accept)
 if (accept) gam.c <-  gam.new.c
 gam.z.c <- gam.c[1:z]
 gam.s.c <- gam.c[(z + 1):(z + s)]
 gam.m.c <- gam.c[(z + s + 1):(z + s + q)]
 n.accept.gam.c <- n.accept.gam.c + accept
 #n.accept.gam.c/n.mcmc
    result$gam.c[it, ] <- gam.c 
 
  # step5: update piecewise constant hazard lam0, conjugate gamma prior?
   
  alpha.lam0.star.t <- rep(alpha.lam0.t, G) + as.vector(t(X[,2]) %*% u)

  temp2.t <- array(0, dim = c(n, G))

  for (j in 1:G) {
    temp3.t <- 0
        if (j < G) {
          for (g in (j + 1):G) {
            temp3.t <- temp3.t + u[, g]*(d[j + 1] - d[j])
          }
        }
    temp2.t[,j] <- exp(cbind(Z, S, M) %*% as.vector(gam.t)+omage)*(u[,j]*(X[,1] - d[j]) + temp3.t)
  }

  beta.lam0.star.t <- rep(beta.lam0.t, G) + colSums(temp2.t)

  for (j in 1:G) {
    lam0.t[j] <- rgamma(1, shape = alpha.lam0.star.t[j], rate = beta.lam0.star.t[j])
  }
 
     result$lam0.t[it, ] <- lam0.t
 

  # step5: update piecewise constant hazard lam0, conjugate gamma prior?
   
  alpha.lam0.star.c <- rep(alpha.lam0.c, G) + as.vector(t(1-X[,2]) %*% u)

  temp2.c <- array(0, dim = c(n, G))

  for (j in 1:G) {
    temp3.c <- 0
        if (j < G) {
          for (g in (j + 1):G) {
            temp3.c <- temp3.c + u[, g]*(d[j + 1] - d[j])
          }
        }
    temp2.c[,j] <- exp(cbind(Z, S, M) %*% as.vector(gam.c)+alphat*omage)*(u[,j]*(X[,1] - d[j]) + temp3.c)
  }

  beta.lam0.star.c <- rep(beta.lam0.c, G) + colSums(temp2.c)

  for (j in 1:G) {
    lam0.c[j] <- rgamma(1, shape = alpha.lam0.star.c[j], rate = beta.lam0.star.c[j])
  }
 
     result$lam0.c[it, ] <- lam0.c

# step7

 alphat.new <- alphat + rnorm(1, 0, 0.25)
 lla1 <- Loglike_alphat(lam0.t, gam.t, lam0.c, gam.c, omage, alphat, osigma,d, u, M, Z, S, X)
 lla2 <- Loglike_alphat(lam0.t, gam.t, lam0.c, gam.c, omage, alphat.new, osigma,d, u, M, Z, S, X)
 ## acceptance ratio
 a.accept <- exp(lla2 - lla1)
 accept <- (runif(1) < a.accept)
 if (accept) alphat <-  alphat.new
  alphat[accept] <- alphat.new[accept]  
 n.accept.alphat <- n.accept.alphat + accept
  result$alphat[it] <- alphat
 #n.accept.alphat/n.mcmc 
 #alphaT 
# step8

  #### Random Walk Metropolis, M_t + N[0,   _mh^2*  _  _t]
  omage.new <- omage + rnorm(n, 0, 0.045)
  lo1 <- Loglike.omage(B, psi, lam0.t, gam.t, lam0.c, gam.c, omage, alphat, osigma,d, u, Delta,PI, psd, V, M, Z, S, X)
  lo2 <- Loglike.omage(B, psi, lam0.t, gam.t, lam0.c, gam.c, omage.new, alphat, osigma,d, u, Delta,PI, psd, V, M, Z, S, X)
  ## acceptance ratio
  o.accept <- exp(lo2 - lo1)
  accept <- (runif(n) < o.accept)
  omage[accept] <- omage.new[accept]  
  n.accept.omage <- n.accept.omage + accept
  result$omage[it,] <- omage
  #n.accept.omage/n.mcmc 
 #omegaT -  omage
#
library(invgamma)
#actovr
aa=0.01
bb=0.01
osigma <- rinvgamma(1, aa+n/2, bb+0.5*t(omage)%*%(omage))
#omage-omegaT

  it=it+1
} # end of mcmc
 
  mean.B[crep, ] <- colMeans(result$B[n.burn:n.mcmc,])
  mean.psi[crep, ] <- colMeans(result$psi[n.burn:n.mcmc,])
  mean.L.se[crep, ] <- colMeans(result$L.se[n.burn:n.mcmc,])
  mean.psd[crep] <- mean(result$psd[n.burn:n.mcmc])
  mean.gam.t[crep, ] <- colMeans(result$gam.t[n.burn:n.mcmc,])
  mean.lam0.t[crep, ] <- colMeans(result$lam0.t[n.burn:n.mcmc,])
  mean.gam.c[crep, ] <- colMeans(result$gam.c[n.burn:n.mcmc,])
  mean.lam0.c[crep, ] <- colMeans(result$lam0.c[n.burn:n.mcmc,])
  mean.omage[crep,] <- colMeans(result$omage[n.burn:n.mcmc,])
  mean.alphat[crep] <- mean(result$alphat[n.burn:n.mcmc])
   
  sd.B[crep, ] <- apply(result$B[n.burn:n.mcmc,], 2, sd)
  sd.psi[crep, ] <- apply(result$psi[n.burn:n.mcmc,],2, sd)
  sd.L.se[crep, ] <- apply(result$L.se[n.burn:n.mcmc,], 2, sd)
  sd.psd[crep, ] <- sd(result$psd[n.burn:n.mcmc])
  sd.gam.t[crep, ] <- apply(result$gam.t[n.burn:n.mcmc,], 2, sd)
  sd.lam0.t[crep, ] <- apply(result$lam0.t[n.burn:n.mcmc,], 2, sd)
  sd.gam.c[crep, ] <- apply(result$gam.c[n.burn:n.mcmc,], 2, sd)
  sd.lam0.c[crep, ] <- apply(result$lam0.c[n.burn:n.mcmc,], 2, sd)

  
  q.B[crep, ] <- as.vector(apply(result$B[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.psi[crep, ] <- as.vector(apply(result$psi[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.L.se[crep, ] <- as.vector(apply(result$L.se[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.psd[crep,] <- quantile(result$psd[n.burn:n.mcmc], c(0.025, 0.975), na.rm = TRUE)
  q.gam.t[crep, ] <- as.vector(apply(result$gam.t[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.lam0.t[crep, ] <- as.vector(apply(result$lam0.t[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.gam.c[crep, ] <- as.vector(apply(result$gam.c[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.lam0.c[crep, ] <- as.vector(apply(result$lam0.c[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))


}
 

B.list <- list(Mean = colMeans(mean.B), Sd = colMeans(sd.B),Se = apply(mean.B,2,sd),
               Bias = colMeans(mean.B) - B.true[Id.B],
               Rms = sqrt(rowSums((t(mean.B) - B.true[Id.B])^2) / rep_num),
               p.CI95 = colMeans((q.B[ , (1:n.B)*2 - 1] - t(replicate(
                 rep_num, B.true[Id.B])) > 0) * (q.B[, (1:n.B)*2] - t(
                   replicate(rep_num, B.true[Id.B])) > 0) == 0))

psi.list <- list(Mean = colMeans(mean.psi), Sd = colMeans(sd.psi),Se = apply(mean.psi,2,sd),
                 Bias = colMeans(mean.psi) - psi.true[Id.psi],
                 Rms = sqrt(rowSums((t(mean.psi) - psi.true[Id.psi])^2) / rep_num),
                 p.CI95 = colMeans((q.psi[ , (1:n.psi)*2 - 1] - t(replicate(
                   rep_num, psi.true[Id.psi])) > 0) * (q.psi[, (1:n.psi)*2] - t(
                     replicate(rep_num, psi.true[Id.psi])) > 0) == 0))

L.se.list <- list(Mean = colMeans(mean.L.se), Sd = colMeans(sd.L.se),Se = apply(mean.L.se,2,sd),
               Bias = colMeans(mean.L.se) - t(L.se.true)[t(Id.se)],
               Rms = sqrt(rowSums((t(mean.L.se) - t(L.se.true)[t(Id.se)])^2) / rep_num),
               p.CI95 = colMeans((q.L.se[ , (1:n.se)*2 - 1] - t(replicate(
                 rep_num, t(L.se.true)[t(Id.se)])) > 0) * (q.L.se[, (1:n.se)*2] - t(
                   replicate(rep_num, t(L.se.true)[t(Id.se)])) > 0) == 0))

psd.list <- list(Mean = colMeans(mean.psd), Sd = colMeans(sd.psd),Se = apply(mean.psd,2,sd),
                 Bias = colMeans(mean.psd) - psd.true[Id.psd],
                 Rms = sqrt(rowSums((t(mean.psd) - psd.true[Id.psd])^2) / rep_num),
                 p.CI95 = rowMeans((q.psd[ , (1:n.psd)*2 - 1] - t(replicate(
                   rep_num, psd.true[Id.psd])) > 0) * (q.psd[, (1:n.psd)*2] - t(
                     replicate(rep_num, psd.true[Id.psd])) > 0) == 0))

gam.t.list <- list(Mean = colMeans(mean.gam.t), Sd = colMeans(sd.gam.t),Se = apply(mean.gam.t,2,sd),
                 Bias = colMeans(mean.gam.t) - gam.true.t,
                 Rms = sqrt(rowSums((t(mean.gam.t) - gam.true.t)^2) / rep_num),
                 p.CI95 = colMeans((q.gam.t[ , (1:n.gam.t)*2 - 1] - t(replicate(
                   rep_num, gam.true.t)) > 0) * (q.gam.t[, (1:n.gam.t)*2] - t(
                     replicate(rep_num, gam.true.t)) > 0) == 0))

gam.c.list <- list(Mean = colMeans(mean.gam.c), Sd = colMeans(sd.gam.c),Se = apply(mean.gam.c,2,sd),
                 Bias = colMeans(mean.gam.c) - gam.true.c,
                 Rms = sqrt(rowSums((t(mean.gam.c) - gam.true.c)^2) / rep_num),
                 p.CI95 = colMeans((q.gam.c[ , (1:n.gam.c)*2 - 1] - t(replicate(
                   rep_num, gam.true.c)) > 0) * (q.gam.c[, (1:n.gam.c)*2] - t(
                     replicate(rep_num, gam.true.c)) > 0) == 0))

lam0.t.list <- list(Mean = colMeans(mean.lam0.t), Sd = colMeans(sd.lam0.t),Se = apply(mean.lam0.t,2,sd),
                  Bias = colMeans(mean.lam0.t) - lam0.true.t,
                  Rms = sqrt(rowSums((t(mean.lam0.t) - lam0.true.t)^2) / rep_num),
                  p.CI95 = colMeans((q.lam0.t[ , (1:G)*2 - 1] - t(replicate(
                    rep_num, lam0.true.t)) > 0) * (q.lam0.t[, (1:G)*2] - t(
                      replicate(rep_num, lam0.true.t)) > 0) == 0))

lam0.c.list <- list(Mean = colMeans(mean.lam0.c), Sd = colMeans(sd.lam0.c),Se = apply(mean.lam0.c,2,sd),
                  Bias = colMeans(mean.lam0.c) - lam0.true.c,
                  Rms = sqrt(rowSums((t(mean.lam0.c) - lam0.true.c)^2) / rep_num),
                  p.CI95 = colMeans((q.lam0.c[ , (1:G)*2 - 1] - t(replicate(
                    rep_num, lam0.true.c)) > 0) * (q.lam0.c[, (1:G)*2] - t(
                      replicate(rep_num, lam0.true.c)) > 0) == 0))


est.list <- list(B = B.list, psi = psi.list, L.se = L.se.list,
                 psd = psd.list, gam.t = gam.t.list, lam0.t = lam0.t.list, gam.c = gam.c.list, lam0.c = lam0.c.list)

names(est.list)[3] <- "Beta"
for (i in 1:length(est.list)) {
  write.csv(est.list[[i]], paste(names(est.list)[i],".csv",sep = ""))
}

Sys.time()-t0; 
save.image(paste(crep, '.Rdata', sep = ''))

 
