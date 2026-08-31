##### MODEL FUNCTIONS FOR BIM-MC #####
##### Author: William Hayes      #####

library(cubature) # for numerical integration of double integrals
library(mvtnorm)  # for multivariate normal random number generation


# Model log-likelihood
# Expects discrete confidence ratings (JOLs) on an n-point scale: 1, 2, . . ., n
# Memory outcomes should be 1 = correct retrieval, 0 = incorrect retrieval (lure), or -1 = no memory
#   data:   N x 2 matrix with confidence ratings in column 1 and memory outcomes in column 2
#   k:      number of discrete categories on the confidence / JOL scale
#   n_alts: number of choice alternatives in the multiple-choice memory test
BIM_MC_loglik <- function(params, data, k, n_alts) {
  P_exp <- params[1]
  mu_m <- params[2]
  M_conf <- params[3]
  rho <- params[4]
  crit <- params[5]
  
  # enforce parameter bounds
  if (
    (P_exp < 0 | P_exp > 1) |
    (M_conf < 0 | M_conf > 1) |
    (rho < -1 | rho > 1)
  ) {return (-Inf)}
  
  # define additional parameters
  P_belief <- 1 - P_exp
  A <- P_exp / sqrt(P_belief)
  B <- qnorm(M_conf) * sqrt( (1 - P_exp * P_belief) / P_belief )
  Sigma <- rbind(c(1, rho), c(rho, 1))
  
  # confidence rating cutoffs on probit scale (z scores)
  cutoffs <- 1:(k-1)
  z_cutoffs <- (qnorm(cutoffs/k) - B) / A
  
  # confidence x memory contingency table
  tab <- table(factor(data[,1], levels=1:k), 
               factor(data[,2], levels=c(-1, 0, 1)))
  
  # use for non-answer / no memory responses
  integrand1 <- function(inputs) {
    zm <- inputs[1]; ze <- inputs[2]
    pnorm(crit)^(n_alts-1) * dnorm(ze) * dnorm(zm, rho * ze, sqrt(1 - rho^2))
  }
  
  # use for incorrect responses
  integrand2 <- function(inputs) {
    zm <- inputs[1]; ze <- inputs[2]
    (1 - pnorm(max(crit, zm + mu_m))^(n_alts-1)) * dnorm(ze) * dnorm(zm, rho * ze, sqrt(1 - rho^2))
  }
  
  # use for correct responses
  integrand3 <- function(inputs) {
    zm <- inputs[1]; ze <- inputs[2]
    pnorm(zm + mu_m)^(n_alts-1) * dnorm(ze) * dnorm(zm, rho * ze, sqrt(1 - rho^2))
  }
  
  # compute log-likelihood for each JOL / memory combination
  LL <- 0
  for ( i in 1:k ) {
    n_1 <- tab[i, 1]; n0 <- tab[i, 2]; n1 <- tab[i, 3]
    if ( i == 1 ) {
      if ( n_1 > 0 ) {
        LL <- LL + n_1 * log(cubature::hcubature(f=integrand1,
                                                 lowerLimit=c(-Inf, -Inf),
                                                 upperLimit=c(crit - mu_m, z_cutoffs[i]))$integral)
      }
      if ( n0 > 0 ) {
        LL <- LL + n0 * log(cubature::hcubature(f=integrand2,
                                                lowerLimit=c(-Inf, -Inf),
                                                upperLimit=c(Inf, z_cutoffs[i]))$integral)
      }
      if ( n1 > 0 ) {
        LL <- LL + n1 * log(cubature::hcubature(f=integrand3,
                                                lowerLimit=c(crit - mu_m, -Inf),
                                                upperLimit=c(Inf, z_cutoffs[i]))$integral)
      }
    } else if ( (i > 1) & (i < k) ) {
      if ( n_1 > 0 ) {
        LL <- LL + n_1 * log(cubature::hcubature(f=integrand1,
                                                 lowerLimit=c(-Inf, z_cutoffs[i-1]),
                                                 upperLimit=c(crit - mu_m, z_cutoffs[i]))$integral)
      }
      if ( n0 > 0 ) {
        LL <- LL + n0 * log(cubature::hcubature(f=integrand2,
                                                lowerLimit=c(-Inf, z_cutoffs[i-1]),
                                                upperLimit=c(Inf, z_cutoffs[i]))$integral)
      }
      if ( n1 > 0 ) {
        LL <- LL + n1 * log(cubature::hcubature(f=integrand3,
                                                lowerLimit=c(crit - mu_m, z_cutoffs[i-1]),
                                                upperLimit=c(Inf, z_cutoffs[i]))$integral)
      }
    } else if ( i == k ) {
      if ( n_1 > 0 ) {
        LL <- LL + n_1 * log(cubature::hcubature(f=integrand1,
                                                 lowerLimit=c(-Inf, z_cutoffs[i-1]),
                                                 upperLimit=c(crit - mu_m, Inf))$integral)
      }
      if ( n0 > 0 ) {
        LL <- LL + n0 * log(cubature::hcubature(f=integrand2,
                                                lowerLimit=c(-Inf,z_cutoffs[i-1]),
                                                upperLimit=c(Inf, Inf))$integral)
      }
      if ( n1 > 0 ) {
        LL <- LL + n1 * log(cubature::hcubature(f=integrand3,
                                                lowerLimit=c(crit - mu_m, z_cutoffs[i-1]),
                                                upperLimit=c(Inf, Inf))$integral)
      }
    }
  }
  
  output <- ifelse(is.nan(LL), -Inf, LL)
  return ( output )
}


# Model log-prior 
BIM_MC_logprior <- function(params) {
  P_exp <- params[1]
  mu_m <- params[2]
  M_conf <- params[3]
  rho <- params[4]
  crit <- params[5]
  
  return (
    dunif(P_exp, log=TRUE) +
      dnorm(mu_m, log=TRUE) +
      dunif(M_conf, log=TRUE) +
      dunif(rho, min=-1, max=1, log=TRUE) +
      dnorm(crit, log=TRUE)
  )
}


# Model simulation function
#   N:   number of trials to simulate
#   rating: whether confidence ratings are continuous or discrete
#   memory: whether memory outcomes are based on 
#   n_alts: number of choice alternatives in the multiple-choice memory test
BIM_MC_simulate <- function(params, N, n_ratings, n_alts, return_table=FALSE) {
  P_exp <- params[1]
  mu_m <- params[2]
  M_conf <- params[3]
  rho <- params[4]
  crit <- params[5]
  
  output <- matrix(nrow=N, ncol=2)
  
  P_belief <- 1 - P_exp
  A <- P_exp / sqrt(P_belief)
  B <- qnorm(M_conf) * sqrt( (1 - P_exp * P_belief) / P_belief )
  
  for ( i in 1:N ) {
    
    # sample an objective memory strength and subjective processing experience
    # from their joint distribution
    X <- mvtnorm::rmvnorm(n=1, mean=c(0, 0), sigma=rbind(c(1, rho), c(rho, 1)))
    zm <- X[1]; ze <- X[2]
    m <- zm + mu_m
    
    # memory outcome (multiple-choice format)
    lures <- rnorm(n=n_alts-1)
    mem <- 1 * (m > crit) * (m > max(lures)) +    # 1 = correct choice
      0 * any((lures > m) * (lures > crit)) +     # 0 = incorrect choice
      -1 * (m < crit) * (max(lures) < crit)       # -1 = not remembered

    # confidence rating
    pred_conf <- pnorm(A * ze + B)
    criteria <- 1:(n_ratings - 1)/n_ratings
    conf <- sum(pred_conf >= criteria) + 1
    
    output[i,] <- c(conf, mem)
  }
  
  if (return_table) {
    return (table(factor(output[,1], levels=1:n_ratings),
                  factor(output[,2], levels=c(-1,0,1))))
  }
  return (output)
}


