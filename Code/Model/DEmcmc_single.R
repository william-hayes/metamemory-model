##### FUNCTIONS FOR DIFFERENTIAL EVOLUTION MCMC #####
##### Author: William Hayes                     #####


#Function for performing crossover step
crossover <- function(X, subpop.ids, ll.fn, lp.fn, data, current.ll, gam, sig, ...) {
  n.chains <- nrow(X)
  chain.ids <- 1:n.chains
  accept <- rep(0, n.chains)
  for(i in 1:n.chains) {
    theta <- X[i,]
    mate.ids <- setdiff(chain.ids[subpop.ids == subpop.ids[i]], chain.ids[i])
    rands <- sample(mate.ids, 2)
    proposal <- theta + gam * (X[rands[1],] - X[rands[2],]) + rnorm(length(theta), 0, sig)
    proposal.ll <- ll.fn(proposal, data, ...)
    target.theta <- current.ll[i] + lp.fn(theta)
    target.proposal <- proposal.ll + lp.fn(proposal)
    if(is.nan(target.proposal - target.theta)) {
      next
    }
    if(runif(1) <= exp(target.proposal - target.theta)) {
      X[i,] <- proposal
      accept[i] <- 1
      current.ll[i] <- proposal.ll
    }
  }
  output <- list(X = X, accept = accept, loglik = current.ll)
  return(output)
}


#Function for performing migration step
migration <- function(X, subpop.ids, current.ll) {
  n.chains <- nrow(X)
  chain.ids <- 1:n.chains
  accept <- rep(0, n.chains)
  n.subpops <- length(unique(subpop.ids))
  # permute sub-population indices
  subpops <- sample(unique(subpop.ids))
  emigrant.ids <- integer(n.subpops)
  # randomly sample a chain from each sub-population to migrate
  for( i in 1:n.subpops ) {
    emigrant.ids[i] <- sample(chain.ids[subpop.ids == subpops[i]], 1)
  }
  # migration ring (left-to-right)
  ring <- c(emigrant.ids, emigrant.ids[1])
  # update population 
  X[ring[2:(n.subpops+1)],] <- X[ring[1:n.subpops],]
  # update current log-likelihood values
  current.ll[ring[2:(n.subpops+1)]] <- current.ll[ring[1:n.subpops]]
  accept[emigrant.ids] <- 1
  output <- list(X = X, accept = accept, loglik = current.ll)
  return(output)
}


#function for updating parameters
updateParams <- function(chains, subpop.ids, ll.fn, lp.fn, data, current.ll, gam, sig, p.migration, ...) {
  
  if(gam == "auto") {
    gam <- 2.38 / sqrt(2 * ncol(chains))
  }else if(gam == "random") {
    gam <- runif(1, .5, 1)
  }
  
  # crossover or migration?
  if(runif(1) <= p.migration) {
    output <- migration(chains, subpop.ids, current.ll)
  }else {
    output <- crossover(chains, subpop.ids, ll.fn, lp.fn, data, current.ll, gam, sig, ...)
  }
  
  return(output)
  
}

#Function for DE-MCMC sampling for an individual fit
# end.migration=floor(n.burnin/2)
DEmcmc <- function(model, data, ..., init.chains, par.names, n.burnin, n.samples, 
                   gam="auto", sig=.001, subpop.ids=rep(1,nrow(init.chains)), 
                   p.migration=.10, save.burnin=FALSE, progress=TRUE) {
  
  X <- init.chains; colnames(X) <- par.names
  n.chains <- nrow(X)
  n.params <- ncol(X)
  n.subpops <- length(unique(subpop.ids))
  n.iter <- n.burnin + n.samples
  
  samples.out <- array(dim=c(n.samples, n.chains, n.params),
                       dimnames=list(iter=1:n.samples, chain=1:n.chains, param=par.names))
  loglik.out <- array(dim=c(n.samples, n.chains),
                      dimnames=list(iter=1:n.samples, chain=1:n.chains))
  accept.out <- array(dim=c(n.samples, n.chains),
                      dimnames=list(iter=1:n.samples, chain=1:n.chains))
  
  if(save.burnin) { 
    burnin.out <- array(dim=c(n.burnin, n.chains, n.params), 
                        dimnames=list(iter=1:n.burnin, chain=1:n.chains, param=par.names))
  }else {
    burnin.out <- NULL
  }
  
  #Compute log-likelihoods given the initial state of each chain
  LL <- numeric(n.chains)
  for(i in 1:n.chains) {
    LL[i] <- model$LL(params=X[i,], data=data, ...)
  }
  
  #Initialize progress bar
  if(progress) {
    pbar <- txtProgressBar(min=1, max=n.iter, style=3)
    chkpts <- seq(0, 1, .01) * n.iter
    cat(paste("#chains:", n.chains, "#subpops:", n.subpops, "\n"))
    cat(paste("Sampling underway...", "\n"))
  }
  
  for(i in 1:n.iter) {
    
    temp <- updateParams(X, subpop.ids, model$LL, model$Prior, data, LL, gam, sig, p.migration, ...)
                     
    X <- temp$X
    LL <- temp$loglik
    if(i > n.burnin) { accept.out[i-n.burnin,] <- temp$accept }
  
    #Save samples
    if((i <= n.burnin) & save.burnin) {burnin.out[i,,] <- X}
    if(i > n.burnin) {
      samples.out[i-n.burnin,,] <- X
      loglik.out[i-n.burnin,] <- LL
    }
    
    #Update progress bar 
    if(progress) {
      if(i >= chkpts[1]) {
        setTxtProgressBar(pbar, i)
        chkpts <- chkpts[-1]
      }
    }
    
  }
  
  if(progress) {
    cat(paste("\n", "Sampling complete.", "\n"))
    close(pbar)
  }
  
  output <- list(burnin = burnin.out, samples = samples.out, loglik = loglik.out, 
                 accept_rates = apply(accept.out, 2, mean))
  
  return(output)   
  
}