##### PARAMETER RECOVERY SCRIPT  #####
##### Author: William Hayes #####

source('model_utils.R')
source('DEmcmc_single.R')

library(parallel)
library(posterior)

##### Model functions #####

model = list(
  LL = BIM_MC_loglik,
  Prior = BIM_MC_logprior
)


##### Full parameter recovery test #####

set.seed(735911)

# sample data-generating parameter sets
P_exp <- runif(n=100, .05, .95)
mu_m <- runif(n=100, -3, 3)
M_conf <- runif(n=100, .05, .95)
rho <- runif(n=100, -.95, .95)
crit <- runif(n=100, -3, 3)
params <- cbind(P_exp=P_exp, mu_m=mu_m, M_conf=M_conf, rho=rho, crit=crit)


recover_BIM_MC <- function(params, model, simulator, N, n_ratings, n_alts, progress_bar=FALSE) {
  # simulate data using true parameter values
  true_params <- unname(params) 
  sims <- simulator(params=true_params, N=N, n_ratings=n_ratings, n_alts=n_alts)
  
  # set up DE-MCMC parameters
  param_names <- c('P_exp', 'mu_m', 'M_conf', 'rho', 'crit')
  m <- length(param_names) * 2
  k <- 3
  n_chains <- m * k
  subpop.ids <- rep(1:k, each=m)
  inits <- cbind(P_exp = runif(n_chains),
                 mu_m = rnorm(n_chains),
                 M_conf = runif(n_chains),
                 rho = runif(n_chains, min=-1, max=1),
                 crit = rnorm(n_chains))
  
  # sample from posterior
  results <- DEmcmc(model = model,
                    data = sims,
                    k = n_ratings,
                    n_alts = n_alts,
                    init.chains = inits,
                    par.names = param_names,
                    n.burnin = 1000,
                    n.samples = 2000,
                    subpop.ids = subpop.ids,
                    p.migration = 0.05,
                    progress = progress_bar)
  
  # compute convergence statistic
  rhats <- apply(results$samples, MARGIN=3, FUN=rhat_basic)
  
  # compute posterior means
  posteriorMeans <- apply(results$samples, MARGIN=3, FUN=mean)
  
  # 95% credible intervals
  credibleIntervals <- apply(results$samples, MARGIN=3, FUN=quantile, prob=c(.025, .975))
  
  # simulated data
  data_table <- table(JOL=factor(sims[,1], levels=1:n_ratings),
                      MEM=factor(sims[,2], levels=c(-1,0,1)))
  
  return(list(true=params,
              recovered=posteriorMeans,
              CIs=credibleIntervals,
              Rhats=rhats,
              data=data_table))
}


# set up cluster
n_cores <- min(14, detectCores())
cat('No. of cores: ', n_cores)
cl <- makeCluster(n_cores)
clusterExport(cl, varlist = c('rhat_basic', 'DEmcmc', 'updateParams', 'crossover', 'migration'))

system.time(
  paramRecoveryResults <- parApply(
    cl = cl,
    X = params,
    MARGIN = 1,
    FUN = recover_BIM_MC,
    model = model,
    simulator = BIM_MC_simulate,
    N = 60,
    n_ratings = 3,
    n_alts = 4
  )
)

saveRDS(paramRecoveryResults, file='../../Results/BIM_MC_recovery.RDS')

stopCluster(cl)
