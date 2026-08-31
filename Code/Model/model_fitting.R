##### MODEL FITTING SCRIPT  #####
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

##### Load data #####

dat <- read.csv('../../Data/modeling_data_v2.csv', header=T)


##### Individual-level fits #####

set.seed(92474)

fit_BIM_MC <- function(index, model, data, n_ratings, n_alts, progress_bar=FALSE) {
  # individual data
  d <- with(subset(data, idnum == index), cbind(JOL, MEM))
  
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
                    data = d,
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
  
  # calculate DIC
  LL_pointEstimate <- model$LL(posteriorMeans, d, n_ratings, n_alts)
  pDIC <- 2 * (LL_pointEstimate - mean(results$loglik))
  dic <- -2 * LL_pointEstimate + 2 * pDIC
  
  return(list(idnum=index,
              point_estimates=posteriorMeans,
              CIs=credibleIntervals,
              Rhats=rhats,
              DIC=dic))
}


#fit <- fit_BIM_MC(1, model, dat, n_ratings=3, n_alts=4, progress_bar=TRUE)


# set up cluster
n_cores <- min(14, detectCores())
cat('No. of cores: ', n_cores)
cl <- makeCluster(n_cores)
clusterExport(cl, varlist = c('rhat_basic', 'DEmcmc', 'updateParams', 'crossover', 'migration'))


IDs <- unique(dat$idnum)

system.time(
  fitResults <- parLapply(
    cl = cl,
    X = IDs,
    fun = fit_BIM_MC,
    model = model,
    data = dat,
    n_ratings = 3,
    n_alts = 4
  )
)

saveRDS(fitResults, file='../../Results/BIM_MC_fits.RDS')

stopCluster(cl)



