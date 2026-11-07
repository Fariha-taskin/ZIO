# Simulation 4.1: ZIPBN benchmark
# Script location:ZIO/code/simulation/4.1
project_dir <- "ZIO"
if (!dir.exists(file.path(project_dir, "code", "helper"))) {
  if (dir.exists(file.path("code", "helper"))) {
    project_dir <- "."
  } else {
    stop("Cannot locate ZIO/code/helper. Run from ZIO/ or its parent directory.")
  }
}
helper_dir <- file.path(project_dir, "code", "helper")
output_dir <- file.path(project_dir, "output", "simulation", "4.1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# load necessary libraries
library(foreach)
library(doParallel)
library(doRNG)
library(igraph)
library(pscl)

# load necessary functions
source(file.path(helper_dir, "data_generator_for_given_order.R"))
source(file.path(helper_dir, "ZIPBNfunctions.R"))

#--------------------------------------------------------------------------------------------
# simulate data
n_rep = 30

# set true parameters
n = 1000   # sample size
p = 50     # number of node 25, 50, 75
K = 1
data <- vector("list", n_rep)
dat  <- vector("list", n_rep)

for(i in 1:n_rep){
  data[[i]] <- generate_data_from_order(n, p, K, topo = NULL, seed = i)
  dat[[i]]  <- data[[i]]$dat[[1]]
}

# set random seed
set.seed(05012020)


#----------------------------------------------------------------------------------------
# run the parallel-tempered MCMC 
cluster = makeCluster(30)
registerDoParallel(cluster)

zipbn = foreach(rep = 1 : n_rep, .packages = c("igraph", "pscl")) %dorng%
  {
    data = dat[[rep]]
    
    # necessary for the MCMC implementation
    MCMC      = 3000
    BURNIN    = 1500
    n_chains  = 10
    prob_swap = 0.1
    n = nrow(data); p = ncol(data)
    
    # sequence of temperatures for the parallel tempering
    temps = exp(seq(0, 1, length.out = n_chains))
    
    # hyperparameters for priors
    priors <- list(
      tau = c(0.01, 0.01),
      rho = c(0.5, 0.5)
    )
    
    # Initial proposal tuning parameters
    tuning <- lapply(seq_len(n_chains), function(m) {
      
      scale_m <- 2^((m - 1) / (n_chains - 1))
      
      list(
        sigma_alpha = 0.65 * scale_m,
        sigma_beta = 0.325 * scale_m,
        sigma_delta = 0.475 * scale_m,
        sigma_gamma = 0.125 * scale_m,
        sigma_E = c(0.05, 0.05, 0.475, 0.125),
        sigma_E_rev = c(0.5, 0.5, 0.95, 0.25)
      )
    })
    
    # starting values for MCMC
    starting = list()
    starting$alpha = matrix(0, p, p)
    starting$beta  = matrix(0, p, p)
    starting$delta = rep(0, p)
    starting$gamma = rep(0, p)
    starting$E     = matrix(0, p, p)
    starting$tau   = c(10, 10, 1, 1)
    starting$rho   = 1 / p
    for (j in 1 : p)
    {
      mle = zeroinfl(data[ , j] ~ 1 | 1, dist = "poisson")
      starting$delta[j] = mle$coefficients$zero
      starting$gamma[j] = mle$coefficients$count
    }
    
    # calculate logitPi and logLambda for each group with starting values
    starting$logitPi   = tcrossprod(data, starting$alpha) + 
      matrix(starting$delta, n, p, byrow = TRUE)
    starting$logLambda = tcrossprod(data, starting$beta) + 
      matrix(starting$gamma, n, p, byrow = TRUE)
    
    # initialize parameters
    param = list()
    for (m in 1 : n_chains)
    {
      param[[m]] = starting
    }
    
    # initialize MCMC samples
    mcmc_samples = list()
    mcmc_samples$alpha = array(NA, dim = c(p, p, MCMC))
    mcmc_samples$beta  = array(NA, dim = c(p, p, MCMC))
    mcmc_samples$delta = matrix(NA, p, MCMC)
    mcmc_samples$gamma = matrix(NA, p, MCMC)
    mcmc_samples$E     = array(NA, dim = c(p, p, MCMC))
    mcmc_samples$tau   = matrix(NA, 4, MCMC)
    mcmc_samples$rho   = rep(NA, MCMC)
    
    # run MCMC iterations
    for (t in 1 : MCMC)
    {
      if (runif(1) > prob_swap)
      {
        # perform one-step update for all chains 
        for (m in 1 : n_chains)
        {
          # update each chain
          out = update_chain(data, param[[m]], tuning[[m]], priors, temps[m])
          param[[m]] = out$param
        }
      } else
      {
        # propose a swapping move
        # randomly choose chains to swap
        rand = sort(sample(1 : n_chains, 2))
        m1   = rand[1]
        m2   = rand[2]
        
        # calculate MH ratio for swapping
        ratio_swap = log_dZIPBN(data, param[[m1]], priors, temps[m2]) + 
          log_dZIPBN(data, param[[m2]], priors, temps[m1]) -
          log_dZIPBN(data, param[[m1]], priors, temps[m1]) - 
          log_dZIPBN(data, param[[m2]], priors, temps[m2])
        ratio_swap = exp(ratio_swap)
        
        # accept swapping
        if (is.nan(ratio_swap)) ratio_swap = 0
        if (runif(1) < min(1, ratio_swap))
        {
          cat("swap the states of chains", m1, "and", m2, "\n")
          param_m1    = param[[m1]]
          param[[m1]] = param[[m2]]
          param[[m2]] = param_m1
        }   
      }
      
      # save MCMC samples for the cold chain
      mcmc_samples$alpha[ , , t] = param[[1]]$alpha
      mcmc_samples$beta[ , , t]  = param[[1]]$beta
      mcmc_samples$delta[ , t]   = param[[1]]$delta
      mcmc_samples$gamma[ , t]   = param[[1]]$gamma
      mcmc_samples$E[ , , t]     = param[[1]]$E
      mcmc_samples$tau[ , t]     = param[[1]]$tau
      mcmc_samples$rho[t]        = param[[1]]$rho
      
      # print progress
      if (t %% 100 == 0)
        cat("iter =", t, "\n")
    }
    
    # return MCMC samples for the cold chain
    mcmc_samples
  }

# Stop workers and save the benchmark result in the shared output folder.
stopCluster(cluster)
save(zipbn, data, file = file.path(output_dir, sprintf("zipbn_z%d.RData", p)))

