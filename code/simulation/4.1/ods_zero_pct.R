
# Simulation 4.1: ODS benchmark
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

source(file.path(helper_dir, "data_generator_for_given_order.R"))
# set random seed
#set.seed(05012020)
# load necessary libraries
library(foreach)
library(doParallel)
library(doRNG)
library(igraph)
library(glmnet)

#----------------------------------------------------------------------------------------------------#
# simulate data
n_rep = 30

# set true parameters
n = 1000   # sample size
p = 50    # number of node 25, 50, 75
K = 1
data0 <- vector("list", n_rep)
dat  <- vector("list", n_rep)
# E_true <- vector("list", n_rep)
# graph <- vector("list", n_rep)

for(i in 1:n_rep){
  data0[[i]] <- generate_data_from_order(n, p, K, topo = NULL, seed = i)
  dat[[i]]  <- data0[[i]]$dat[[1]]
  # E_true[[i]] <- data[[i]]$E_true[[1]]
  # graph[[i]] <- data[[i]]$graph[[1]]
}
 
# set random seed
set.seed(05012020)

# # set true parameters
# E_true     = matrix(0, p, p)
# alpha_true = matrix(0, p, p)
# beta_true  = matrix(0, p, p)
# delta_true = rep(0, p)
# gamma_true = rep(0, p)
# 
# # generate a DAG with p edges
# n_edge = p
# while (sum(E_true == 1) < n_edge)
# {
#   id_edge = matrix(sample(1 : p, 2), ncol = 2)
#   E_true[id_edge] = 1
#   g_true = graph_from_adjacency_matrix(t(E_true))
#   if (!(is_dag(g_true)))
#     E_true[id_edge] = 0
# }
# 
# alpha_true[E_true == 1] = runif(n_edge, 0.3, 1)  # zero-inflation effects
# beta_true[E_true == 1]  = runif(n_edge, -1, -0.3) # Poisson log-rate effects
# delta_true[ ] = runif(p, -2.0, -1.0)  # zero-inflation intercepts
# gamma_true[ ] = runif(p,  1.0,  2.0)  # Poisson log-rate intercepts
# 
# 
# # generate n_rep datasets
# dat = list()
# for (rep in 1 : n_rep)
# {
#   data = matrix(0, n, p)
#   order_nodes = as_ids(topo_sort(g_true))
#   for (j in order_nodes)
#   {
#     pi = exp(data %*% alpha_true[j, ] + delta_true[j])
#     pi = pi / (1 + pi)
#     pi[is.nan(pi)] = 1
#     lambda  = exp(data %*% beta_true[j, ] + gamma_true[j])
#     data[ , j] = rpois(n, lambda) * (1 - rbinom(n, 1, pi))
#   }
#   
#   dat[[rep]] = data
# }
#----------------------------------------------------------------------------------------------------#
# implementation of [Park & Raskutti, 2015]
# tuning parameter
c0 = 0.004

cluster = makeCluster(10)
registerDoParallel(cluster)
ods = foreach(iter = 1 : n_rep, .packages = "glmnet") %dorng%
  {
    data = dat[[iter]]
    
    # step 1: estimate the undirected edges corresponding to the moralized graph with neighborhood sets
    moral_graph = matrix(0, p, p)
    for (j in 1 : p)
    {
      cv  = cv.glmnet(data[ , -j], data[ , j], family = "poisson")
      out = glmnet(data[ , -j], data[ , j], family = "poisson", lambda = cv$lambda.1se)
      moral_graph[j, -j] = as.vector(out$beta != 0) + 0
    }
    moral_graph = ((moral_graph != 0) & (t(moral_graph) != 0)) + 0
    
    # step 2: estimate causal ordering using overdispersion score
    dag = matrix(0, p, p)
    ordering = rep(NA, p)
    
    # 1st element of the causal ordering
    score = apply(data, 2, var) - colMeans(data)
    ordering[1] = which.min(score)
    
    # 2nd to (p-1)-th elements of the causal ordering
    for (j in 2 : (p - 1))
    {
      score = rep(0, p)
      score[ordering[1 : (j - 1)]] = Inf
      
      for (k in (1 : p)[-ordering[1 : (j - 1)]])
      {
        x_k = data[ , k]
        parents = intersect(which(moral_graph[k, ] != 0), ordering[1 : (j - 1)])
        
        if (length(parents) == 0)
        {
          score[k] = var(x_k) - mean(x_k)
        } else
        {
          id      = apply(cbind(rep("x", n), data[ , parents]) , 1 , paste, collapse = "-")
          uniq_id = unique(id)
          num = denom = 0
          
          for (l in 1 : length(uniq_id))
          {
            n_Xs = sum(id == uniq_id[l])
            if (n_Xs < c0 * n)
              next
            num   = num + n_Xs * (var(x_k[id == uniq_id[l]]) - mean(x_k[id == uniq_id[l]]))
            denom = denom + n_Xs 
          }
          
          score[k] = num / denom
        }
      }
      
      pi_j = which.min(score)
      ordering[j] = pi_j
      
      # step 3: estimate directed edges toward the j-th element of the causal ordering
      candidates = intersect(which(moral_graph[pi_j, ] != 0), ordering[1 : (j - 1)])
      
      if (length(candidates) <= 1)
      {
        dag[pi_j, candidates] <- 1
      } else 
      {
        cv  = cv.glmnet(data[ , candidates, drop = FALSE], data[ , pi_j], family = "poisson")
        out = glmnet(data[ , candidates], data[ , pi_j], family = "poisson", lambda = cv$lambda.1se)
        dag[pi_j, candidates] = as.vector(out$beta != 0) + 0
      }
    }
    
    # p-th element of the causal ordering
    pi_p = (1 : p)[-ordering[1 : (p - 1)]]
    ordering[p] = pi_p
    
    # directed edges toward the p-th element of the causal ordering
    candidates = intersect(which(moral_graph[pi_p, ] != 0), ordering[1 : (p - 1)])
    if (length(candidates) <= 1)
    {
      dag[pi_p, candidates] = 1
    } else 
    {
      cv  = cv.glmnet(data[ , candidates], data[ , pi_p], family = "poisson")
      out = glmnet(data[ , candidates], data[ , pi_p], family = "poisson", lambda = cv$lambda.1se)
      dag[pi_p, candidates] = as.vector(out$beta != 0) + 0
    }
    
    # print progress
    cat("iter=", iter, "\n")
    
    # return the estimated dag
    dag
  }

stopCluster(cluster)

save(ods, data0, file = file.path(output_dir, sprintf("ods_z%d.RData", p)))
