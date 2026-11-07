# Simulation 4.1: MRS benchmark
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

# load necessary libraries
library(igraph)
library(MXM)

#----------------------------------------------------------------------------------------------------#
# simulate data
n_rep = 30

# set true parameters
n = 1000   # sample size
p = 50     # number of node 25, 50, 75
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

#----------------------------------------------------------------------------------------------------#
# implementation of [Park & Park, 2019]
# assume hyper-Poisson 
# tuning parameter
Nmin = 2

# define the CMR function
CMRfunction = function(x, a, b)
{
  x * x * ((a + 1) / a) * (b / (b + 1))
}

# MRS algorithm with r = 2
mrs = list()
for (iter in 1 : n_rep)
{
  data = dat[[iter]]
  
  # step 1: estimate the skeleton of the graph 
  skeleton = mmhc.skel(data, test = "testIndPois", nc = 10)$G
  
  # step 2: estimate an ordering of the graph using the second moments ratio scores
  dag      = matrix(0, p, p)
  ordering = rep(NA, p)
  
  # 1st element of the ordering
  score = rep(0, p) 
  b     = apply(data, 2, var) / colMeans(data)
  for (j in 1 : p)
  {
    moment1  = mean(data[ , j])
    moment2  = mean(data[ , j] * data[ , j])
    score[j] = moment2 / (CMRfunction(moment1, 1, b[j]) + moment1)
  }
  ordering[1] = which.min(score)
  
  # 2nd to (p-1)-th elements of the ordering
  for (m in 2 : (p - 1))
  {
    score = rep(0, p)
    score[ordering[1 : (m - 1)]] = Inf
    
    for (j in (1 : p)[-ordering[1 : (m - 1)]])
    {
      x_j     = data[ , j]
      parents = intersect(which(skeleton[j, ] != 0), ordering[1 : (m - 1)])
      
      if (length(parents) == 0)
      {
        moment1  = mean(x_j)
        moment2  = mean(x_j * x_j)
        score[j] = moment2 / (CMRfunction(moment1, 1, b[j]) + moment1)
      } else
      {
        id      = apply(cbind(rep("x", n), data[ , parents]) , 1 , paste, collapse = "-")
        uniq_id = unique(id)
        num = denom = 0
        
        for (s in 1 : length(uniq_id))
        {
          n_Xs = sum(id == uniq_id[s])
          if (n_Xs < Nmin)
            next
          
          denom   = denom + n_Xs 
          moment1 = mean(x_j[id == uniq_id[s]])
          moment2 = mean(x_j[id == uniq_id[s]] * x_j[id == uniq_id[s]])
          if (moment1 == 0)
            num = num + n_Xs
          else
            num = num + n_Xs * moment2 / (CMRfunction(moment1, 1, b[j]) + moment1)
        }
        
        score[j] = num / denom
      }
    }
    
    ordering[m] = which.min(score)
  }
  
  # last element of the ordering
  ordering[p] = (1 : p)[-ordering[1 : (p - 1)]]
  
  # estimate the edge sets
  for (m in 2 : p)
  {
    edges = intersect(which(skeleton[ordering[m], ] != 0), ordering[1 : (m - 1)])
    dag[ordering[m], edges] = 1
  }
  
  # save the estimated dag
  mrs[[iter]] = dag
  
  # print progress
  cat("iter=", iter, "\n")
}

save(mrs, data0, file = file.path(output_dir, sprintf("mrs_z%d.RData", p)))
