library(igraph)

generate_data_from_order <- function(
    n,
    p          = NULL,
    K,
    topo       = NULL,      # e.g., c(1,3,4,2,6,5,8,9,7,10)
    #edge_prob  = 0.3,       # used if 'indeg' is NULL
    #indeg      = NULL,      # set to force exactly k parents for i>1
    seed       = NULL       # optional base seed for reproducibility
) {
  # Resolve topo / p
  if (is.null(topo)) {
    if (is.null(p)) stop("Provide either 'p' or a 'topo' vector.")
    topo <- 1:p
  }
  p <- length(topo)
  n_edge <- p
  # outputs
  dat    <- vector("list", K)
  E_true <- vector("list", K)
  g_true <- vector("list", K)
  
  for (k in seq_len(K)) {
    if (!is.null(seed)) set.seed(seed + k - 1L)  # optional: reproducible per dataset
    
    # 1) adjacency consistent with topo
    E <- make_E_from_order(topo, n_edge)
    #E <- make_E_from_order(topo, edge_prob = edge_prob, indeg = indeg)
    stopifnot(is_E_consistent_with_order(E, topo))
    
    # 2) igraph DAG (transpose so edges are parent -> child)
    g <- graph_from_adjacency_matrix(t(E), mode = "directed", diag = FALSE)
    stopifnot(is_dag(g))
    
    # 3) simulate ZIP-Poisson BN data
    X <- matrix(0, nrow = n, ncol = p)
    
    n_edge <- sum(E == 1L)
    alpha_true <- matrix(0, p, p)
    beta_true  <- matrix(0, p, p)
    
    #----------------------- ~25% zeros in the data ----------------------------
    # if (n_edge > 0) {
    #   alpha_true[E == 1L] <- runif(n_edge, -1.0, -0.3)  # zero-inflation effects
    #   beta_true [E == 1L] <- runif(n_edge, -0.4, -0.2)  # Poisson log-rate effects
    # }
    # delta_true <- runif(p, -2.0, -1.5)  # zero-inflation intercepts
    # gamma_true <- runif(p,  1.5,  2.0)  # Poisson log-rate intercepts
    #----------------------- ~50% zeros in the data ----------------------------
    if (n_edge > 0) {
      alpha_true[E == 1L] <- runif(n_edge, 0.3, 1)  # zero-inflation effects
      beta_true [E == 1L] <- runif(n_edge, -1, -0.3)  # Poisson log-rate effects
    }
    delta_true <- runif(p, -2.0, -1.0)  # zero-inflation intercepts
    gamma_true <- runif(p,  1.0,  2.0)  # Poisson log-rate intercepts
    
    #----------------------- ~75% zeros in the data ----------------------------
    # if (n_edge > 0) {
    #   alpha_true[E == 1L] <- runif(n_edge, 0.3, 1)  # zero-inflation effects
    #   beta_true [E == 1L] <- runif(n_edge, -1, -0.3)  # Poisson log-rate effects
    # }
    # delta_true <- runif(p, 0.5, 1)  # zero-inflation intercepts
    # gamma_true <- runif(p, 1, 1.5)  # Poisson log-rate intercepts
    
    # generate in topo order so parents are filled first
    for (j in topo) {
      eta_p <- as.vector(X %*% alpha_true[j, ] + delta_true[j])  # logit(pi)
      eta_l <- as.vector(X %*% beta_true [j, ] + gamma_true[j])  # log(lambda)
      p_bin  <- plogis(eta_p)
      lambda <- exp(eta_l)
      X[, j] <- rpois(n, lambda) * (1 - rbinom(n, 1, p_bin))    # ZIP: zero with prob pi
    }
    
    dat[[k]]    <- X
    E_true[[k]] <- E
    g_true[[k]] <- g
  }
  
  list(E_true = E_true, order = topo, dat = dat, graph = g_true)
}

# ---- helpers ----

make_E_from_order <- function(topo, n_edge) {
  p <- length(topo)
  E <- matrix(0L, p, p)
  
  while (sum(E == 1L) < n_edge) {
    # sample a parent (earlier in topo) and a child (later in topo)
    pos <- sort(sample(p, 2))         # ensures topo[pos[1]] is earlier
    parent <- topo[pos[1]]
    child  <- topo[pos[2]]
    
    # avoid duplicates
    if (E[child, parent] == 0L) {
      E[child, parent] <- 1L
    }
  }
  E
}



# make_E_from_order <- function(topo, edge_prob = 0.3, indeg = NULL) {
#   p <- length(topo)
#   E <- matrix(0L, p, p)
#   for (i in 2:p) {
#     child   <- topo[i]
#     parents <- topo[seq_len(i - 1)]
#     if (length(parents) == 0) next
#     if (is.null(indeg)) {
#       E[child, parents] <- rbinom(length(parents), 1, edge_prob)
#     } else {
#       k <- min(indeg, length(parents))
#       if (k > 0) {
#         chosen <- sample(parents, k, replace = FALSE)
#         E[child, chosen] <- 1L
#       }
#     }
#   }
#   E
# }

is_E_consistent_with_order <- function(E, topo) {
  p <- nrow(E); stopifnot(ncol(E) == p, length(topo) == p)
  pos <- integer(p); pos[topo] <- seq_len(p)
  idx <- which(E == 1L, arr.ind = TRUE)   # rows=j (child), cols=k (parent)
  if (nrow(idx) == 0) return(TRUE)
  all(pos[idx[, 2]] < pos[idx[, 1]])
}
