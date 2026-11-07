library(igraph)

# =====================================================================
# Calibrate the zero-inflation intercept
#
# For a zero-inflated negative-binomial variable:
#
# P(Y = 0) = pi + (1 - pi) * P_NB(Y = 0)
#
# This function chooses the intercept delta so that the average
# expected zero proportion is approximately target_zero.
# =====================================================================
calibrate_zinb_intercept <- function(zero_linear_predictor, mu, theta, 
                                     target_zero = 0.50, lower = -25, upper = 25) {
  
  # Probability of zero from the negative-binomial component
  nb_zero_probability <- dnbinom(x = 0, size = theta, mu = mu)
  
  objective <- function(delta) {
    
    pi_zero <- plogis(zero_linear_predictor + delta)
    
    expected_zero_probability <- pi_zero + (1 - pi_zero) * nb_zero_probability
    
    mean(expected_zero_probability) - target_zero
  }
  
  lower_value <- objective(lower)
  upper_value <- objective(upper)
  
  # The NB component may already produce at least target_zero zeros.
  # In that case, use almost no additional zero inflation.
  if (lower_value >= 0) {return(lower)}
  
  if (upper_value <= 0) {return(upper)}
  
  uniroot(f = objective, interval = c(lower, upper), tol = 1e-8)$root
}


# =====================================================================
# Generate data from a zero-inflated negative-binomial Bayesian network
# =====================================================================
generate_data_from_order <- function(n, p = NULL, K, topo = NULL, seed = NULL,
                                     target_zero = 0.50, dispersion_range = c(1, 5), 
                                     mu_floor = 1.05) {
  
  # -------------------------------------------------------------------
  # Validate inputs
  # -------------------------------------------------------------------
  if (is.null(topo)) {
    if (is.null(p)) {stop("Provide either 'p' or a 'topo' vector.")}
    topo <- seq_len(p)
  }
  
  p <- length(topo)
  
  if (!setequal(topo, seq_len(p))) {stop("topo must be a permutation of 1, ..., p.")}
  
  if (length(dispersion_range) != 2L) {stop("dispersion_range must contain two values.")}
  
  if (dispersion_range[1] <= 0 ||
      dispersion_range[2] <= dispersion_range[1]) {
    stop("Use a valid positive dispersion interval.")
  }
  
  if (target_zero <= 0 || target_zero >= 1) {stop("target_zero must be between 0 and 1.")}
  
  # Exactly p edges, as in the original generator
  number_of_edges <- p
  
  # -------------------------------------------------------------------
  # Output storage
  # -------------------------------------------------------------------
  dat <- vector("list", K)
  E_true <- vector("list", K)
  g_true <- vector("list", K)
  
  alpha_list <- vector("list", K)
  beta_list <- vector("list", K)
  delta_list <- vector("list", K)
  gamma_list <- vector("list", K)
  theta_list <- vector("list", K)
  
  observed_zero_rate <- numeric(K)
  
  # -------------------------------------------------------------------
  # Generate K datasets
  # -------------------------------------------------------------------
  for (k in seq_len(K)) {
    
    if (!is.null(seed)) {set.seed(seed + k - 1L)}
    
    # ---------------------------------------------------------------
    # 1. Generate a DAG consistent with the supplied order
    # ---------------------------------------------------------------
    E <- make_E_from_order(topo = topo, n_edge = number_of_edges)
    
    stopifnot(is_E_consistent_with_order(E, topo))
    
    # E[child, parent] = 1.
    # Transpose for igraph so that edges are parent -> child.
    g <- graph_from_adjacency_matrix(t(E), mode = "directed", diag = FALSE)
    
    stopifnot(is_dag(g))
    
    # ---------------------------------------------------------------
    # 2. Initialize the count-data matrix
    # ---------------------------------------------------------------
    X <- matrix(0L, nrow = n, ncol = p)
    colnames(X) <- paste0("X", seq_len(p))
    
    # ---------------------------------------------------------------
    # 3. Generate regression parameters
    # ---------------------------------------------------------------
    actual_number_of_edges <- sum(E == 1L)
    
    # Zero-inflation coefficients
    alpha_true <- matrix(
      0,
      nrow = p,
      ncol = p
    )
    
    # Negative-binomial mean coefficients
    beta_true <- matrix(
      0,
      nrow = p,
      ncol = p
    )
    
    if (actual_number_of_edges > 0L) {
      
      alpha_true[E == 1L] <- runif(
        actual_number_of_edges,
        min = 0.3,
        max = 1.0
      )
      
      beta_true[E == 1L] <- runif(
        actual_number_of_edges,
        min = -1.0,
        max = -0.3
      )
    }
    
    # Mean-model intercepts
    gamma_true <- runif(
      p,
      min = 1.0,
      max = 2.0
    )
    
    # Dispersion parameters:
    # theta_j ~ Uniform(1, 5)
    theta_true <- runif(
      p,
      min = dispersion_range[1],
      max = dispersion_range[2]
    )
    
    # These intercepts are calibrated below to obtain about 50% zeros.
    delta_true <- numeric(p)
    
    # ---------------------------------------------------------------
    # 4. Generate variables in topological order
    # ---------------------------------------------------------------
    for (j in topo) {
      
      # Zero-inflation linear predictor without its intercept
      zero_linear_predictor <- as.vector(X %*% alpha_true[j, ])
      
      # Negative-binomial log-mean
      eta_mu <- as.vector(X %*% beta_true[j, ] + gamma_true[j])
      mu <- exp(eta_mu)
      
      # Prevent the NB component itself from creating much more than
      # 50% zeros before zero inflation is applied.
      mu <- pmax(mu, mu_floor)
      
      # Select the zero-inflation intercept that gives an expected
      # zero proportion near target_zero.
      delta_true[j] <- calibrate_zinb_intercept(
        zero_linear_predictor = zero_linear_predictor,
        mu = mu,
        theta = theta_true[j],
        target_zero = target_zero
      )
      
      # Structural-zero probability
      pi_zero <- plogis(zero_linear_predictor + delta_true[j])
      
      # Draw the negative-binomial count component
      nb_count <- rnbinom(n = n, size = theta_true[j], mu = mu)
      
      # Draw structural-zero indicators
      structural_zero <- rbinom(n = n, size = 1L, prob = pi_zero)
      
      # Zero-inflated negative-binomial observation
      X[, j] <- ifelse(structural_zero == 1L, 0L, nb_count)
    }
    
    # ---------------------------------------------------------------
    # 5. Store dataset and true parameters
    # ---------------------------------------------------------------
    dat[[k]] <- X
    E_true[[k]] <- E
    g_true[[k]] <- g
    
    alpha_list[[k]] <- alpha_true
    beta_list[[k]] <- beta_true
    delta_list[[k]] <- delta_true
    gamma_list[[k]] <- gamma_true
    theta_list[[k]] <- theta_true
    
    observed_zero_rate[k] <- mean(X == 0)
    
    # message(
    #   sprintf(
    #     "Dataset %d: observed zero proportion = %.3f; theta range = [%.3f, %.3f]",
    #     k,
    #     observed_zero_rate[k],
    #     min(theta_true),
    #     max(theta_true)
    #   )
    # )
  }
  
  # -------------------------------------------------------------------
  # Return generated data and true parameter values
  # -------------------------------------------------------------------
list(E_true = E_true, order = topo, dat = dat, graph = g_true, alpha_true = alpha_list,
       beta_true = beta_list, delta_true = delta_list, gamma_true = gamma_list,
       theta_true = theta_list, zero_rate = observed_zero_rate, target_zero = target_zero,
       dispersion_range = dispersion_range)
}


# =====================================================================
# Generate an adjacency matrix compatible with a topological order
#
# Convention:
# E[child, parent] = 1
# =====================================================================
make_E_from_order <- function(topo, n_edge) {
  
  p <- length(topo)
  maximum_edges <- p * (p - 1L) / 2L
  
  if (n_edge > maximum_edges) {stop("n_edge cannot exceed p * (p - 1) / 2.")}
  
  E <- matrix(0L, nrow = p, ncol = p)
  
  while (sum(E) < n_edge) {
    
    # The first sampled position is earlier in the topological order
    position <- sort(sample.int(p, size = 2L, replace = FALSE))
    
    parent <- topo[position[1L]]
    child <- topo[position[2L]]
    
    if (E[child, parent] == 0L) {E[child, parent] <- 1L}
  }
  
  E
}


# =====================================================================
# Check adjacency compatibility with a topological order
# =====================================================================
is_E_consistent_with_order <- function(E, topo) {
  
  p <- nrow(E)
  
  stopifnot(ncol(E) == p, length(topo) == p)
  position <- integer(p)
  position[topo] <- seq_len(p)
  
  edge_indices <- which(E == 1L, arr.ind = TRUE)
  
  if (nrow(edge_indices) == 0L) {return(TRUE)}
  
  # Columns contain parents; rows contain children
  all(position[edge_indices[, 2L]] < position[edge_indices[, 1L]])
}


# data <- generate_data_from_order(n = 100, p = 5, K = 5, topo = NULL, seed = NULL,
#                                      target_zero = 0.50, dispersion_range = c(1, 5), 
#                                      mu_floor = 1.05)
# is_E_consistent_with_order(data$E_true[[2]], data$order)
