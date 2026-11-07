############################################################
## 1. Helper functions
############################################################

E_to_adj_parent_child <- function(E) {
  t(E)
}

is_valid_order_E <- function(ordering, E) {
  is_E_consistent_with_order(E, ordering)
}

all_topological_orders <- function(E) {
  p <- nrow(E)
  children_of <- vector("list", p)
  indeg <- integer(p)
  
  for (child in 1:p) {
    pa <- which(E[child, ] == 1L)
    indeg[child] <- length(pa)
    for (parent in pa) {
      children_of[[parent]] <- c(children_of[[parent]], child)
    }
  }
  
  out <- list()
  
  backtrack <- function(current_order, indeg_now, used) {
    if (length(current_order) == p) {
      out[[length(out) + 1L]] <<- current_order
      return(NULL)
    }
    
    available <- which(indeg_now == 0 & !used)
    
    for (v in available) {
      indeg_next <- indeg_now
      used_next <- used
      used_next[v] <- TRUE
      
      if (length(children_of[[v]]) > 0) {
        indeg_next[children_of[[v]]] <- indeg_next[children_of[[v]]] - 1L
      }
      
      backtrack(
        current_order = c(current_order, v),
        indeg_now = indeg_next,
        used = used_next
      )
    }
  }
  
  backtrack(
    current_order = integer(0),
    indeg_now = indeg,
    used = rep(FALSE, p)
  )
  
  out
}

simulate_fixed_zio <- function(n, E, alpha_true, beta_true, delta_true, gamma_true, topo) {
  p <- nrow(E)
  X <- matrix(0, nrow = n, ncol = p)
  
  for (j in topo) {
    eta_p <- as.vector(X %*% alpha_true[j, ] + delta_true[j])
    eta_l <- as.vector(X %*% beta_true[j, ] + gamma_true[j])
    
    p_bin  <- plogis(eta_p)
    lambda <- exp(eta_l)
    
    X[, j] <- rpois(n, lambda) * (1 - rbinom(n, 1, p_bin))
  }
  
  colnames(X) <- paste0("V", 1:p)
  X
}

fit_graph_given_order <- function(X, ordering, c0 = 1, criterion = "BIC", dist = "poisson") {
  fit <- zip_stepwise_chain(
    X = X,
    topo = ordering,
    c0 = c0,
    direction = "both",
    criterion = criterion,
    dist = dist,
    trace = FALSE,
    parallel = FALSE
  )
  
  list(
    adjacency = fit$adjacency,
    IC = fit$IC,
    score_sum = sum(fit$IC, na.rm = TRUE),
    PEN = -sum(fit$IC, na.rm = TRUE)
  )
}


############################################################
## 2. Sampling helpers
############################################################

sample_valid_invalid_orders <- function(L_G_mat, E_true, n_valid = 500, n_invalid = 500) {
  p <- nrow(E_true)
  
  ## ---- sample valid orders ----
  n_valid_use <- min(n_valid, nrow(L_G_mat))
  valid_idx <- sample(seq_len(nrow(L_G_mat)), size = n_valid_use, replace = FALSE)
  valid_mat <- L_G_mat[valid_idx, , drop = FALSE]
  
  ## ---- sample invalid orders ----
  invalid_list <- vector("list", n_invalid)
  invalid_strings <- character(0)
  count <- 0L
  
  while (count < n_invalid) {
    ord <- sample.int(p, size = p, replace = FALSE)
    ord_string <- paste(ord, collapse = "-")
    
    if (!is_valid_order_E(ord, E_true) && !(ord_string %in% invalid_strings)) {
      count <- count + 1L
      invalid_list[[count]] <- ord
      invalid_strings <- c(invalid_strings, ord_string)
    }
  }
  
  invalid_mat <- do.call(rbind, invalid_list)
  
  sampled_orders <- rbind(valid_mat, invalid_mat)
  valid_flag <- c(rep(TRUE, nrow(valid_mat)), rep(FALSE, nrow(invalid_mat)))
  
  list(
    orders = sampled_orders,
    valid = valid_flag
  )
}

############################################################
## 3. One-run evaluation on sampled valid/invalid orders
############################################################

one_run_theorem_sampled <- function(n,
                                    E_true,
                                    alpha_true,
                                    beta_true,
                                    delta_true,
                                    gamma_true,
                                    true_topo,
                                    L_G_mat,
                                    n_valid = 500,
                                    n_invalid = 500,
                                    c0 = 1,
                                    criterion = "BIC",
                                    dist = "poisson") {
  
  X <- simulate_fixed_zio(
    n = n,
    E = E_true,
    alpha_true = alpha_true,
    beta_true = beta_true,
    delta_true = delta_true,
    gamma_true = gamma_true,
    topo = true_topo
  )
  
  sampled <- sample_valid_invalid_orders(
    L_G_mat = L_G_mat,
    E_true = E_true,
    n_valid = n_valid,
    n_invalid = n_invalid
  )
  
  ord_mat <- sampled$orders
  valid_vec <- sampled$valid
  
  adj_true_mat <- unname(as.matrix(t(E_true)))
  
  out_list <- vector("list", nrow(ord_mat))
  
  for (i in seq_len(nrow(ord_mat))) {
    ord <- ord_mat[i, ]
    fit <- fit_graph_given_order(
      X = X,
      ordering = ord,
      c0 = c0,
      criterion = criterion,
      dist = dist
    )
    
    adj_fit_mat <- unname(as.matrix(fit$adjacency))
    
    recovered_exact <- isTRUE(all.equal(
      adj_fit_mat,
      adj_true_mat,
      check.attributes = FALSE
    ))
    
    out_list[[i]] <- data.frame(
      n = n,
      ordering = paste(ord, collapse = "-"),
      valid = valid_vec[i],
      PEN = fit$PEN,
      IC_sum = fit$score_sum,
      recovered_exact = recovered_exact,
      edge_diff = sum(abs(adj_fit_mat - adj_true_mat)),
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(out_list)
}