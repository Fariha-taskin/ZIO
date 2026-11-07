## =====================================================================
## ZIPBN order-based hill-climbing: zero-inflation simulations
##
## Key change:
##   Use ONE shared load-balanced task queue over (dataset k, position j)
##   fits instead of nested parallelism.
##
## Why this is better on a 128-core node:
##   - no oversubscription from "parallel over K" and then again inside each k
##   - dynamic load balancing automatically mixes heavy and light regressions
##   - the same pool is reused for every hill-climbing iteration
##   - only the affected block is recomputed after an R2R move
## =====================================================================

# ----------------------------------------------------------------------
# 1. Project paths and required helpers
# ----------------------------------------------------------------------
# From the parent of the repository use ZIO/; from inside ZIO/ use ./.
project_dir <- "ZIO"
if (!dir.exists(file.path(project_dir, "code", "helper"))) {
  if (dir.exists(file.path("code", "helper"))) {
    project_dir <- "."
  } else {
    stop("Cannot find ZIO/code/helper. Run from ZIO/ or its parent.")
  }
}

helper_dir <- file.path(project_dir, "code", "helper")
output_dir <- file.path(project_dir, "output", "simulation", "zero_proportions")

source(file.path(helper_dir, "helper_zio.R"))
source(file.path(helper_dir, "data_generator_from_order_zip_zeros.R"))
source(file.path(helper_dir, "metric.R"))

# ----------------------------------------------------------------------
# 2. Command-line arguments and output paths
# ----------------------------------------------------------------------
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 5L) {
    stop(paste(
      "Usage: Rscript Simulation_ZIO_zeros.R n p K zeros seed",
      "Example: Rscript Simulation_ZIO_zeros.R 1000 20 1 25% 1",
      sep = "\n"
    ))
  }

  n <- suppressWarnings(as.integer(args[1L]))
  p <- suppressWarnings(as.integer(args[2L]))
  K <- suppressWarnings(as.integer(args[3L]))
  seed <- suppressWarnings(as.integer(args[5L]))
  zero_pct <- suppressWarnings(as.numeric(sub("%", "", trimws(args[4L]), fixed = TRUE)))

  if (anyNA(c(n, p, K, seed)) || n < 1L || p < 2L || K < 1L || seed < 1L) {
    stop("n, p, K and seed must be valid positive integers (p >= 2).")
  }
  if (is.na(zero_pct) || !is.finite(zero_pct) || zero_pct < 0 || zero_pct > 100) {
    stop("zeros must be a percentage, e.g. 25%, 50% or 75%.")
  }
  zeros <- paste0(format(zero_pct, trim = TRUE, scientific = FALSE), "%")

  N <- 3500L
  requested_workers <- 20L
  checkpoint_every <- 100L
  leave_free <- 2L

  top_res_dir <- file.path(output_dir, "final_simulation_results")
  log_path <- file.path(
    top_res_dir,
    sprintf("log_result_n_%d_K_%d_P_%d_zeros_%s", n, K, p, zeros)
  )
  final_res_dir <- file.path(
    top_res_dir,
    sprintf("simulation_result_n_%d_K_%d_P_%d_zeros_%s", n, K, p, zeros)
  )
  dir.create(log_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(final_res_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------
## 3. Main loop
## ---------------------------------------------------------------------
nseed <- seed
set.seed(nseed)
    
    ## ---------------------------------------------------------------
    ## 1. Data generation
    ## ---------------------------------------------------------------
    out <- generate_data_from_order(n = n, p = p, K = K, topo = 1:p, zeros = zeros)
    out$dat <- zip_normalize_data_list(out$dat)
    
    K_now <- length(out$dat)
    p_now <- ncol(out$dat[[1L]])
    vars <- colnames(out$dat[[1L]])
    true_edges <- sum(out$E_true[[1L]])
    
    ## ---------------------------------------------------------------
    ## 2. Parallel setup
    ## ---------------------------------------------------------------
    max_init_tasks <- K_now * p_now
    workers <- zip_choose_workers(
      total_tasks = max_init_tasks,
      requested_workers = requested_workers,
      leave_free = leave_free,
      logical = FALSE
    )
    
    cl <- zip_make_cluster(workers = workers, dat_list = out$dat)
    
    ## make sure cluster stops even if error happens later
    on.exit({
      try(zip_stop_cluster(cl), silent = TRUE)
    }, add = TRUE)
    
    ## ---------------------------------------------------------------
    ## 3. Model / search settings
    ## ---------------------------------------------------------------
    c0 <- 1
    criterion <- "BIC"
    dist <- "poisson"
    direction <- "both"
    max_steps <- 1000L
    
    curr_order <- sample.int(p_now, size = p_now, replace = FALSE)
    
    ## ---------------------------------------------------------------
    ## 4. Initial full evaluation
    ## ---------------------------------------------------------------
    start_time <- Sys.time()
    
    init_fit <- zip_evaluate_order(
      dat_list = out$dat,
      topo = curr_order,
      c0 = c0,
      cl = cl,
      direction = direction,
      criterion = criterion,
      dist = dist,
      trace = FALSE,
      max_steps = max_steps
    )
    
    curr_IC <- init_fit$IC_mat
    curr_adj <- init_fit$adjacency
    curr_selected <- init_fit$selected
    curr_bic <- init_fit$total_bic
    
    ## ---------------------------------------------------------------
    ## 5. Storage
    ## ---------------------------------------------------------------
    score <- rep(NA_real_, N)
    BIC_score <- vector("list", N)
    result_IC <- vector("list", N)
    ordering_mat <- matrix(0L, nrow = N, ncol = p_now)
    diff <- vector("list", N)
    Adjacency_matrix <- vector("list", N)
    
    log_file <- file.path(log_path, sprintf("result_p%d_K%d_seed%d.txt", p_now, K_now, nseed))
    cat(
      sprintf("seed: %d | p: %d | K: %d | true_edges(first DAG): %d | workers: %d\n",
              nseed, p_now, K_now, true_edges, workers),
      file = log_file,
      append = TRUE
    )
    
    ## ---------------------------------------------------------------
    ## 6. Hill-climbing over orderings
    ## ---------------------------------------------------------------
    for (t in seq_len(N)) {
      
      cat(sprintf("seed %d | iteration %d / %d\r", nseed, t, N))
      
      iter_start <- Sys.time()
      
      move_pos <- sample.int(p_now, 2L, replace = FALSE)
      
      next_fit <- zip_evaluate_move(
        dat_list = out$dat,
        curr_topo = curr_order,
        pos_i = move_pos[1L],
        pos_j = move_pos[2L],
        curr_IC = curr_IC,
        curr_adj = curr_adj,
        curr_selected = curr_selected,
        c0 = c0,
        cl = cl,
        direction = direction,
        criterion = criterion,
        dist = dist,
        trace = FALSE,
        max_steps = max_steps
      )
      
      next_bic <- next_fit$total_bic
      
      ## greedy accept/reject
      if (next_bic + 1e-8 < curr_bic) {
        curr_order <- next_fit$topo
        curr_IC <- next_fit$IC_mat
        curr_adj <- next_fit$adjacency
        curr_selected <- next_fit$selected
        curr_bic <- next_bic
      }
      
      ## diagnostics
      diff_list <- numeric(K_now)
      FP_k <- numeric(K_now)
      FN_k <- numeric(K_now)
      edge_est_k <- numeric(K_now)
      
      for (k in seq_len(K_now)) {
        adj_est <- t(curr_adj[[k]])
        diff_list[k] <- sum(abs(adj_est - out$E_true[[k]]))
        FP_k[k] <- FP(adj_est, out$E_true[[k]])
        FN_k[k] <- FN(adj_est, out$E_true[[k]])
        edge_est_k[k] <- sum(curr_adj[[k]])
      }
      
      tau <- suppressWarnings(cor(curr_order, 1:p_now, method = "kendall"))
      ave_diff <- mean(diff_list)
      FP_t <- mean(FP_k)
      FN_t <- mean(FN_k)
      edge_est <- mean(edge_est_k)
      iter_seconds <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))
      
      ## store current state
      score[t] <- curr_bic
      BIC_score[[t]] <- curr_bic
      result_IC[[t]] <- curr_IC
      ordering_mat[t, ] <- curr_order
      diff[[t]] <- diff_list
      Adjacency_matrix[[t]] <- curr_adj
      
      ## checkpoint
      if (t %% checkpoint_every == 0L || t == N) {
        checkpoint <- list(
          seed = nseed,
          iteration = t,
          N = N,
          p = p_now,
          K = K_now,
          workers = workers,
          criterion = criterion,
          dist = dist,
          c0 = c0,
          current_order = curr_order,
          current_BIC = curr_bic,
          tau = tau,
          ave_diff = ave_diff,
          FP = FP_t,
          FN = FN_t,
          score = score[seq_len(t)],
          BIC_score = BIC_score[seq_len(t)],
          ordering_mat = ordering_mat[seq_len(t), , drop = FALSE],
          diff = diff[seq_len(t)],
          Adjacency_matrix = Adjacency_matrix[seq_len(t)],
          result_IC = result_IC[seq_len(t)]
        )
        
        checkpoint_file <- file.path(
          log_path,
          sprintf("ZIPBN_checkpoint_P_%d_K_%d_seed%d_iter%d.RData", p_now, K_now, nseed, t)
        )
        save(checkpoint, file = checkpoint_file)
      }
      
      if (t %% 100L == 0L || t == N) {
        cat(
          sprintf(
            "iter %4d / %4d | move=(%d,%d) | affected=%d | time=%7.2f s | tau=%6.3f | FP=%8.3f | FN=%8.3f | Delta=%8.3f | est_edges=%8.3f | BIC=%12.4f\n",
            t, N,
            move_pos[1L], move_pos[2L],
            length(next_fit$affected_positions),
            iter_seconds,
            tau, FP_t, FN_t, ave_diff, edge_est, curr_bic
          ),
          file = log_file,
          append = TRUE
        )
      }
    }
    
    elapsed <- Sys.time() - start_time
    
    ## ---------------------------------------------------------------
    ## 7. Final result object
    ## ---------------------------------------------------------------
    ZIPBN_DAG <- list(
      workers = workers,
      Adjacency_matrix = Adjacency_matrix[[N]],
      result_IC = result_IC[[N]],
      BIC_score = BIC_score[[N]],
      diff = diff,
      ave_diff = ave_diff,
      sampled_ordering = ordering_mat[N, ],
      timeit = elapsed,
      out = out,
      tau = tau,
      FN = FN_t,
      FP = FP_t
    )
    
    final_rdata <- file.path(final_res_dir, sprintf("p%d_K%d_seed%d.RData", p_now, K_now, nseed))
    save(ZIPBN_DAG, file = final_rdata)
    
    summary_file <- file.path(final_res_dir, sprintf("simulation_result_%d_%d_%d.txt", n, p_now, K_now))
    cat(sprintf("\nsim result for seed: %d", nseed), file = summary_file, append = TRUE)
    cat(sprintf("\niteration number: %d", N), file = summary_file, append = TRUE)
    cat("\npermutation:", file = summary_file, append = TRUE)
    cat(ZIPBN_DAG$sampled_ordering, sep = " ", file = summary_file, append = TRUE)
    cat("\ndifference between true DAGs:", file = summary_file, append = TRUE)
    cat(ZIPBN_DAG$ave_diff, sep = " ", file = summary_file, append = TRUE)
    cat(sprintf("\ntau: %s", ZIPBN_DAG$tau), file = summary_file, append = TRUE)
    cat(sprintf("\ntime: %s", ZIPBN_DAG$timeit), file = summary_file, append = TRUE)
    cat(sprintf("\ntrajectory: %s", paste(unlist(ZIPBN_DAG$BIC_score), collapse = " ")),
        file = summary_file, append = TRUE)
    cat(sprintf("\nFP: %s", ZIPBN_DAG$FP), file = summary_file, append = TRUE)
    cat(sprintf("\nFN: %s", ZIPBN_DAG$FN), file = summary_file, append = TRUE)
    cat(sprintf("\nworkers: %d\n", workers), file = summary_file, append = TRUE)
}

main()
