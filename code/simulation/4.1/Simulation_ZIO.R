## =====================================================================
## ZIO order-based hill-climbing
## Single global worker-pool version
##
## Required command-line arguments:
##   args[1] = p
##   args[2] = K
##   args[3] = seed
##
## Optional command-line arguments:
##   args[4] = N
##   args[5] = workers
##   args[6] = checkpoint_every
##   args[7] = leave_free
##   args[8] = output_dir
## =====================================================================

# Simulation 4.1: order-based ZIO
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

source(file.path(helper_dir, "helper_zio.R"))
source(file.path(helper_dir, "data_generator_for_given_order.R"))
source(file.path(helper_dir, "metric.R"))


main <- function() {

  # -------------------------------------------------------------------
  # 0. Command-line arguments
  # -------------------------------------------------------------------
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 3L) {
    stop(
      paste(
        "Usage:",
        "Rscript Simulation_ZIO.R",
        "p K seed [N] [workers] [checkpoint_every] [leave_free] [output_dir]",
        "\nExample:",
        "Rscript Simulation_ZIO.R 40 20 1 3000 120 100 4 ZIO/output/simulation/4.1/final_simulation_results"
      )
    )
  }

  p <- suppressWarnings(as.integer(args[1]))
  K <- suppressWarnings(as.integer(args[2]))
  seed <- suppressWarnings(as.integer(args[3]))

  N <- if (length(args) >= 4L) {
    suppressWarnings(as.integer(args[4]))
  } else {
    30L
  }

  requested_workers <- if (length(args) >= 5L) {
    suppressWarnings(as.integer(args[5]))
  } else {
    4L
  }

  checkpoint_every <- if (length(args) >= 6L) {
    suppressWarnings(as.integer(args[6]))
  } else {
    100L
  }

  leave_free <- if (length(args) >= 7L) {
    suppressWarnings(as.integer(args[7]))
  } else {
    1L
  }

  top_res_dir <- if (length(args) >= 8L) {
    args[8]
  } else {
    file.path(
      output_dir,
      "final_simulation_results"
    )
  }

  # Validate inputs
  if (is.na(p) || p < 2L) {
    stop("p must be an integer greater than or equal to 2.")
  }

  if (is.na(K) || K < 1L) {
    stop("K must be a positive integer.")
  }

  if (is.na(seed) || seed < 1L) {
    stop("seed must be a positive integer.")
  }

  if (is.na(N) || N < 1L) {
    stop("N must be a positive integer.")
  }

  if (is.na(requested_workers) || requested_workers < 1L) {
    stop("workers must be a positive integer.")
  }

  if (is.na(checkpoint_every) || checkpoint_every < 1L) {
    stop("checkpoint_every must be a positive integer.")
  }

  if (is.na(leave_free) || leave_free < 0L) {
    stop("leave_free must be a nonnegative integer.")
  }

  dir.create(
    top_res_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  cat("===== Simulation settings =====\n")
  cat("p:                  ", p, "\n")
  cat("K:                  ", K, "\n")
  cat("fixed seed:         ", seed, "\n")
  cat("iterations:         ", N, "\n")
  cat("requested workers:  ", requested_workers, "\n")
  cat("checkpoint interval:", checkpoint_every, "\n")
  cat("leave free:         ", leave_free, "\n")
  cat("output directory:   ", top_res_dir, "\n")
  cat("===============================\n")


  # -------------------------------------------------------------------
  # 1. Data generation
  # -------------------------------------------------------------------
  nseed <- seed
  set.seed(nseed)

  out <- generate_data_from_order(
    n = 1000,
    p = p,
    K = K,
    topo = seq_len(p)
  )

  out$dat <- zip_normalize_data_list(out$dat)

  K <- length(out$dat)
  p <- ncol(out$dat[[1L]])
  vars <- colnames(out$dat[[1L]])
  true_edges <- sum(out$E_true[[1L]])


  # -------------------------------------------------------------------
  # 2. Parallel setup: one persistent shared cluster
  # -------------------------------------------------------------------
  max_init_tasks <- K * p

  workers <- zip_choose_workers(
    total_tasks = max_init_tasks,
    requested_workers = requested_workers,
    leave_free = leave_free,
    logical = FALSE
  )

  cat("Workers actually used:", workers, "\n")

  cl <- NULL

  on.exit(
    {
      if (!is.null(cl)) {
        zip_stop_cluster(cl)
      }
    },
    add = TRUE
  )

  cl <- zip_make_cluster(
    workers = workers,
    dat_list = out$dat
  )


  # -------------------------------------------------------------------
  # 3. Model and search settings
  # -------------------------------------------------------------------
  c0 <- 1
  criterion <- "BIC"
  dist <- "poisson"
  direction <- "both"
  max_steps <- 1000L

  # Random initial order
  curr_order <- sample.int(
    n = p,
    size = p,
    replace = FALSE
  )


  # -------------------------------------------------------------------
  # 4. Initial full evaluation
  # -------------------------------------------------------------------
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

  init_diff <- numeric(K)

  for (k in seq_len(K)) {
    init_diff[k] <- sum(
      abs(t(curr_adj[[k]]) - out$E_true[[k]])
    )
  }


  # -------------------------------------------------------------------
  # 5. Storage
  # -------------------------------------------------------------------
  score <- rep(NA_real_, N)
  BIC_score <- vector("list", N)
  result_IC <- vector("list", N)

  ordering_mat <- matrix(
    0L,
    nrow = N,
    ncol = p
  )

  diff <- vector("list", N)
  Adjacency_matrix <- vector("list", N)

  log_path <- file.path(
    top_res_dir,
    sprintf("log_result_K_%d_P_%d", K, p)
  )

  dir.create(
    log_path,
    showWarnings = FALSE,
    recursive = TRUE
  )

  log_file <- file.path(
    log_path,
    sprintf(
      "result_p%d_K%d_seed%d.txt",
      p,
      K,
      nseed
    )
  )

  cat(
    sprintf(
      paste0(
        "seed: %d | p: %d | K: %d | ",
        "true_edges(first DAG): %d | workers: %d\n"
      ),
      nseed,
      p,
      K,
      true_edges,
      workers
    ),
    file = log_file,
    append = FALSE
  )


  # -------------------------------------------------------------------
  # 6. Hill-climbing over orderings
  # -------------------------------------------------------------------
  for (t in seq_len(N)) {

    cat(sprintf("iteration %d / %d\r", t, N))
    flush.console()

    iter_start <- Sys.time()

    # Propose one R2R move in position space
    move_pos <- sample.int(
      n = p,
      size = 2L,
      replace = FALSE
    )

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

    # Greedy accept/reject
    if (next_bic + 1e-8 < curr_bic) {
      curr_order <- next_fit$topo
      curr_IC <- next_fit$IC_mat
      curr_adj <- next_fit$adjacency
      curr_selected <- next_fit$selected
      curr_bic <- next_bic
    }

    # Diagnostics
    diff_list <- numeric(K)
    FP_k <- numeric(K)
    FN_k <- numeric(K)
    edge_est_k <- numeric(K)

    for (k in seq_len(K)) {
      adj_est <- t(curr_adj[[k]])

      diff_list[k] <- sum(
        abs(adj_est - out$E_true[[k]])
      )

      FP_k[k] <- FP(
        adj_est,
        out$E_true[[k]]
      )

      FN_k[k] <- FN(
        adj_est,
        out$E_true[[k]]
      )

      edge_est_k[k] <- sum(curr_adj[[k]])
    }

    tau <- suppressWarnings(
      cor(
        curr_order,
        seq_len(p),
        method = "kendall"
      )
    )

    ave_diff <- mean(diff_list)
    FP_t <- mean(FP_k)
    FN_t <- mean(FN_k)
    edge_est <- mean(edge_est_k)

    iter_seconds <- as.numeric(
      difftime(
        Sys.time(),
        iter_start,
        units = "secs"
      )
    )

    # Store trajectory
    score[t] <- curr_bic
    BIC_score[[t]] <- curr_bic
    result_IC[[t]] <- curr_IC
    ordering_mat[t, ] <- curr_order
    diff[[t]] <- diff_list
    Adjacency_matrix[[t]] <- curr_adj

    # Periodic checkpoint
    if (t %% checkpoint_every == 0L || t == N) {

      checkpoint <- list(
        seed = nseed,
        iteration = t,
        N = N,
        p = p,
        K = K,
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
        ordering_mat = ordering_mat[
          seq_len(t),
          ,
          drop = FALSE
        ],
        diff = diff[seq_len(t)],
        Adjacency_matrix = Adjacency_matrix[seq_len(t)],
        result_IC = result_IC[seq_len(t)]
      )

      checkpoint_file <- file.path(
        log_path,
        sprintf(
          "ZIO_checkpoint_P_%d_K_%d_seed%d_iter%d.RData",
          p,
          K,
          nseed,
          t
        )
      )

      save(
        checkpoint,
        file = checkpoint_file
      )
    }

    # Progress log
    if (t %% 100L == 0L || t == N) {

      cat(
        sprintf(
          paste0(
            "iter %4d / %4d | move=(%d,%d) | affected=%d | ",
            "time=%7.2f s | tau=%6.3f | FP=%8.3f | ",
            "FN=%8.3f | Delta=%8.3f | est_edges=%8.3f | ",
            "BIC=%12.4f\n"
          ),
          t,
          N,
          move_pos[1L],
          move_pos[2L],
          length(next_fit$affected_positions),
          iter_seconds,
          tau,
          FP_t,
          FN_t,
          ave_diff,
          edge_est,
          curr_bic
        ),
        file = log_file,
        append = TRUE
      )
    }
  }

  cat("\n")

  elapsed <- Sys.time() - start_time


  # -------------------------------------------------------------------
  # 7. Final result object
  # -------------------------------------------------------------------
  ZIO_DAG <- list(
    seed = nseed,
    p = p,
    K = K,
    N = N,
    workers = workers,
    Adjacency_matrix = Adjacency_matrix[[N]],
    result_IC = result_IC[[N]],
    BIC_score = score,
    diff = diff,
    ave_diff = ave_diff,
    sampled_ordering = ordering_mat[N, ],
    timeit = elapsed,
    out = out,
    tau = tau,
    FN = FN_t,
    FP = FP_t
  )

  final_res_dir <- file.path(
    top_res_dir,
    sprintf(
      "simulation_result_K_%d_P_%d",
      K,
      p
    )
  )

  dir.create(
    final_res_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  final_rdata <- file.path(
    final_res_dir,
    sprintf(
      "p%d_K%d_seed%d.RData",
      p,
      K,
      nseed
    )
  )

  save(
    ZIO_DAG,
    file = final_rdata
  )

  summary_file <- file.path(
    final_res_dir,
    sprintf(
      "simulation_result_%d_%d_%d.txt",
      1000,
      p,
      K
    )
  )

  cat(
    sprintf("Simulation result for seed: %d\n", nseed),
    file = summary_file,
    append = FALSE
  )

  cat(
    sprintf("Iteration number: %d\n", N),
    file = summary_file,
    append = TRUE
  )

  cat(
    paste0(
      "Permutation: ",
      paste(ZIO_DAG$sampled_ordering, collapse = " "),
      "\n"
    ),
    file = summary_file,
    append = TRUE
  )

  cat(
    sprintf(
      "Average difference between true and estimated DAGs: %s\n",
      ZIO_DAG$ave_diff
    ),
    file = summary_file,
    append = TRUE
  )

  cat(
    sprintf("Tau: %s\n", ZIO_DAG$tau),
    file = summary_file,
    append = TRUE
  )

  cat(
    sprintf("Time: %s\n", ZIO_DAG$timeit),
    file = summary_file,
    append = TRUE
  )

  cat(
    paste0(
      "BIC trajectory: ",
      paste(ZIO_DAG$BIC_score, collapse = " "),
      "\n"
    ),
    file = summary_file,
    append = TRUE
  )

  cat(
    sprintf("FP: %s\n", ZIO_DAG$FP),
    file = summary_file,
    append = TRUE
  )

  cat(
    sprintf("FN: %s\n", ZIO_DAG$FN),
    file = summary_file,
    append = TRUE
  )

  cat(
    sprintf("Workers: %d\n", workers),
    file = summary_file,
    append = TRUE
  )

  cat("Simulation completed successfully.\n")
  cat("Final RData:", final_rdata, "\n")
  cat("Summary file:", summary_file, "\n")
}


main()