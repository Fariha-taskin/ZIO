## =====================================================================
## ZIO order-based hill-climbing (single global worker pool version)
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
# ------------------------------------------------------------
# Portable ZIO repository paths (run from repo root or parent)
# ------------------------------------------------------------
project_dir <- "ZIO"
if (!dir.exists(file.path(project_dir, "code", "helper")) &&
    dir.exists(file.path("code", "helper"))) {
  project_dir <- "."
}
helper_dir <- file.path(project_dir, "code", "helper")
if (!dir.exists(helper_dir)) {
  stop("Cannot find ZIO/code/helper. Run from the ZIO repository or its parent.")
}

source(file.path(helper_dir, "helper_zio_adj.R"))
source(file.path(helper_dir, "data_generator_for_given_order.R"))
source(file.path(helper_dir, "metric.R"))

# ---------------------------------------------------------------------
# 0. Command-line arguments
# ---------------------------------------------------------------------
# Required:
#   args[1] = p
#   args[2] = K
#   args[3] = seed
#
# Optional:
#   args[4] = N                 (number of hill-climbing iterations)
#   args[5] = workers           (size of shared worker pool)
#   args[6] = checkpoint_every  (save every this many iterations)
#   args[7] = leave_free        (reserve this many cores for the OS)
#   args[8] = output_dir        (top-level result directory)
main <- function() {

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3L) {
  stop(
    paste(
      "Usage: Rscript Simulation_ZIO_fixed_data_ADJ.R p K nseed [N] [workers] [checkpoint_every] [leave_free] [output_dir]",
      "Example: Rscript Simulation_ZIO_fixed_data_ADJ.R 10 20 1 2000 120 100 4 /path/to/results"
    )
  )
}

p <- as.integer(args[1])
K <- as.integer(args[2])
seed <- as.integer(args[3])

if (anyNA(c(p, K, seed)) || p < 2L || K < 1L) {
  stop("p, K, and seed must be valid integers (p >= 2, K >= 1).")
}

N <- if (length(args) >= 4L) as.integer(args[4]) else 3333  #  (10/3)*p^3
requested_workers <- if (length(args) >= 5L) as.integer(args[5]) else NA_integer_
checkpoint_every <- if (length(args) >= 6L) as.integer(args[6]) else 100L
leave_free <- if (length(args) >= 7L) as.integer(args[7]) else 4L

if (length(args) >= 8L) {
  top_res_dir <- args[8]
} else {
  top_res_dir <- file.path(project_dir, "output", "simulation", "operator_comparison", "final_simulation_results_fixed_data_adj")
}

if (!dir.exists(top_res_dir)) {
  dir.create(top_res_dir, recursive = TRUE)
}

# ---------------------------------------------------------------------
# 1. Data generation
# ---------------------------------------------------------------------

set.seed(2026)
out <- generate_data_from_order(n = 1000, p = p, K = K, topo = 1:p)
out$dat <- zip_normalize_data_list(out$dat)

K <- length(out$dat)
p <- ncol(out$dat[[1L]])
vars <- colnames(out$dat[[1L]])
true_edges <- sum(out$E_true[[1L]])

nseed = seed + 30
set.seed(nseed)

# ---------------------------------------------------------------------
# 2. Parallel setup: one persistent shared cluster
# ---------------------------------------------------------------------
max_init_tasks <- K * p
workers <- zip_choose_workers(
  total_tasks = max_init_tasks,
  requested_workers = requested_workers,
  leave_free = leave_free,
  logical = FALSE
)

cl <- zip_make_cluster(workers = workers, dat_list = out$dat)
on.exit(zip_stop_cluster(cl), add = TRUE)

# ---------------------------------------------------------------------
# 3. Model / search settings
# ---------------------------------------------------------------------
c0 <- 1
criterion <- "BIC"
dist <- "poisson"
direction <- "both"
max_steps <- 1000L

# Random initial order
curr_order <- sample.int(p, size = p, replace = FALSE)

# ---------------------------------------------------------------------
# 4. Initial full evaluation
# ---------------------------------------------------------------------
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
  init_diff[k] <- sum(abs(t(curr_adj[[k]]) - out$E_true[[k]]))
}

# ---------------------------------------------------------------------
# 5. Storage
# ---------------------------------------------------------------------
score <- rep(NA_real_, N)
BIC_score <- vector("list", N)
result_IC <- vector("list", N)
ordering_mat <- matrix(0L, nrow = N, ncol = p)
diff <- vector("list", N)
Adjacency_matrix <- vector("list", N)

# lightweight progress log
log_path <- file.path(top_res_dir, sprintf("log_result_K_%d_P_%d", K, p))
dir.create(log_path, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(log_path, sprintf("result_p%d_K%d_seed%d.txt", p, K, nseed))
cat(
  sprintf("seed: %d | p: %d | K: %d | true_edges(first DAG): %d | workers: %d\n",
          nseed, p, K, true_edges, workers),
  file = log_file,
  append = TRUE
)

# ---------------------------------------------------------------------
# 6. Hill-climbing over orderings
# ---------------------------------------------------------------------
for (t in seq_len(N)) {
  cat(sprintf("iteration %d / %d\r", t, N))
  iter_start <- Sys.time()

  #  ADJ move in POSITION space
    move_pos <- sample.int(p - 1L, 1L)

    next_fit <- zip_evaluate_move(
      dat_list = out$dat,
      curr_topo = curr_order,
      pos_i = move_pos,
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
    diff_list[k] <- sum(abs(adj_est - out$E_true[[k]]))
    FP_k[k] <- FP(adj_est, out$E_true[[k]])
    FN_k[k] <- FN(adj_est, out$E_true[[k]])
    edge_est_k[k] <- sum(curr_adj[[k]])
  }

  tau <- suppressWarnings(cor(curr_order, 1:p, method = "kendall"))
  ave_diff <- mean(diff_list)
  FP_t <- mean(FP_k)
  FN_t <- mean(FN_k)
  edge_est <- mean(edge_est_k)
  iter_seconds <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))

  # Store trajectory state
  score[t] <- curr_bic
  BIC_score[[t]] <- curr_bic
  result_IC[[t]] <- curr_IC
  ordering_mat[t, ] <- curr_order
  diff[[t]] <- diff_list
  Adjacency_matrix[[t]] <- curr_adj

  # Periodic checkpoint (avoid saving the full state every iteration)
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
      ordering_mat = ordering_mat[seq_len(t), , drop = FALSE],
      diff = diff[seq_len(t)],
      Adjacency_matrix = Adjacency_matrix[seq_len(t)],
      result_IC = result_IC[seq_len(t)]
    )

    checkpoint_file <- file.path(
      log_path,
      sprintf("ZIO_checkpoint_P_%d_K_%d_seed%d_iter%d.RData", p, K, nseed, t)
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

# ---------------------------------------------------------------------
# 7. Final result object
# ---------------------------------------------------------------------
ZIO_DAG <- list(
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

final_res_dir <- file.path(top_res_dir, sprintf("simulation_result_K_%d_P_%d", K, p))
dir.create(final_res_dir, showWarnings = FALSE, recursive = TRUE)

final_rdata <- file.path(final_res_dir, sprintf("p%d_K%d_seed%d.RData", p, K, nseed))
save(ZIO_DAG, file = final_rdata)

summary_file <- file.path(final_res_dir, sprintf("simulation_result_%d_%d_%d.txt", 1000, p, K))
cat(sprintf("\nsim result for seed: %d", nseed), file = summary_file, append = TRUE)
cat(sprintf("\niteration number: %d", N), file = summary_file, append = TRUE)
cat("\npermutation:", file = summary_file, append = TRUE)
cat(ZIO_DAG$sampled_ordering, sep = " ", file = summary_file, append = TRUE)
cat("\ndifference between true DAGs:", file = summary_file, append = TRUE)
cat(ZIO_DAG$ave_diff, sep = " ", file = summary_file, append = TRUE)
cat(sprintf("\ntau: %s", ZIO_DAG$tau), file = summary_file, append = TRUE)
cat(sprintf("\ntime: %s", ZIO_DAG$timeit), file = summary_file, append = TRUE)
cat(sprintf("\ntrajectory: %s", ZIO_DAG$BIC_score), file = summary_file, append = TRUE)
cat(sprintf("\nFP: %s", ZIO_DAG$FP), file = summary_file, append = TRUE)
cat(sprintf("\nFN: %s", ZIO_DAG$FN), file = summary_file, append = TRUE)
cat(sprintf("\nworkers: %d\n", workers), file = summary_file, append = TRUE)

}

main()
