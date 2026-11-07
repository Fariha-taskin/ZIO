## =====================================================================
## ZIPBN fixed-tau ordering violation simulation
##
## Command line:
##   Rscript Simulation_ZIO_fixed_tau.R n p K seed tau_target [N] [workers] [checkpoint_every] [leave_free] [output_dir]
##
## Example:
##   Rscript Simulation_ZIO_fixed_tau.R 1000 10 20 1 0.5 500 4 100 2
##
## Required:
##   args[1] = n
##   args[2] = p
##   args[3] = K
##   args[4] = seed
##   args[5] = tau_target
##
## Optional:
##   args[6]  = N
##   args[7]  = workers
##   args[8]  = checkpoint_every
##   args[9]  = leave_free
##   args[10] = output_dir
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

source(file.path(helper_dir, "helper_zio.R"))
source(file.path(helper_dir, "data_generator_from_order_zip_zeros.R"))
source(file.path(helper_dir, "metric.R"))
source(file.path(helper_dir, "order_sampling_for_fixed_tau.R"))

## ---------------------------------------------------------------------
## 0. Command-line arguments
## ---------------------------------------------------------------------
main <- function() {

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5L) {
  stop(
    paste(
      "Usage: Rscript Simulation_ZIO_fixed_tau.R n p K seed tau_target [N] [workers] [checkpoint_every] [leave_free] [output_dir]",
      "Example: Rscript Simulation_ZIO_fixed_tau.R 1000 10 20 1 0.5 500 4 100 2"
    )
  )
}

n <- as.integer(args[1])
p <- as.integer(args[2])
K <- as.integer(args[3])
seed <- as.integer(args[4])
tau_target <- as.numeric(args[5])

if (anyNA(c(n, p, K, seed, tau_target)) ||
    n < 1L || p < 2L || K < 1L || abs(tau_target) > 1) {
  stop("Invalid n, p, K, seed, or tau_target.")
}

N <- if (length(args) >= 6L) as.integer(args[6]) else 500L
requested_workers <- if (length(args) >= 7L) as.integer(args[7]) else NA_integer_
checkpoint_every <- if (length(args) >= 8L) as.integer(args[8]) else 100L
leave_free <- if (length(args) >= 9L) as.integer(args[9]) else 4L

if (length(args) >= 10L) {
  top_res_dir <- args[10]
} else {
  top_res_dir <- file.path(project_dir, "output", "simulation", "fixed_tau", "final_simulation_results_fixed_tau")
}

if (!dir.exists(top_res_dir)) {
  dir.create(top_res_dir, recursive = TRUE)
}

## ---------------------------------------------------------------------
## 1. Data generation
## ---------------------------------------------------------------------
zeros <- "50%"

nseed <- seed
set.seed(nseed)

out <- generate_data_from_order(
  n = n,
  p = p,
  K = K,
  topo = 1:p,
  zeros = zeros
)

## Get true DAG object
if (!is.null(out$E_true)) {
  G <- out$E_true
} else if (!is.null(out$graph)) {
  G <- out$graph
} else {
  stop("Cannot find true DAGs. Expected out$E_true or out$graph.")
}

X <- out$dat

K <- length(X)
p <- ncol(X[[1L]])

vars <- colnames(X[[1L]])
true_edges <- sum(G[[1L]])

## ---------------------------------------------------------------------
## 2. Sample one ordering per dataset k at fixed Kendall tau
## ---------------------------------------------------------------------
perm_bank <- lapply(seq_len(K), function(k) {
  sample_perm_fixed_tau(
    n = p,
    tau = tau_target,
    baseline = 1:p
  )
})

## Inverse permutations
inv_perm <- lapply(seq_len(K), function(k) {
  order(perm_bank[[k]]$perm)
})

## ---------------------------------------------------------------------
## 3. Apply permutation to each DAG and data matrix
## ---------------------------------------------------------------------
G_perm <- vector("list", K)
X_perm <- vector("list", K)

for (k in seq_len(K)) {
  idx <- inv_perm[[k]]

  G_perm[[k]] <- G[[k]][idx, idx]
  X_perm[[k]] <- X[[k]][, idx]
}

## ZIPBN uses the permuted data matrices directly
dat_list <- X_perm
E_true_perm <- G_perm

## Normalize after permutation
dat_list <- zip_normalize_data_list(dat_list)

## Save actual imposed permutations
perm <- vector("list", K)

for (k in seq_len(K)) {
  perm[[k]] <- perm_bank[[k]]$perm
}

perm_mat <- do.call(cbind, perm)
colnames(perm_mat) <- paste0("sigma", seq_len(K))

if (K >= 2) {
  cor_mat_kendall <- round(cor(perm_mat, method = "kendall"), 3)
  mean_upper_tau <- mean(cor_mat_kendall[upper.tri(cor_mat_kendall)])
} else {
  cor_mat_kendall <- NA
  mean_upper_tau <- NA
}

## Replace out$dat and out$E_true by the permuted versions
out_original <- out
out$dat <- dat_list
out$E_true <- E_true_perm
out$graph <- E_true_perm

K <- length(out$dat)
p <- ncol(out$dat[[1L]])
vars <- colnames(out$dat[[1L]])
true_edges <- sum(out$E_true[[1L]])

## ---------------------------------------------------------------------
## 4. Parallel setup: one persistent shared cluster
## ---------------------------------------------------------------------
max_init_tasks <- K * p

workers <- zip_choose_workers(
  total_tasks = max_init_tasks,
  requested_workers = requested_workers,
  leave_free = leave_free,
  logical = FALSE
)

cl <- zip_make_cluster(
  workers = workers,
  dat_list = out$dat
)

on.exit(zip_stop_cluster(cl), add = TRUE)

## ---------------------------------------------------------------------
## 5. Model / search settings
## ---------------------------------------------------------------------
c0 <- 1
criterion <- "BIC"
dist <- "poisson"
direction <- "both"
max_steps <- 1000L

## Random initial order
curr_order <- sample.int(p, size = p, replace = FALSE)

## ---------------------------------------------------------------------
## 6. Initial full evaluation
## ---------------------------------------------------------------------
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
init_FP <- numeric(K)
init_FN <- numeric(K)

for (k in seq_len(K)) {
  adj_est <- t(curr_adj[[k]])

  init_diff[k] <- sum(abs(adj_est - out$E_true[[k]]))
  init_FP[k] <- FP(adj_est, out$E_true[[k]])
  init_FN[k] <- FN(adj_est, out$E_true[[k]])
}

init_ave_diff <- mean(init_diff)

## ---------------------------------------------------------------------
## 7. Storage
## ---------------------------------------------------------------------
score <- rep(NA_real_, N)
BIC_score <- vector("list", N)
result_IC <- vector("list", N)
ordering_mat <- matrix(0L, nrow = N, ncol = p)
diff <- vector("list", N)
Adjacency_matrix <- vector("list", N)

FP_store <- rep(NA_real_, N)
FN_store <- rep(NA_real_, N)
tau_store <- rep(NA_real_, N)
edge_est_store <- rep(NA_real_, N)

## ---------------------------------------------------------------------
## 8. Output setup
## ---------------------------------------------------------------------
log_path <- file.path(
  top_res_dir,
  sprintf("log_result_K_%d_P_%d_tau_%s", K, p, tau_target)
)

dir.create(log_path, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(
  log_path,
  sprintf("result_p%d_K%d_seed%d_tau%s.txt", p, K, nseed, tau_target)
)

cat(
  sprintf(
    "seed: %d | n: %d | p: %d | K: %d | tau_target: %.3f | mean_upper_tau: %s | true_edges(first DAG): %d | workers: %d\n",
    nseed, n, p, K, tau_target, mean_upper_tau, true_edges, workers
  ),
  file = log_file,
  append = TRUE
)

cat("\nInitial ordering:\n", file = log_file, append = TRUE)
cat(curr_order, sep = " ", file = log_file, append = TRUE)

cat(
  sprintf(
    "\nInitial Delta: %.4f | Initial FP: %.4f | Initial FN: %.4f | Initial BIC: %.4f\n\n",
    init_ave_diff,
    mean(init_FP),
    mean(init_FN),
    curr_bic
  ),
  file = log_file,
  append = TRUE
)

## ---------------------------------------------------------------------
## 9. Hill-climbing over orderings
## ---------------------------------------------------------------------
for (t in seq_len(N)) {
  cat(sprintf("iteration %d / %d\r", t, N))

  iter_start <- Sys.time()

  ## Propose one R2R move in POSITION space
  move_pos <- sample.int(p, 2L, replace = FALSE)

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

  ## Greedy accept/reject
  if (next_bic + 1e-8 < curr_bic) {
    curr_order <- next_fit$topo
    curr_IC <- next_fit$IC_mat
    curr_adj <- next_fit$adjacency
    curr_selected <- next_fit$selected
    curr_bic <- next_bic
  }

  ## Diagnostics
  diff_list <- numeric(K)
  FP_k <- numeric(K)
  FN_k <- numeric(K)
  edge_est_k <- numeric(K)

  for (k in seq_len(K)) {
    adj_est <- t(curr_adj[[k]])

    diff_list[k] <- sum(abs(adj_est - out$E_true[[k]]))
    FP_k[k] <- FP(adj_est, out$E_true[[k]])
    FN_k[k] <- FN(adj_est, out$E_true[[k]])
    edge_est_k[k] <- sum(adj_est)
  }

  tau <- suppressWarnings(cor(curr_order, 1:p, method = "kendall"))

  ave_diff <- mean(diff_list)
  FP_t <- mean(FP_k)
  FN_t <- mean(FN_k)
  edge_est <- mean(edge_est_k)

  iter_seconds <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))

  ## Store trajectory state
  score[t] <- curr_bic
  BIC_score[[t]] <- curr_bic
  result_IC[[t]] <- curr_IC
  ordering_mat[t, ] <- curr_order
  diff[[t]] <- diff_list
  Adjacency_matrix[[t]] <- curr_adj

  FP_store[t] <- FP_t
  FN_store[t] <- FN_t
  tau_store[t] <- tau
  edge_est_store[t] <- edge_est

  ## Periodic checkpoint
  if (t %% checkpoint_every == 0L || t == N) {
    checkpoint <- list(
      seed = nseed,
      n = n,
      iteration = t,
      N = N,
      p = p,
      K = K,
      tau_target = tau_target,
      mean_upper_tau = mean_upper_tau,
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
      result_IC = result_IC[seq_len(t)],
      FP_store = FP_store[seq_len(t)],
      FN_store = FN_store[seq_len(t)],
      tau_store = tau_store[seq_len(t)],
      edge_est_store = edge_est_store[seq_len(t)],
      perm_bank = perm_bank,
      perm_mat = perm_mat,
      cor_mat_kendall = cor_mat_kendall
    )

    checkpoint_file <- file.path(
      log_path,
      sprintf(
        "ZIPBN_checkpoint_P_%d_K_%d_seed%d_tau%s_iter%d.RData",
        p, K, nseed, tau_target, t
      )
    )

    save(checkpoint, file = checkpoint_file)
  }

  ## Periodic log
  if (t %% 100L == 0L || t == N) {
    affected_len <- if (!is.null(next_fit$affected_positions)) {
      length(next_fit$affected_positions)
    } else {
      NA_integer_
    }

    cat(
      sprintf(
        "iter %4d / %4d | move=(%d,%d) | affected=%s | time=%7.2f s | tau=%6.3f | FP=%8.3f | FN=%8.3f | Delta=%8.3f | est_edges=%8.3f | BIC=%12.4f\n",
        t, N,
        move_pos[1L], move_pos[2L],
        affected_len,
        iter_seconds,
        tau, FP_t, FN_t, ave_diff, edge_est, curr_bic
      ),
      file = log_file,
      append = TRUE
    )
  }
}

elapsed <- Sys.time() - start_time

## ---------------------------------------------------------------------
## 10. Final result object
## ---------------------------------------------------------------------
ZIPBN_DAG <- list(
  seed = nseed,
  n = n,
  p = p,
  K = K,
  tau_target = tau_target,
  mean_upper_tau = mean_upper_tau,
  zeros = zeros,

  workers = workers,
  criterion = criterion,
  dist = dist,
  direction = direction,
  c0 = c0,
  max_steps = max_steps,

  Adjacency_matrix = Adjacency_matrix[[N]],
  result_IC = result_IC[[N]],
  BIC_score = BIC_score[[N]],
  score = score,

  diff = diff,
  ave_diff = ave_diff,
  sampled_ordering = ordering_mat[N, ],
  ordering_mat = ordering_mat,

  timeit = elapsed,

  out = out,
  out_original = out_original,

  G_perm = G_perm,
  X_perm = X_perm,
  perm_bank = perm_bank,
  perm_mat = perm_mat,
  cor_mat_kendall = cor_mat_kendall,

  tau = tau,
  FN = FN_t,
  FP = FP_t,
  edge_est = edge_est,

  FP_store = FP_store,
  FN_store = FN_store,
  tau_store = tau_store,
  edge_est_store = edge_est_store
)

## ---------------------------------------------------------------------
## 11. Save final result
## ---------------------------------------------------------------------
final_res_dir <- file.path(
  top_res_dir,
  sprintf("simulation_result_K_%d_P_%d_tau_%s", K, p, tau_target)
)

dir.create(final_res_dir, showWarnings = FALSE, recursive = TRUE)

final_rdata <- file.path(
  final_res_dir,
  sprintf("p%d_K%d_seed%d_tau%s.RData", p, K, nseed, tau_target)
)

save(ZIPBN_DAG, file = final_rdata)

summary_file <- file.path(
  final_res_dir,
  sprintf("simulation_result_%d_%d_%d_tau%s.txt", n, p, K, tau_target)
)

cat(sprintf("\nsim result for seed: %d", nseed), file = summary_file, append = TRUE)
cat(sprintf("\niteration number: %d", N), file = summary_file, append = TRUE)
cat(sprintf("\ntau_target: %s", tau_target), file = summary_file, append = TRUE)
cat(sprintf("\nmean imposed pairwise tau: %s", mean_upper_tau), file = summary_file, append = TRUE)

cat("\npermutation:", file = summary_file, append = TRUE)
cat(ZIPBN_DAG$sampled_ordering, sep = " ", file = summary_file, append = TRUE)

cat("\ndifference between true DAGs:", file = summary_file, append = TRUE)
cat(ZIPBN_DAG$ave_diff, sep = " ", file = summary_file, append = TRUE)

cat(sprintf("\ntau: %s", ZIPBN_DAG$tau), file = summary_file, append = TRUE)
cat(sprintf("\ntime: %s", ZIPBN_DAG$timeit), file = summary_file, append = TRUE)
cat(sprintf("\ntrajectory: %s", ZIPBN_DAG$BIC_score), file = summary_file, append = TRUE)
cat(sprintf("\nFP: %s", ZIPBN_DAG$FP), file = summary_file, append = TRUE)
cat(sprintf("\nFN: %s", ZIPBN_DAG$FN), file = summary_file, append = TRUE)
cat(sprintf("\nworkers: %d\n", workers), file = summary_file, append = TRUE)

cat("\nFinished successfully.\n", file = log_file, append = TRUE)

}

main()
