# Joint ZIO real-data analysis: one final .RData file
# Script location: ZIO/code/real_data/5/

# ------------------------------------------------------------
# Paths and settings
# ------------------------------------------------------------

project_dir <- "ZIO"
if (!dir.exists(file.path(project_dir, "code", "helper"))) {
  if (dir.exists(file.path("code", "helper"))) {
    project_dir <- "."
  } else {
    stop("Cannot locate ZIO/code/helper. Run from ZIO/ or its parent directory.")
  }
}
helper_file <- file.path(project_dir, "code", "helper", "helper_zio.R")
data_dir <- file.path(project_dir, "data")
res_dir <- file.path(project_dir, "output")

p <- 41L
seed <- 1L
N <- 3000L
c0 <- 1
model_distribution <- "poisson"
remove_original_index <- 14L
parallel_min_tasks <- 4L

control_file <- file.path(data_dir, paste0("control_list_", p, ".rds"))
case_file <- file.path(data_dir, paste0("suicide_list_", p, ".rds"))

# Paths are absolute, so the script works from any R working directory.
for (path in c(helper_file, control_file, case_file)) {
  if (!file.exists(path)) stop("Required file not found: ", path)
}

library(pscl)
library(MASS)
source(helper_file)

# Prevent excessive threading within parallel workers.
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)

# ------------------------------------------------------------
# Load and prepare datasets
# ------------------------------------------------------------

control <- readRDS(control_file)
case <- readRDS(case_file)
stopifnot(is.list(control), is.list(case))

fill_names <- function(x, prefix) {
  nm <- names(x)
  if (is.null(nm)) nm <- rep("", length(x))
  missing <- is.na(nm) | !nzchar(trimws(nm))
  nm[missing] <- paste0(prefix, which(missing))
  nm
}

names(control) <- fill_names(control, "control_")
names(case) <- fill_names(case, "case_")

X <- c(unname(control), unname(case))
dataset_info <- data.frame(
  original_index = seq_along(X),
  group = c(rep("control", length(control)), rep("case", length(case))),
  dataset_name = c(names(control), names(case)),
  stringsAsFactors = FALSE
)

if (remove_original_index > length(X) || remove_original_index < 1L) {
  stop("remove_original_index is outside the available datasets.")
}
X <- X[-remove_original_index]
dataset_info <- dataset_info[-remove_original_index, , drop = FALSE]
dataset_info$dataset_id <- seq_along(X)
K <- length(X)
if (K == 0L) stop("No datasets remain.")

prepare_matrix <- function(x, k) {
  x <- as.matrix(x)
  if (ncol(x) != p && nrow(x) == p) x <- t(x)
  if (ncol(x) != p) stop("Dataset ", k, ": expected ", p, " genes.")
  storage.mode(x) <- "numeric"
  x
}
X <- lapply(seq_along(X), function(k) prepare_matrix(X[[k]], k))

# Align all datasets to the gene ordering in the first dataset.
clean_genes <- function(x) {
  nm <- colnames(x)
  if (is.null(nm)) nm <- paste0("V", seq_len(p))
  missing <- is.na(nm) | !nzchar(trimws(nm))
  nm[missing] <- paste0("V", which(missing))
  if (anyDuplicated(nm)) stop("Duplicated gene names found.")
  nm
}

genes <- clean_genes(X[[1L]])
safe_genes <- make.names(genes, unique = TRUE)

for (k in seq_len(K)) {
  current_genes <- clean_genes(X[[k]])
  if (!setequal(current_genes, genes)) {
    stop("Dataset ", k, " does not contain the same genes as dataset 1.")
  }
  X[[k]] <- X[[k]][, match(genes, current_genes), drop = FALSE]
  colnames(X[[k]]) <- safe_genes
  X[[k]] <- as.data.frame(X[[k]], check.names = FALSE)
}

names(X) <- dataset_info$dataset_name
n_sample <- vapply(X, nrow, integer(1))
dataset_info$sample_size <- n_sample
rm(control, case)

# ------------------------------------------------------------
# Parallel setup
# ------------------------------------------------------------

main <- function() {
  start_time <- Sys.time()

  n_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
  if (is.na(n_cpus) || n_cpus < 1L) {
    n_cpus <- suppressWarnings(parallel::detectCores(logical = FALSE))
  }
  if (is.na(n_cpus) || n_cpus < 1L) n_cpus <- 1L

  n_workers <- max(1L, min(K * p, n_cpus - 1L))
  cl <- if (n_workers > 1L) {
    zip_make_cluster(workers = n_workers, dat_list = X)
  } else {
    NULL
  }
  on.exit(zip_stop_cluster(cl), add = TRUE)

  # ----------------------------------------------------------
  # Initial joint fit
  # ----------------------------------------------------------

  set.seed(seed)
  current_order <- sample.int(p)
  fit <- zip_evaluate_order(
    dat_list = X, topo = current_order, c0 = c0, cl = cl,
    direction = "both", criterion = "BIC",
    dist = model_distribution, trace = FALSE, max_steps = 1000L
  )

  current_IC <- fit$IC_mat
  current_BIC <- fit$total_bic
  current_adj <- fit$adjacency
  current_selected <- fit$selected

  if (!identical(dim(current_IC), c(K, p)) ||
      length(current_adj) != K || length(current_selected) != K ||
      !is.finite(current_BIC)) {
    stop("Initial joint fit returned invalid results.")
  }

  # Store distinct accepted states; reconstruct every iteration at the end.
  compact_adjacency <- list(current_adj)
  compact_IC <- list(current_IC)
  state_index <- integer(N)
  sampled_ordering <- matrix(NA_integer_, nrow = N, ncol = p)
  BIC_score <- numeric(N)
  accepted <- logical(N)

  # ----------------------------------------------------------
  # Joint random-to-random hill-climbing
  # ----------------------------------------------------------

  for (t in seq_len(N)) {
    idx <- sample.int(p, 2L)
    n_tasks <- K * (abs(diff(idx)) + 1L)
    move_cl <- if (!is.null(cl) && n_tasks >= parallel_min_tasks) cl else NULL

    proposal <- zip_evaluate_move(
      dat_list = X, curr_topo = current_order,
      pos_i = idx[1L], pos_j = idx[2L],
      curr_IC = current_IC, curr_adj = current_adj,
      curr_selected = current_selected, c0 = c0, cl = move_cl,
      direction = "both", criterion = "BIC",
      dist = model_distribution, trace = FALSE, max_steps = 1000L
    )

    if (is.finite(proposal$total_bic) && proposal$total_bic < current_BIC) {
      current_order <- proposal$topo
      current_IC <- proposal$IC_mat
      current_BIC <- proposal$total_bic
      current_adj <- proposal$adjacency
      current_selected <- proposal$selected
      compact_adjacency[[length(compact_adjacency) + 1L]] <- current_adj
      compact_IC[[length(compact_IC) + 1L]] <- current_IC
      accepted[t] <- TRUE
    }

    sampled_ordering[t, ] <- current_order
    BIC_score[t] <- current_BIC
    state_index[t] <- length(compact_adjacency)

    if (t %% 100L == 0L || t == N) {
      cat(sprintf("Iteration %d/%d | BIC = %.6f | accepted = %s\n",
                  t, N, current_BIC, accepted[t]))
    }
  }

  zip_stop_cluster(cl)
  cl <- NULL
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # ----------------------------------------------------------
  # Save the complete results in ONE .RData file
  # ----------------------------------------------------------

  ZIO_DAG <- list(
    K = K, p = p, c0 = c0, distribution = model_distribution,
    dataset_info = dataset_info, n_sample = n_sample,
    removed_index = remove_original_index,
    allocated_cpus = n_cpus, node_workers = n_workers,
    completed = N, target_iterations = N,
    Adjacency_matrix = compact_adjacency[state_index],
    result_IC = compact_IC[state_index],
    BIC_score = BIC_score, score = BIC_score,
    sampled_ordering = sampled_ordering,
    accepted = accepted, state_index = state_index,
    history_mode = "full iteration history",
    final_ordering = current_order, final_IC = current_IC,
    final_BIC = current_BIC, final_adjacency = current_adj,
    E_estimate = current_adj,
    timeit = as.difftime(elapsed, units = "secs"),
    timeit_seconds = elapsed
  )

  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(res_dir, sprintf("ZIO_DAG_P_%d_K_%d.RData", p, K))
  save(ZIO_DAG, file = output_file)
  cat("Saved:", output_file, "\n")
  invisible(ZIO_DAG)
}

main()
