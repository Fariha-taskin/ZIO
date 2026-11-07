# Independent ZIPBN real-data analysis (K = 1)
# Script location: ZIO/code/supplementary/C2/
# Usage: Rscript ZIO_realdata_individual_server_parallel.R <dataset_id>
# Under Slurm, SLURM_ARRAY_TASK_ID is used if no dataset ID is passed.
# Each run saves ONE .RData file; no checkpoints, logs, or CSVs.

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
K <- 1L
c0 <- 1
model_distribution <- "poisson"
N_target <- 300L                      # Same as the supplied script
base_seed <- 1L
remove_original_index <- 14L
parallel_min_tasks <- 4L

initial_order <- c(
  30, 16, 21, 12, 37, 33, 24, 28, 41, 40, 18, 27, 5, 1, 39, 34,
  3, 11, 31, 13, 7, 6, 23, 29, 25, 8, 36, 26, 38, 15, 14, 2,
  10, 17, 35, 22, 4, 20, 9, 32, 19
)

control_file <- file.path(data_dir, paste0("control_list_", p, ".rds"))
case_file <- file.path(data_dir, paste0("suicide_list_", p, ".rds"))

for (path in c(helper_file, control_file, case_file)) {
  if (!file.exists(path)) stop("Required file not found: ", path)
}
if (!identical(sort(as.integer(initial_order)), seq_len(p))) {
  stop("initial_order must be a permutation of 1:p.")
}

args <- commandArgs(trailingOnly = TRUE)
id_text <- if (length(args)) args[1L] else Sys.getenv("SLURM_ARRAY_TASK_ID", "")
dataset_id <- suppressWarnings(as.integer(id_text))
if (is.na(dataset_id) || !grepl("^[0-9]+$", id_text)) {
  stop("Specify an integer dataset ID: Rscript ZIO_realdata_individual_server_parallel.R 1")
}

suppressPackageStartupMessages({
  library(pscl)
  library(MASS)
})
sys.source(helper_file, envir = .GlobalEnv)

required_functions <- c("zip_make_cluster", "zip_stop_cluster",
                        "zip_evaluate_order", "zip_evaluate_move")
missing_functions <- required_functions[!vapply(required_functions, exists,
  logical(1L), mode = "function", envir = .GlobalEnv, inherits = FALSE)]
if (length(missing_functions)) {
  stop("Missing helper functions: ", paste(missing_functions, collapse = ", "))
}

# Keep each parallel R process single-threaded for BLAS/OpenMP.
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)

# ------------------------------------------------------------
# Load and select one dataset
# ------------------------------------------------------------

control <- readRDS(control_file)
case <- readRDS(case_file)
if (!is.list(control) || !is.list(case)) {
  stop("Control and case RDS objects must both be lists.")
}

fill_names <- function(x, prefix) {
  nm <- names(x)
  if (is.null(nm)) nm <- rep("", length(x))
  missing <- is.na(nm) | !nzchar(trimws(nm))
  nm[missing] <- paste0(prefix, which(missing))
  nm
}

names(control) <- fill_names(control, "control_")
names(case) <- fill_names(case, "case_")

all_data <- c(unname(control), unname(case))
dataset_info <- data.frame(
  original_index = seq_along(all_data),
  group = c(rep("control", length(control)), rep("case", length(case))),
  dataset_name = c(names(control), names(case)),
  stringsAsFactors = FALSE
)

if (remove_original_index < 1L || remove_original_index > length(all_data)) {
  stop("remove_original_index is outside the available datasets.")
}
all_data <- all_data[-remove_original_index]
dataset_info <- dataset_info[-remove_original_index, , drop = FALSE]
dataset_info$dataset_id <- seq_along(all_data)

if (dataset_id < 1L || dataset_id > length(all_data)) {
  stop("dataset_id must be between 1 and ", length(all_data), ".")
}

current_info <- dataset_info[dataset_id, , drop = FALSE]
X <- as.matrix(all_data[[dataset_id]])
if (ncol(X) != p && nrow(X) == p) X <- t(X)
if (ncol(X) != p) {
  stop("Dataset ", dataset_id, " must have ", p, " gene columns.")
}
storage.mode(X) <- "numeric"

if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(p))
missing <- is.na(colnames(X)) | !nzchar(trimws(colnames(X)))
colnames(X)[missing] <- paste0("V", which(missing))
colnames(X) <- make.names(colnames(X), unique = TRUE)
X <- as.data.frame(X, check.names = FALSE)
n_sample <- nrow(X)
dat_list <- list(X)

rm(control, case, all_data, dataset_info)

# ------------------------------------------------------------
# Parallel setup and independent hill-climbing
# ------------------------------------------------------------

main <- function() {
  start_time <- Sys.time()

  n_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
  if (is.na(n_cpus) || n_cpus < 1L) {
    n_cpus <- suppressWarnings(parallel::detectCores(logical = FALSE))
  }
  if (is.na(n_cpus) || n_cpus < 1L) n_cpus <- 1L
  n_workers <- max(1L, min(p, n_cpus - 1L))

  cl <- if (n_workers > 1L) {
    zip_make_cluster(workers = n_workers, dat_list = dat_list)
  } else {
    NULL
  }
  on.exit(zip_stop_cluster(cl), add = TRUE)

  set.seed(base_seed + dataset_id)
  current_order <- initial_order

  fit <- zip_evaluate_order(
    dat_list = dat_list, topo = current_order, c0 = c0, cl = cl,
    direction = "both", criterion = "BIC",
    dist = model_distribution, trace = FALSE, max_steps = 1000L
  )

  current_IC <- fit$IC_mat
  current_BIC <- fit$total_bic
  current_adj <- fit$adjacency[[1L]]
  current_selected <- fit$selected[[1L]]
  if (!is.finite(current_BIC)) stop("Initial fit returned a non-finite BIC.")

  sampled_ordering <- matrix(NA_integer_, nrow = N_target, ncol = p)
  BIC_score <- numeric(N_target)
  Adjacency_matrix <- vector("list", N_target)
  result_IC <- vector("list", N_target)
  accepted <- logical(N_target)

  for (t in seq_len(N_target)) {
    idx <- sample.int(p, 2L)
    n_affected <- abs(idx[1L] - idx[2L]) + 1L
    move_cl <- if (!is.null(cl) && n_affected >= parallel_min_tasks) cl else NULL

    proposal <- zip_evaluate_move(
      dat_list = dat_list, curr_topo = current_order,
      pos_i = idx[1L], pos_j = idx[2L],
      curr_IC = current_IC, curr_adj = list(current_adj),
      curr_selected = list(current_selected), c0 = c0, cl = move_cl,
      direction = "both", criterion = "BIC",
      dist = model_distribution, trace = FALSE, max_steps = 1000L
    )

    if (is.finite(proposal$total_bic) && proposal$total_bic < current_BIC) {
      current_order <- proposal$topo
      current_IC <- proposal$IC_mat
      current_BIC <- proposal$total_bic
      current_adj <- proposal$adjacency[[1L]]
      current_selected <- proposal$selected[[1L]]
      accepted[t] <- TRUE
    }

    sampled_ordering[t, ] <- current_order
    BIC_score[t] <- current_BIC
    result_IC[[t]] <- current_IC
    Adjacency_matrix[[t]] <- list(current_adj)  # Original nesting retained

    if (t %% 100L == 0L || t == N_target) {
      cat(sprintf("Dataset %02d | iteration %d/%d | BIC = %.6f\n",
                  dataset_id, t, N_target, current_BIC))
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # ----------------------------------------------------------
  # Save the full results in ONE .RData file per dataset
  # ----------------------------------------------------------

  ZIPBN_DAG <- list(
    dataset_id = dataset_id, dataset_info = current_info,
    n_sample = n_sample, removed_index = remove_original_index,
    K = K, p = p, c0 = c0, distribution = model_distribution,
    allocated_cpus = n_cpus, node_workers = n_workers,
    completed = N_target, target_iterations = N_target,
    Adjacency_matrix = Adjacency_matrix, result_IC = result_IC,
    BIC_score = BIC_score, score = BIC_score,
    sampled_ordering = sampled_ordering, accepted = accepted,
    final_ordering = current_order, final_IC = current_IC,
    final_BIC = current_BIC, E_estimate = list(current_adj),
    timeit = as.difftime(elapsed, units = "secs"),
    timeit_seconds = elapsed
  )

  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(
    res_dir, sprintf("ZIPBN_DAG_P_%d_K_1_dataset_%02d.RData", p, dataset_id)
  )
  save(ZIPBN_DAG, file = output_file, compress = FALSE)
  cat("Saved:", output_file, "\n")
  invisible(ZIPBN_DAG)
}

main()
