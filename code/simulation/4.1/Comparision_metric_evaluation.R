# ==============================================================================
# Simulation 4.1: ODS, MRS and ZIPBN performance comparison
# ==============================================================================

library(dplyr)

project_dir <- "ZIO"
if (!dir.exists(file.path(project_dir, "code", "helper"))) {
  if (dir.exists(file.path("code", "helper"))) {
    project_dir <- "."
  } else {
    stop("Cannot locate ZIO/code/helper. Run from ZIO/ or its parent directory.")
  }
}
output_dir <- file.path(project_dir, "output", "simulation", "4.1")

# All three benchmark scripts currently use p = 50. change to 25, and 75
# Change p consistently in the baseline scripts and here for other settings.
p <- 50L
n_sim <- 30L
zipbn_burnin <- 1500L
zipbn_cutoff <- 0.5  # Posterior inclusion-probability threshold; adjust if needed.

# Load each method into a separate environment to avoid overwriting data0.
load_method <- function(method, p) {
  path <- file.path(output_dir, sprintf("%s_z%d.RData", method, p))
  if (!file.exists(path)) {
    stop("Missing ", method, " results for p = ", p, ": ", path,
         "\nAll three methods must be run with the same p.")
  }
  env <- new.env(parent = baseenv())
  load(path, envir = env)
  env
}

ods_env <- load_method("ods", p)
mrs_env <- load_method("mrs", p)
zip_env <- load_method("zipbn", p)

ods <- ods_env$ods
mrs <- mrs_env$mrs
zipbn <- zip_env$zipbn

stopifnot(length(ods) >= n_sim,
          length(mrs) >= n_sim,
          length(zipbn) >= n_sim,
          length(ods_env$data0) >= n_sim,
          length(mrs_env$data0) >= n_sim,
          length(zip_env$data) >= n_sim)

# Compute all metrics from binary adjacency matrices.
get_metrics <- function(truth, estimated) {
  truth <- as.matrix(truth)
  estimated <- as.matrix(estimated)
  if (!identical(dim(truth), c(p, p)) ||
      !identical(dim(estimated), c(p, p)) ||
      anyNA(truth) || anyNA(estimated) ||
      !all(truth %in% 0:1) || !all(estimated %in% 0:1)) {
    stop("Both graphs must be binary p x p matrices without missing values.")
  }

  TP <- sum(estimated == 1 & truth == 1)
  TN <- sum(estimated == 0 & truth == 0)
  FP <- sum(estimated == 1 & truth == 0)
  FN <- sum(estimated == 0 & truth == 1)
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

  c(TPR = if (TP + FN) TP / (TP + FN) else NA_real_,
    FDR = if (TP + FP) FP / (TP + FP) else NA_real_,
    MCC = if (denom > 0) (TP * TN - FP * FN) / denom else NA_real_,
    HD = FP + FN)
}

get_zipbn_dag <- function(samples) {
  edge_samples <- samples$E
  if (length(dim(edge_samples)) != 3L ||
      !identical(dim(edge_samples)[1:2], c(p, p)) ||
      dim(edge_samples)[3] <= zipbn_burnin) {
    stop("Unexpected ZIPBN posterior edge-array dimensions.")
  }
  kept <- seq.int(zipbn_burnin + 1L, dim(edge_samples)[3])
  1L * (apply(edge_samples[, , kept, drop = FALSE], c(1, 2), mean) >= zipbn_cutoff)
}

# Confirm that all benchmark methods generated the same true DAG per replicate.
compare_rep <- function(i) {
  truth <- as.matrix(ods_env$data0[[i]]$E_true[[1]])
  truth_mrs <- as.matrix(mrs_env$data0[[i]]$E_true[[1]])
  truth_zipbn <- as.matrix(zip_env$data[[i]]$E_true[[1]])
  if (!isTRUE(all.equal(truth, truth_mrs)) ||
      !isTRUE(all.equal(truth, truth_zipbn))) {
    stop("True DAGs differ among methods for replicate ", i, ".")
  }
  rbind(
    data.frame(METHOD = "ODS", Replicate = i,
               as.list(get_metrics(truth, ods[[i]]))),
    data.frame(METHOD = "MRS", Replicate = i,
               as.list(get_metrics(truth, mrs[[i]]))),
    data.frame(METHOD = "ZIPBN", Replicate = i,
               as.list(get_metrics(truth, get_zipbn_dag(zipbn[[i]]))))
  )
}

OC <- bind_rows(lapply(seq_len(n_sim), compare_rep))
OC$METHOD <- factor(OC$METHOD, levels = c("MRS", "ODS", "ZIPBN"))

summary_OC <- OC %>%
  group_by(METHOD) %>%
  summarise(
    across(c(TPR, FDR, MCC, HD),
           list(mean = ~mean(.x, na.rm = TRUE),
                se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))),
           .names = "{.fn}_{.col}"),
    .groups = "drop"
  )

print(summary_OC)
