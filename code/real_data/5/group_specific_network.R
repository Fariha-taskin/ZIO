################################################################################
# Shared and group-enriched edge-frequency summaries: Control vs MDD
################################################################################

# ==============================================================================
# 1. Settings and data
# ==============================================================================

# Each independent result stores its final DAG in E_estimate.
top_n <- 5
project_dir <- "ZIO"
if (!dir.exists(file.path(project_dir, "code", "helper"))) {
  if (dir.exists(file.path("code", "helper"))) {
    project_dir <- "."
  } else {
    stop("Cannot locate ZIO/code/helper. Run from ZIO/ or its parent directory.")
  }
}

data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "output")

load(file.path(output_dir, "ZIPBN_DAG_P_41_K_33.RData"))
control <- readRDS(file.path(data_dir, "control_list_41.rds"))
genes <- colnames(control[["control_2"]])
stopifnot(length(genes) == 41L,length(ZIPBN_DAG$E_estimate) == 33L)

save_dir <- file.path(output_dir, "real_data_41")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
dags <- ZIPBN_DAG$E_estimate

# Assumption: adjacency[i, j] = 1 represents gene i -> gene j.
dags <- lapply(seq_along(dags), function(k) {
  G <- as.matrix(dags[[k]])
  
  if (!identical(dim(G), c(41L, 41L)) ||
      is.null(rownames(G)) || is.null(colnames(G)) ||
      anyDuplicated(rownames(G)) || anyDuplicated(colnames(G)) ||
      !setequal(rownames(G), genes) ||
      !setequal(colnames(G), genes)) {
    stop("Invalid dimensions or gene names in subject DAG ", k)
  }
  
  G <- G[genes, genes, drop = FALSE]
  
  if (anyNA(G) || !all(G %in% c(0, 1)) ||
      any(diag(G) != 0)) {
    stop("Subject DAG ", k, " must have binary entries and zero diagonal.")
  }
  
  G
})

# These indices assume that subjects are ordered as 16 controls, then 17 MDD.
control_ids <- 1:16
mdd_ids <- 17:33
K_control <- length(control_ids)
K_mdd <- length(mdd_ids)

# ==============================================================================
# 2. Cross-subject edge-selection frequencies
# ==============================================================================

control_count <- Reduce("+", dags[control_ids])
mdd_count <- Reduce("+", dags[mdd_ids])

control_freq <- control_count / K_control
mdd_freq <- mdd_count / K_mdd

edges <- expand.grid(
  from = seq_along(genes),
  to = seq_along(genes),
  KEEP.OUT.ATTRS = FALSE
)
edges <- edges[edges$from != edges$to, , drop = FALSE]

ij <- as.matrix(edges[, c("from", "to")])
edges$n_control <- control_count[ij]
edges$n_mdd <- mdd_count[ij]
edges$pi_control <- control_freq[ij]
edges$pi_mdd <- mdd_freq[ij]

edges$shared_score <- pmin(edges$pi_control, edges$pi_mdd)
edges$control_score <- pmax(edges$pi_control - edges$pi_mdd, 0)
edges$mdd_score <- pmax(edges$pi_mdd - edges$pi_control, 0)

# ==============================================================================
# 3. Top edges within each component
# ==============================================================================

get_top_edges <- function(category, score, n_top = top_n) {
  idx <- which(score > 0)
  idx <- idx[order(-score[idx], edges$from[idx], edges$to[idx])]
  idx <- head(idx, n_top)
  
  result <- data.frame(
    Category = rep(category, length(idx)),
    Rank = seq_along(idx),
    Directed_edge = sprintf(
      "%s -> %s", genes[edges$from[idx]], genes[edges$to[idx]]
    ),
    pi_control = sprintf(
      "%.3f (%d/%d)",
      edges$pi_control[idx], edges$n_control[idx], K_control
    ),
    pi_mdd = sprintf(
      "%.3f (%d/%d)",
      edges$pi_mdd[idx], edges$n_mdd[idx], K_mdd
    ),
    Score = sprintf("%.3f", score[idx])
  )
  names(result) <- c(
    "Category", "Rank", "Directed edge", "pi(ctrl)", "pi(MDD)", "Score"
  )
  result
}

edge_table <- rbind(
  get_top_edges("Shared", edges$shared_score),
  get_top_edges("Control-enriched", edges$control_score),
  get_top_edges("MDD-enriched", edges$mdd_score)
)

cat("\nTop edges in each component:\n")
print(edge_table, row.names = FALSE)

# ==============================================================================
# 4. All edges present in every control and every MDD subject
# ==============================================================================

shared_in_all <- (
  edges$n_control == K_control &
    edges$n_mdd == K_mdd
)

shared_all_table <- get_top_edges(
  category = "Shared in all subjects",
  score = ifelse(shared_in_all, edges$shared_score, 0),
  n_top = sum(shared_in_all)
)

cat("\nEdges present in all 33 subject-specific DAGs:\n")
print(shared_all_table, row.names = FALSE)

# ==============================================================================
# 5. Overall shared-network summaries
# ==============================================================================

union_weight <- sum(pmax(control_freq, mdd_freq))

weighted_jaccard <- if (union_weight > 0) {
  sum(pmin(control_freq, mdd_freq)) / union_weight
} else {
  NA_real_
}

shared_in_all_33 <- sum(shared_in_all)
shared_in_both_groups <- sum(edges$shared_score > 0)

cat("\nWeighted Jaccard:", round(weighted_jaccard, 3), "\n")
cat("Edges present in all 33 DAGs:", shared_in_all_33, "\n")
cat("Edges selected in at least one subject in each group:",
    shared_in_both_groups, "\n")

# Save the printed summaries in the same output directory as the figures.
summary_lines <- c(
  "Top edges in each component:",
  capture.output(print(edge_table, row.names = FALSE)),
  "",
  "Edges present in all 33 subject-specific DAGs:",
  capture.output(print(shared_all_table, row.names = FALSE)),
  "",
  paste0("Weighted Jaccard: ", round(weighted_jaccard, 3)),
  paste0("Edges present in all 33 DAGs: ", shared_in_all_33),
  paste0("Edges selected in both groups: ", shared_in_both_groups)
)
writeLines(
  summary_lines,
  file.path(save_dir, "group_specific_network_summary_P41.txt")
)
