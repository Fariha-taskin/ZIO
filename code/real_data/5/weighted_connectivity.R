library(dplyr)
library(ggplot2)
library(tidyr)
library(ggbump)
library(shadowtext)

# ==============================================================================
# 1. Load data
# ==============================================================================

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

load(file.path(output_dir, "ZIO_DAG_P_41_K_33.RData"))
control <- readRDS(file.path(data_dir, "control_list_41.rds"))
genes <- colnames(control[["control_2"]])
stopifnot(length(genes) == 41L,length(ZIO_DAG$E_estimate) == 33L)

save_dir <- file.path(output_dir, "real_data_41")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

p <- length(genes)

control_ids <- 1:16
mdd_ids <- 17:33
cutoff <- 0.9999

# ==============================================================================
# 2. Construct subject-specific graphs and group edge frequencies
# ==============================================================================

extract_dag <- function(A) {
  A <- as.matrix(A)
  
  if (!all(dim(A) == c(p, p)) || anyNA(A)) {
    stop("Check adjacency-matrix dimensions and missing values.")
  }
  
  # Align named matrices with the gene vector.
  if (!is.null(rownames(A)) && !is.null(colnames(A))) {
    if (!setequal(rownames(A), genes) ||
        !setequal(colnames(A), genes)) {
      stop("Adjacency-matrix gene names do not match 'genes'.")
    }
    
    A <- A[genes, genes, drop = FALSE]
  }
  
  A <- 1L * (A >= cutoff)
  diag(A) <- 0L
  dimnames(A) <- list(genes, genes)
  
  A
}

control_dags <- lapply(
  ZIO_DAG$E_estimate[control_ids], extract_dag
)

mdd_dags <- lapply(
  ZIO_DAG$E_estimate[mdd_ids], extract_dag
)

# control_freq <- Reduce(`+`, control_dags) / length(control_dags)
# mdd_freq <- Reduce(`+`, mdd_dags) / length(mdd_dags)

union_dag_control_subset <- matrix(0, p, p)
for (i in 1:p){
  for (j in 1:p){
    for (k in 1:length(control_dags)){
      union_dag_control_subset[i, j] <- union_dag_control_subset[i, j] + control_dags[[k]][i,j]
    }
  }
}

control_freq <- union_dag_control_subset / length(control_dags)

union_dag_case_subset <- matrix(0, p, p)
for (i in 1:p){
  for (j in 1:p){
    for (k in 1:length(mdd_dags)){
      union_dag_case_subset[i, j] <- union_dag_case_subset[i, j] + mdd_dags[[k]][i,j]
    }
  }
}
mdd_freq <- union_dag_case_subset / length(mdd_dags)

# ==============================================================================
# 3. Weighted total connectivity: outgoing + incoming frequencies
# ==============================================================================

connectivity <- data.frame(
  Gene = genes,
  Control = rowSums(control_freq) + colSums(control_freq),
  MDD = rowSums(mdd_freq) + colSums(mdd_freq)
) %>%
  mutate(
    Difference = Control - MDD,
    Absolute_difference = abs(Difference)
  )

top_genes <- connectivity %>%
  slice_max(Absolute_difference, n = 20, with_ties = FALSE) %>%
  arrange(desc(Absolute_difference), Gene)

top_genes$Gene <- factor(
  top_genes$Gene,
  levels = rev(top_genes$Gene)
)




# ==============================================================================
# 4. Rank all 41 genes, then retain the selected top 20 for plotting
# ==============================================================================

ranked <- connectivity %>%
  mutate(
    # Round only to prevent numerical precision from breaking ties.
    Control_rank = dense_rank(desc(round(Control, 10))),
    MDD_rank = dense_rank(desc(round(MDD, 10)))
  ) %>%
  filter(Gene %in% as.character(top_genes$Gene)) %>%
  select(Gene, Control_rank, MDD_rank) %>%
  pivot_longer(
    cols = c(Control_rank, MDD_rank),
    names_to = "condition",
    names_pattern = "(Control|MDD)_rank",
    values_to = "rank"
  ) %>%
  mutate(
    gene = Gene,
    x_num = match(condition, c("Control", "MDD"))
  )

# Keep exact rank positions, including ties.
labels <- ranked %>%
  filter(condition == "MDD") %>%
  arrange(rank, gene) %>%
  group_by(rank) %>%
  mutate(
    x_step = 2.015,
    x_label = 2.10 + if_else(rank %% 2 == 0, 0.70, 0),
    y_label = rank +
      (row_number() - (n() + 1) / 2) * 0.45
  ) %>%
  ungroup()

# Separate labels within each label column.
labels <- labels %>%
  group_by(x_label) %>%
  arrange(y_label, gene, .by_group = TRUE) %>%
  mutate(
    original_y = y_label,
    offset = (row_number() - 1) * 0.65,
    y_label = cummax(original_y - offset) + offset,
    y_label = y_label - mean(y_label - original_y)
  ) %>%
  ungroup()

# ==============================================================================
# 5. Plot weighted total connectivity
# ==============================================================================

plot_data <- connectivity %>%
  slice_max(Absolute_difference, n = 10, with_ties = FALSE) %>%
  arrange(desc(Absolute_difference), Gene) %>%
  mutate(
    Gene = factor(Gene, levels = rev(Gene)),
    midpoint = (Control + MDD) / 2,
    difference_label = sprintf("%.2f", Absolute_difference)
  )

point_data <- plot_data %>%
  pivot_longer(
    cols = c(Control, MDD),
    names_to = "Condition",
    values_to = "Connectivity"
  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("Control", "MDD"))
  )

connectivity_plot <- ggplot(plot_data, aes(y = Gene)) +
  # Connect the Control and MDD values.
  geom_segment(
    aes(x = Control, xend = MDD, yend = Gene),
    color = "grey65",
    linewidth = 0.45
  ) +
  geom_point(
    data = point_data,
    aes(x = Connectivity, color = Condition),
    size = 2.4
  ) +
  # Absolute difference above the midpoint of each pair.
  geom_text(
    aes(x = midpoint, label = difference_label),
    nudge_y = 0.24,
    size = 3.2,
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      Control = "#56B4E9",
      MDD = "#B07D52"
    )
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 5),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  scale_y_discrete(
    expand = expansion(add = c(0.55, 0.65))
  ) +
  labs(
    x = "Weighted total connectivity",
    y = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.35
    ),
    panel.grid.minor.x = element_line(
      color = "grey94",
      linewidth = 0.25
    ),
    axis.text = element_text(color = "black"),
    axis.ticks.y = element_blank(),
    legend.position = c(0.87, 0.12),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.height = grid::unit(0.35, "cm"),
    plot.margin = margin(8, 10, 5, 5)
  )

print(connectivity_plot)

ggsave(
  filename = file.path(save_dir, "weighted_total_connectivity_top10.pdf"),
  plot = connectivity_plot,
  width = 4.3,
  height = 3.1,
  units = "in"
)