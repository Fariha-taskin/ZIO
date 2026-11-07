# ==============================================================================
# Consensus-network edge counts across voting thresholds (41 genes, 33 subjects)
# ==============================================================================

library(ggplot2)
library(tidyr)

# ------------------------------------------------------------
# 1. Directories and data
# ------------------------------------------------------------

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
save_dir <- file.path(output_dir, "real_data_41")

p <- 41L
cutoff <- 0.9999
control_ids <- 1:16
case_ids <- 17:33
thresholds <- 1:17

load(file.path(output_dir, "ZIO_DAG_P_41_K_33.RData"))
control <- readRDS(file.path(data_dir, "control_list_41.rds"))
genes <- colnames(control[["control_2"]])

stopifnot(
  length(genes) == p,
  !anyNA(genes),
  !anyDuplicated(genes),
  length(ZIO_DAG$E_estimate) == 33L
)

dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Threshold subject-specific adjacency matrices
# ------------------------------------------------------------

extract_dag <- function(A) {
  A <- as.matrix(A)
  if (!identical(dim(A), c(p, p)) || anyNA(A)) {
    stop("Invalid 41 x 41 subject adjacency matrix.")
  }

  # Use matrix names when available; otherwise retain the existing gene order.
  if (!is.null(rownames(A)) && !is.null(colnames(A))) {
    if (!setequal(rownames(A), genes) || !setequal(colnames(A), genes)) {
      stop("Adjacency-matrix gene names do not match the control data.")
    }
    A <- A[genes, genes, drop = FALSE]
  }

  A <- 1L * (A >= cutoff)
  diag(A) <- 0L
  A
}

dags <- lapply(ZIO_DAG$E_estimate, extract_dag)
control_counts <- Reduce(`+`, dags[control_ids])
case_counts <- Reduce(`+`, dags[case_ids])

# ------------------------------------------------------------
# 3. Count edges for each consensus voting threshold
# ------------------------------------------------------------

ktable <- do.call(rbind, lapply(thresholds, function(ell) {
  control_adj <- control_counts > ell
  case_adj <- case_counts > ell

  data.frame(
    Threshold = ell,
    Control = sum(control_adj),
    Case = sum(case_adj),
    Shared = sum(control_adj & case_adj),
    Difference = sum(xor(control_adj, case_adj))
  )
}))

ktable_long <- tidyr::pivot_longer(
  ktable,
  cols = c("Control", "Case", "Shared", "Difference"),
  names_to = "Type",
  values_to = "EdgeCount"
)
ktable_long$Type <- factor(
  ktable_long$Type,
  levels = c("Control", "Case", "Shared", "Difference")
)

# ------------------------------------------------------------
# 4. Plot and save
# ------------------------------------------------------------

edge_plot <- ggplot(ktable_long, aes(x = Threshold, y = EdgeCount, color = Type)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = c(
    Control = "#1b9e77",
    Case = "#d95f02",
    Shared = "#e7298a",
    Difference = "#7570b3"
  )) +
  scale_x_continuous(
    breaks = c(1, 8, 16),
    labels = c("1", "8", "16"),
    expand = expansion(mult = c(0.04, 0.02))
  ) +
  labs(x = "Threshold(ℓ)", y = "Number of Edges", color = NULL) +
  theme_minimal(base_size = 30) +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = grid::unit(6, "pt"),
    plot.margin = margin(t = 6, r = 50, b = 6, l = 6),
    legend.position = c(0.50, 0.65),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 30),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

print(edge_plot)

# Only the PDF is saved (no extra RData, CSV, or text outputs).
ggsave(
  filename = file.path(save_dir, "Edge_number.pdf"),
  plot = edge_plot,
  device = grDevices::cairo_pdf,
  width = 8,
  height = 7,
  dpi = 300
)
