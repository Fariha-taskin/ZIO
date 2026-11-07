library(igraph)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbump)

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

ell <- 8
cutoff <- 0.9999
p <- 41L

# Construct consensus networks and calculate degrees.
get_degrees <- function(ids, group) {
  dags <- lapply(
    ZIO_DAG$E_estimate[ids],
    function(A) 1L * (A >= cutoff)
  )
  A <- 1L * (Reduce(`+`, dags) > ell)
  dimnames(A) <- list(genes, genes)
  
  g <- graph_from_adjacency_matrix(
    A, mode = "directed", diag = FALSE
  )
  
  result <- data.frame(
    gene = V(g)$name,
    in_deg = degree(g, mode = "in"),
    out_deg = degree(g, mode = "out"),
    total_deg = degree(g, mode = "all")
  )
  names(result)[-1] <- paste0(group, "_", names(result)[-1])
  result
}

# Exclude genes isolated in both groups, then calculate dense ranks.
ranked <- inner_join(
  get_degrees(1:16, "control"),
  get_degrees(17:33, "case"),
  by = "gene"
) %>%
  filter(if_any(-gene, ~ .x > 0)) %>%
  mutate(across(-gene, ~ dense_rank(desc(.x)))) %>%
  pivot_longer(
    -gene,
    names_to = c("condition", "metric"),
    names_pattern = "(control|case)_(in|out|total)_deg",
    values_to = "rank"
  ) %>%
  group_by(metric, condition, rank) %>%
  mutate(rank_adj = rank + (row_number() - 1) * 0.05) %>%
  ungroup()

plot_ranks <- function(metric_name, label_spread = 0.45) {
  df <- ranked %>%
    filter(metric == metric_name) %>%
    mutate(x_num = match(condition, c("control", "case")))
  
  labels <- df %>%
    filter(condition == "case") %>%
    arrange(rank, rank_adj, gene) %>%
    group_by(rank) %>%
    mutate(
      x_step = 2.015,
      x_label = 2.10 + if_else(rank %% 2 == 0, 0.70, 0),
      y_label = rank_adj +
        (row_number() - (n() + 1) / 2) * label_spread
    ) %>%
    ungroup()
  
  # Prevent overlaps within each label column.
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
  
  ggplot(df, aes(x_num, rank_adj, group = gene, colour = gene)) +
    geom_bump(linewidth = 1, smooth = 6, alpha = 0.9) +
    geom_point(data = labels, size = 1.6) +
    geom_segment(
      data = labels, inherit.aes = FALSE,
      aes(
        x = 2, y = rank_adj, xend = x_step,
        yend = rank_adj, colour = gene
      ),
      linewidth = 0.8, alpha = 0.9
    ) +
    geom_segment(
      data = labels, inherit.aes = FALSE,
      aes(
        x = x_step, y = rank_adj,
        xend = x_label - 0.03, yend = y_label, colour = gene
      ),
      linewidth = 0.8, alpha = 0.95,
      arrow = grid::arrow(
        length = grid::unit(0.08, "in"), type = "closed"
      )
    ) +
    # geom_text(
    #   data = labels, inherit.aes = FALSE,
    #   aes(x_label, y_label, label = gene, colour = gene),
    #   hjust = 0, size = 12
    # ) +
    shadowtext::geom_shadowtext(
      data = labels,
      inherit.aes = FALSE,
      aes(x_label, y_label, label = gene, colour = gene),
      hjust = 0,
      size = 12,
      bg.colour = "white",
      bg.r = 0.15
    ) +
    scale_x_continuous(
      breaks = 1:2,
      labels = c("Control", "Case"),
      limits = c(1, max(labels$x_label) + 0.7),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_reverse(
      breaks = seq_len(floor(max(df$rank_adj)))
    ) +
    scale_colour_viridis_d(option = "turbo") +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 40) +
    theme(
      legend.position = "none",
      axis.ticks = element_line(colour = "black", linewidth = 0.8),
      axis.ticks.length = grid::unit(6, "pt"),
      axis.text = element_text(size = 40),
      panel.grid.minor = element_blank(),
      plot.margin = margin(6, 20, 6, 20),
      panel.border = element_rect(
        colour = "black", fill = NA, linewidth = 0.6
      )
    )
}

p_total <- plot_ranks("total")
p_out <- plot_ranks("out")
p_in <- plot_ranks("in")

p_total
p_out

ggsave(
  file.path(save_dir, "rank_plot_8_total.pdf"),
  plot = p_total,
  device = cairo_pdf,
  width = 12,
  height = 15
)

ggsave(
  file.path(save_dir, "rank_plot_8_out.pdf"),
  plot = p_out,
  device = cairo_pdf,
  width = 12,
  height = 15
)


ggsave(file.path(save_dir, "rank_plot_8_out.png"), p_out,
       width = 10, height = 14, dpi = 200, bg = "white")

ggsave(file.path(save_dir, "rank_plot_8_total.png"), p_total,
       width = 10, height = 14, dpi = 200, bg = "white")



ranked %>%
  filter(metric %in% c("out", "total")) %>%
  select(gene, metric, condition, rank) %>%
  pivot_wider(names_from = condition, values_from = rank) %>%
  arrange(metric, control, case) %>%
  print(n = Inf)
