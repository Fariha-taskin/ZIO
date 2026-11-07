
library(readr)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(viridis)
library(grid)

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


p <- 41L
gene_names <- data.frame(gene = genes)

# The original topological layout assumes an acyclic consensus network.
# Independently fitted subject DAGs need not share a common ordering.
safe_topo_sort <- function(graph, mode = "out") {
  if (!igraph::is_dag(graph)) {
    stop("A consensus/difference network contains a directed cycle. ",
         "A topological layout is unavailable for independently fitted DAGs. ",
         "Use a common-order joint fit or choose a non-topological layout.")
  }
  igraph::topo_sort(graph, mode = mode)
}

control_dag_list <- ZIPBN_DAG$E_estimate[1:16]
case_dag_list    <- ZIPBN_DAG$E_estimate[17:33]

n_control <- length(control_dag_list)
n_case    <- length(case_dag_list)

extract_dags <- function(adj_matrix, cutoff) {
  ifelse(adj_matrix >= cutoff, 1, 0)
}

control_dags <- lapply(
  control_dag_list,
  extract_dags,
  cutoff = 0.9999
)

case_dags <- lapply(
  case_dag_list,
  extract_dags,
  cutoff = 0.9999
)

union_dag_control_subset <- Reduce("+", control_dags)
union_dag_case_subset    <- Reduce("+", case_dags)

union_control_dag <- 1 * (union_dag_control_subset > 8)
union_case_dag    <- 1 * (union_dag_case_subset > 8)

################################################################################
# CONTROL
################################################################################

adj_control <- union_control_dag

rownames(adj_control) <- gene_names[, 1]
colnames(adj_control) <- gene_names[, 1]

g_control <- graph_from_adjacency_matrix(
  adj_control,
  mode = "directed",
  weighted = TRUE,
  diag = FALSE
)

topo_order_control <- as.numeric(
  safe_topo_sort(g_control, mode = "out")
)

V(g_control)$topo_rank <- match(
  V(g_control),
  topo_order_control
)

tg_control <- as_tbl_graph(g_control) %>%
  mutate(
    outdeg = centrality_degree(mode = "out"),
    indeg = centrality_degree(mode = "in"),
    total_deg = indeg + outdeg,
    topo_order = topo_rank
  )

nodes_df_control <- as.data.frame(tg_control) %>%
  select(name, total_deg, topo_order)

zero <- nodes_df_control %>%
  filter(total_deg == 0) %>%
  arrange(topo_order) %>%
  mutate(new_topo = row_number())

nonzero <- nodes_df_control %>%
  filter(total_deg > 0) %>%
  arrange(topo_order) %>%
  mutate(new_topo = row_number() + 4)

nodes_df_control <- bind_rows(zero, nonzero) %>%
  mutate(plot_x = total_deg * 1.40)

edges_df_control <- igraph::as_data_frame(
  g_control,
  what = "edges"
) %>%
  left_join(nodes_df_control, by = c("from" = "name")) %>%
  rename(x = plot_x, y = new_topo) %>%
  left_join(nodes_df_control, by = c("to" = "name")) %>%
  rename(xend = plot_x, yend = new_topo)

nz <- subset(nodes_df_control, total_deg > 0)
z  <- subset(nodes_df_control, total_deg == 0)

p1 <- ggplot(
  nodes_df_control,
  aes(x = plot_x, y = new_topo)
) +
  
  geom_segment(
    data = edges_df_control,
    aes(
      x = x + 0.1 * sign(xend - x),
      y = y + 0.1 * sign(yend - y),
      xend = xend - 0.1 * sign(xend - x),
      yend = yend - 0.1 * sign(yend - y)
    ),
    arrow = arrow(length = unit(4, "mm"), type = "closed"),
    arrow.fill = "gray20",
    color = "gray50",
    alpha = 0.5
  ) +
  
  geom_point(
    aes(size = total_deg, color = total_deg),
    alpha = ifelse(nodes_df_control$total_deg > 0, 0.55, 0.12)
  ) +
  
  geom_text(
    data = nz,
    aes(label = name),
    hjust = 0,
    color = "white",
    nudge_x = 0.45,
    size = 8
  ) +
  
  geom_text(
    data = nz,
    aes(label = name),
    hjust = 0,
    color = "black",
    nudge_x = 0.45,
    size = 7.9
  ) +
  
  geom_text(
    data = z,
    aes(label = name),
    hjust = 0,
    nudge_x = 0.25,
    size = 3,
    color = "grey30",
    alpha = 0.4
  ) +
  
  scale_y_reverse(breaks = vcount(g_control):1) +
  
  scale_size_continuous(range = c(3, 8)) +
  
  scale_color_viridis_c(
    option = "D",
    breaks = c(0, 6, 12, 18, 24),
    labels = c("0", "6", "12", "18", "24")
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.03, 0.10))
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Total Degree (# edges)",
    y = "Topological Ordering",
    size = "Degree",
    color = "Degree"
  ) +
  
  theme_minimal(base_size = 25) +
  
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = NA, color = NA),
    legend.box.just = "right",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    plot.margin = margin(10, 40, 10, 10)
  ) +
  
  guides(size = "none")

p1

ggsave(
  file.path(save_dir, "Real_output_control_l_8.pdf"),
  p1,
  device = cairo_pdf,
  dpi = 300,
  width = 11,
  height = 12,
  units = "in",
  bg = "white"
)

################################################################################
# CASE
################################################################################

adj_case <- union_case_dag

rownames(adj_case) <- gene_names[, 1]
colnames(adj_case) <- gene_names[, 1]

g_case <- graph_from_adjacency_matrix(
  adj_case,
  mode = "directed",
  weighted = TRUE,
  diag = FALSE
)

topo_order_case <- as.numeric(
  safe_topo_sort(g_case, mode = "out")
)

V(g_case)$topo_rank <- match(
  V(g_case),
  topo_order_case
)

tg_case <- as_tbl_graph(g_case) %>%
  mutate(
    outdeg = centrality_degree(mode = "out"),
    indeg = centrality_degree(mode = "in"),
    total_deg = indeg + outdeg,
    topo_order = topo_rank
  )

nodes_df_case <- as.data.frame(tg_case) %>%
  select(name, total_deg, topo_order)

zero <- nodes_df_case %>%
  filter(total_deg == 0) %>%
  arrange(topo_order) %>%
  mutate(new_topo = row_number())

nonzero <- nodes_df_case %>%
  filter(total_deg > 0) %>%
  arrange(topo_order) %>%
  mutate(new_topo = row_number() + 2)

nodes_df_case <- bind_rows(zero, nonzero) %>%
  mutate(plot_x = total_deg * 1.40)

edges_df_case <- igraph::as_data_frame(
  g_case,
  what = "edges"
) %>%
  left_join(nodes_df_case, by = c("from" = "name")) %>%
  rename(x = plot_x, y = new_topo) %>%
  left_join(nodes_df_case, by = c("to" = "name")) %>%
  rename(xend = plot_x, yend = new_topo)

nz <- subset(nodes_df_case, total_deg > 0)
z  <- subset(nodes_df_case, total_deg == 0)

p2 <- ggplot(
  nodes_df_case,
  aes(x = plot_x, y = new_topo)
) +
  
  geom_segment(
    data = edges_df_case,
    aes(
      x = x + 0.1 * sign(xend - x),
      y = y + 0.1 * sign(yend - y),
      xend = xend - 0.1 * sign(xend - x),
      yend = yend - 0.1 * sign(yend - y)
    ),
    arrow = arrow(length = unit(4, "mm"), type = "closed"),
    arrow.fill = "gray20",
    color = "gray50",
    alpha = 0.5
  ) +
  
  geom_point(
    aes(size = total_deg, color = total_deg),
    alpha = ifelse(nodes_df_case$total_deg > 0, 0.55, 0.12)
  ) +
  
  geom_text(
    data = nz,
    aes(label = name),
    hjust = 0,
    color = "white",
    nudge_x = 0.45,
    size = 8
  ) +
  
  geom_text(
    data = nz,
    aes(label = name),
    hjust = 0,
    color = "black",
    nudge_x = 0.45,
    size = 7.9
  ) +
  
  geom_text(
    data = z,
    aes(label = name),
    hjust = 0,
    nudge_x = 0.25,
    size = 3,
    color = "grey30",
    alpha = 0.4
  ) +
  
  scale_y_reverse(breaks = vcount(g_case):1) +
  
  scale_size_continuous(range = c(3, 8)) +
  
  scale_color_viridis_c(
    option = "D",
    breaks = c(0, 6, 12, 18, 22),
    labels = c("0", "6", "12", "18", "24")
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.03, 0.10))
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Total Degree (# edges)",
    y = "Topological Ordering",
    size = "Degree",
    color = "Degree"
  ) +
  
  theme_minimal(base_size = 25) +
  
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = NA, color = NA),
    legend.box.just = "right",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    plot.margin = margin(10, 40, 10, 10)
  ) +
  
  guides(size = "none")

p2

ggsave(
  file.path(save_dir, "Real_output_case_l_8.pdf"),
  p2,
  device = cairo_pdf,
  dpi = 300,
  width = 11,
  height = 12,
  units = "in",
  bg = "white"
)

################################################################################
# INTERSECTION
################################################################################

g_intersect <- igraph::intersection(
  g_control,
  g_case
)

topo_order_intersect <- as.numeric(
  safe_topo_sort(g_intersect, mode = "out")
)

V(g_intersect)$topo_rank <- match(
  V(g_intersect),
  topo_order_intersect
)

tg_intersect <- as_tbl_graph(g_intersect) %>%
  mutate(
    outdeg = centrality_degree(mode = "out"),
    indeg = centrality_degree(mode = "in"),
    total_deg = indeg + outdeg,
    topo_order = topo_rank
  )

nodes_df_intersect <- as.data.frame(tg_intersect) %>%
  select(name, total_deg, topo_order)

zero <- nodes_df_intersect %>%
  filter(total_deg == 0) %>%
  arrange(topo_order) %>%
  mutate(new_topo = row_number())

nonzero <- nodes_df_intersect %>%
  filter(total_deg > 0) %>%
  arrange(topo_order) %>%
  mutate(new_topo = row_number() + 4)

nodes_df_intersect <- bind_rows(zero, nonzero) %>%
  mutate(plot_x = total_deg * 1.40)

edges_df_intersect <- igraph::as_data_frame(
  g_intersect,
  what = "edges"
) %>%
  left_join(nodes_df_intersect, by = c("from" = "name")) %>%
  rename(x = plot_x, y = new_topo) %>%
  left_join(nodes_df_intersect, by = c("to" = "name")) %>%
  rename(xend = plot_x, yend = new_topo)

nz <- subset(nodes_df_intersect, total_deg > 0)
z  <- subset(nodes_df_intersect, total_deg == 0)

p3 <- ggplot(
  nodes_df_intersect,
  aes(x = plot_x, y = new_topo)
) +
  
  geom_segment(
    data = edges_df_intersect,
    aes(
      x = x + 0.1 * sign(xend - x),
      y = y + 0.1 * sign(yend - y),
      xend = xend - 0.1 * sign(xend - x),
      yend = yend - 0.1 * sign(yend - y)
    ),
    arrow = arrow(length = unit(4, "mm"), type = "closed"),
    arrow.fill = "gray20",
    color = "gray50",
    alpha = 0.5
  ) +
  
  geom_point(
    aes(size = total_deg, color = total_deg),
    alpha = ifelse(nodes_df_intersect$total_deg > 0, 0.55, 0.12)
  ) +
  
  geom_text(
    data = nz,
    aes(label = name),
    hjust = 0,
    color = "white",
    nudge_x = 0.45,
    size = 8
  ) +
  
  geom_text(
    data = nz,
    aes(label = name),
    hjust = 0,
    color = "black",
    nudge_x = 0.45,
    size = 7.9
  ) +
  
  geom_text(
    data = z,
    aes(label = name),
    hjust = 0,
    nudge_x = 0.25,
    size = 3,
    color = "grey30",
    alpha = 0.4
  ) +
  
  scale_y_reverse(breaks = vcount(g_intersect):1) +
  
  scale_size_continuous(range = c(3, 8)) +
  
  scale_color_viridis_c(
    option = "D",
    breaks = c(0, 6, 12, 18),
    labels = c("0", "6", "12", "18")
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.03, 0.10))
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Total Degree (# edges)",
    y = "Topological Ordering",
    size = "Degree",
    color = "Degree"
  ) +
  
  theme_minimal(base_size = 25) +
  
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = NA, color = NA),
    legend.box.just = "right",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    plot.margin = margin(10, 40, 10, 10)
  ) +
  
  guides(size = "none")

p3

ggsave(
  file.path(save_dir, "Real_output_intersection_l_8.pdf"),
  p3,
  device = cairo_pdf,
  dpi = 300,
  width = 11,
  height = 12,
  units = "in",
  bg = "white"
)

################################################################################
# DIFFERENCE
################################################################################

g_intersection <- igraph::intersection(
  g_control,
  g_case
)

g_control_only <- igraph::difference(
  g_control,
  g_intersection
)

g_case_only <- igraph::difference(
  g_case,
  g_intersection
)

g_union <- igraph::union(
  g_control,
  g_case
)

all_nodes <- V(g_union)$name

edges_control_only <- igraph::as_data_frame(
  g_control_only,
  what = "edges"
) %>%
  mutate(type = "Control Only")

edges_case_only <- igraph::as_data_frame(
  g_case_only,
  what = "edges"
) %>%
  mutate(type = "Case Only")

edges_combined <- bind_rows(
  edges_control_only,
  edges_case_only
)

g_combined <- graph_from_data_frame(
  edges_combined,
  directed = TRUE,
  vertices = data.frame(name = all_nodes)
)

topo_order <- as.numeric(
  safe_topo_sort(g_combined, mode = "out")
)

V(g_combined)$topo_rank <- match(
  V(g_combined),
  topo_order
)

tg <- as_tbl_graph(g_combined) %>%
  mutate(
    outdeg = centrality_degree(mode = "out"),
    indeg = centrality_degree(mode = "in"),
    total_deg = indeg + outdeg,
    topo_order = topo_rank
  )

nodes_df <- as.data.frame(tg) %>%
  select(name, total_deg, topo_order)

zero <- nodes_df %>%
  filter(total_deg == 0) %>%
  arrange(topo_order) %>%
  mutate(new_topo = row_number())

nonzero <- nodes_df %>%
  filter(total_deg > 0) %>%
  arrange(topo_order) %>%
  mutate(new_topo = row_number() + 3)

nodes_df <- bind_rows(zero, nonzero) %>%
  mutate(plot_x = total_deg * 1.40)

edges_df <- edges_combined %>%
  left_join(nodes_df, by = c("from" = "name")) %>%
  rename(x = plot_x, y = new_topo) %>%
  left_join(nodes_df, by = c("to" = "name")) %>%
  rename(xend = plot_x, yend = new_topo)

nz <- subset(nodes_df, total_deg > 0)
z  <- subset(nodes_df, total_deg == 0)

p4 <- ggplot(
  nodes_df,
  aes(x = plot_x, y = new_topo)
) +
  
  geom_segment(
    data = edges_df,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = type
    ),
    arrow = arrow(length = unit(4, "mm"), type = "closed"),
    alpha = 0.6,
    linewidth = 0.6
  ) +
  
  geom_point(
    aes(size = total_deg, fill = total_deg),
    shape = 21,
    color = "black",
    stroke = 0.2,
    alpha = ifelse(nodes_df$total_deg > 0, 0.55, 0.12)
  ) +
  
  geom_text(
    data = nz,
    aes(label = name),
    hjust = 0,
    color = "white",
    nudge_x = 0.45,
    size = 7.5
  ) +
  
  geom_text(
    data = nz,
    aes(label = name),
    hjust = 0,
    color = "black",
    nudge_x = 0.45,
    size = 7.4
  ) +
  
  geom_text(
    data = z,
    aes(label = name),
    hjust = 0,
    nudge_x = 0.25,
    size = 3,
    color = "grey30",
    alpha = 0.4
  ) +
  
  scale_y_reverse(breaks = vcount(g_combined):1) +
  
  scale_fill_viridis_c(option = "D") +
  
  scale_color_manual(
    values = c(
      "Control Only" = "#E41A1C",
      "Case Only" = "#377EB8"
    )
  ) +
  
  scale_size_continuous(range = c(2.5, 7)) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.03, 0.10))
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Total Degree (# edges)",
    y = "Topological Ordering",
    color = "Edge Source",
    fill = "Degree",
    size = "Degree"
  ) +
  
  theme_minimal(base_size = 25) +
  
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = NA, color = NA),
    legend.box.just = "right",
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    plot.margin = margin(10, 40, 10, 10)
  ) +
  
  guides(
    size = "none",
    fill = "none"
  )

p4

ggsave(
  file.path(save_dir, "Real_output_difference_l_8.pdf"),
  p4,
  device = cairo_pdf,
  dpi = 300,
  width = 11,
  height = 12,
  units = "in",
  bg = "white"
)
