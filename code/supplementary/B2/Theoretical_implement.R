############################################################
## 0. Load existing code
############################################################
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
source(file.path(helper_dir, "data_generator_for_given_order.R"))
source(file.path(helper_dir, "metric.R"))
source(file.path(helper_dir, "Theorem_numeric_sampled_order_helper.R"))

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

############################################################
## 1. Simulation settings
############################################################
set.seed(123)

p_grid <- c(5, 10)
n_grid <- c(250, 500, 1000, 2000)

K <- 1
n_rep <- 30

## For p = 10, we sample this many valid and invalid orderings.
## For p = 5, the code below will override these and use all possible orderings.
n_valid_sample_p10 <- 30
n_invalid_sample_p10 <- 30

c0 <- 1
criterion <- "BIC"
dist <- "poisson"

############################################################
## 2. Output folders and progress log
############################################################
out_dir <- file.path(project_dir, "output", "simulation", "theorem_numeric", "sampled_orders")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

progress_log <- file.path(out_dir, "progress_log.txt")

cat(
  "Starting theorem numerical illustration\n",
  "Start time: ", as.character(Sys.time()), "\n\n",
  file = progress_log,
  append = FALSE
)

############################################################
## 3. Run theoretical numerical illustration
############################################################
res_all <- list()
counter <- 1L

for (p in p_grid) {
  
  cat("\n====================================\n", flush = TRUE)
  cat("Running p =", p, "\n", flush = TRUE)
  cat("====================================\n", flush = TRUE)
  
  cat(
    "\n====================================\n",
    "Running p = ", p, "\n",
    "====================================\n",
    file = progress_log,
    append = TRUE
  )
  
  true_topo <- 1:p
  n_edge <- p
  
  ## Fix one true DAG for this p
  set.seed(1000 + p)
  
  E_true <- make_E_from_order(true_topo, n_edge)
  stopifnot(is_E_consistent_with_order(E_true, true_topo))
  
  alpha_true <- matrix(0, p, p)
  beta_true  <- matrix(0, p, p)
  
  edge_count <- sum(E_true == 1L)
  
  if (edge_count > 0) {
    alpha_true[E_true == 1L] <- runif(edge_count, 0.3, 1.0)
    beta_true[E_true == 1L]  <- runif(edge_count, -1.0, -0.3)
  }
  
  delta_true <- runif(p, -2.0, -1.0)
  gamma_true <- runif(p,  1.0,  2.0)
  
  ## All valid topological orderings
  L_G_list <- all_topological_orders(E_true)
  L_G_mat  <- do.call(rbind, L_G_list)
  
  n_valid_all <- nrow(L_G_mat)
  n_total_all <- factorial(p)
  n_invalid_all <- n_total_all - n_valid_all
  
  cat("Number of valid orderings for p =", p, ":", n_valid_all, "\n", flush = TRUE)
  cat("Number of invalid orderings for p =", p, ":", n_invalid_all, "\n", flush = TRUE)
  
  cat(
    "Number of valid orderings for p = ", p, ": ", n_valid_all, "\n",
    "Number of invalid orderings for p = ", p, ": ", n_invalid_all, "\n",
    file = progress_log,
    append = TRUE
  )
  
  ## Choose number of orderings
  if (p == 5) {
    ## For p = 5, there are at most 120 total orderings.
    ## Use all valid and all invalid orderings.
    n_valid_current <- n_valid_all
    n_invalid_current <- n_invalid_all
  } else {
    ## For p = 10, full enumeration is too expensive.
    ## Use a moderate random sample.
    n_valid_current <- min(n_valid_sample_p10, n_valid_all)
    n_invalid_current <- n_invalid_sample_p10
  }
  
  cat("Using n_valid =", n_valid_current, "\n", flush = TRUE)
  cat("Using n_invalid =", n_invalid_current, "\n", flush = TRUE)
  
  cat(
    "Using n_valid = ", n_valid_current, "\n",
    "Using n_invalid = ", n_invalid_current, "\n",
    file = progress_log,
    append = TRUE
  )
  
  for (n in n_grid) {
    
    cat("\nRunning n =", n, ", p =", p, "\n", flush = TRUE)
    
    cat(
      "\nRunning n = ", n, ", p = ", p, "\n",
      file = progress_log,
      append = TRUE
    )
    
    for (r in seq_len(n_rep)) {
      
      cat("  rep =", r, "started at", as.character(Sys.time()), "\n", flush = TRUE)
      
      cat(
        "  rep = ", r,
        " started at ", as.character(Sys.time()), "\n",
        file = progress_log,
        append = TRUE
      )
      
      set.seed(100000 + 1000 * p + 10 * n + r)
      
      time_one <- system.time({
        
        tmp <- one_run_theorem_sampled(
          n = n,
          E_true = E_true,
          alpha_true = alpha_true,
          beta_true = beta_true,
          delta_true = delta_true,
          gamma_true = gamma_true,
          true_topo = true_topo,
          L_G_mat = L_G_mat,
          n_valid = n_valid_current,
          n_invalid = n_invalid_current,
          c0 = c0,
          criterion = criterion,
          dist = dist
        )
        
      })
      
      cat(
        "  rep =", r,
        "| elapsed =", round(time_one["elapsed"], 2), "sec\n",
        flush = TRUE
      )
      
      cat(
        "  rep = ", r,
        " | elapsed = ", round(time_one["elapsed"], 2), " sec\n",
        file = progress_log,
        append = TRUE
      )
      
      tmp$p <- p
      tmp$K <- K
      tmp$n <- n
      tmp$rep <- r
      tmp$n_valid_used <- n_valid_current
      tmp$n_invalid_used <- n_invalid_current
      
      res_all[[counter]] <- tmp
      counter <- counter + 1L
      
      ## Save intermediate results after every repetition.
      save(
        res_all,
        counter,
        file = file.path(out_dir, "theorem_numeric_intermediate_res_all.RData")
      )
    }
  }
}

res_df <- bind_rows(res_all)

save(
  res_df,
  file = file.path(out_dir, "theorem_numeric_sampled_order_results.RData")
)

cat(
  "\nFinished all simulations at ", as.character(Sys.time()), "\n",
  file = progress_log,
  append = TRUE
)

############################################################
## 4. Statement (1): valid-vs-invalid score separation
############################################################
sep_df <- res_df %>%
  group_by(p, n, rep) %>%
  summarise(
    worst_valid_pen = min(PEN[valid]),
    best_invalid_pen = max(PEN[!valid]),
    gap = worst_valid_pen - best_invalid_pen,
    gap_scaled = gap / n,
    theorem1_success = gap > 0,
    .groups = "drop"
  )

summary_sep <- sep_df %>%
  group_by(p, n) %>%
  summarise(
    mean_gap = mean(gap),
    mean_gap_scaled = mean(gap_scaled),
    median_gap_scaled = median(gap_scaled),
    se_gap_scaled = sd(gap_scaled) / sqrt(n()),
    prob_theorem1_success = mean(theorem1_success),
    .groups = "drop"
  )

cat("\nSummary for score separation:\n")
print(summary_sep)

save(
  sep_df,
  summary_sep,
  file = file.path(out_dir, "theorem_score_gap_summary.RData")
)

############################################################
## 5. Figure 1: boxplot of score gap
############################################################
sep_df <- sep_df %>%
  mutate(
    p_factor = factor(p),
    n_factor = factor(n, levels = n_grid)
  )

P_gap <- ggplot(
  sep_df,
  aes(x = n_factor, y = gap_scaled, color = p_factor)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.8,
    color = "firebrick"
  ) +
  geom_boxplot(
    aes(group = interaction(n_factor, p_factor)),
    position = position_dodge(width = 0.75),
    width = 0.55,
    linewidth = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(group = interaction(n_factor, p_factor)),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 1.4,
    alpha = 0.75
  ) +
  scale_color_manual(
    values = c("5" = "#1f4e99", "10" = "#3a7f5f"),
    name = "p"
  ) +
  labs(
    x = "Sample size (n)",
    y = expression(gap[n] / n)
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = c(0.83, 0.18),
    legend.background = element_rect(color = "black", fill = "white"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

P_gap

ggsave(
  filename = file.path(out_dir, "theorem_score_gap_boxplot.pdf"),
  plot = P_gap,
  device = cairo_pdf,
  width = 7.5,
  height = 5.0,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "theorem_score_gap_boxplot.png"),
  plot = P_gap,
  width = 7.5,
  height = 5.0,
  dpi = 300
)

############################################################
## 6. Statement (2): graph recovery under valid orderings
############################################################
rec_df <- res_df %>%
  filter(valid) %>%
  group_by(p, n, rep) %>%
  summarise(
    recovery_rate = mean(recovered_exact),
    mean_edge_diff = mean(edge_diff),
    all_valid_recovered = all(recovered_exact),
    .groups = "drop"
  )

summary_rec <- rec_df %>%
  group_by(p, n) %>%
  summarise(
    mean_recovery_rate = mean(recovery_rate),
    mean_SHD = mean(mean_edge_diff),
    se_SHD = sd(mean_edge_diff) / sqrt(n()),
    prob_all_valid_recovered = mean(all_valid_recovered),
    .groups = "drop"
  )

cat("\nSummary for recovery under valid orderings:\n")
print(summary_rec)

save(
  rec_df,
  summary_rec,
  file = file.path(out_dir, "theorem_recovery_summary.RData")
)

############################################################
## 7. Figure 2: mean SHD with standard errors
############################################################
summary_rec <- summary_rec %>%
  mutate(
    p_factor = factor(p),
    n_factor = factor(n, levels = n_grid)
  )

P_shd <- ggplot(
  summary_rec,
  aes(x = n, y = mean_SHD, color = p_factor, group = p_factor)
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.4) +
  geom_errorbar(
    aes(
      ymin = mean_SHD - se_SHD,
      ymax = mean_SHD + se_SHD
    ),
    width = 50,
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = c("5" = "#1f4e99", "10" = "#3a7f5f"),
    name = "p"
  ) +
  scale_x_continuous(
    breaks = n_grid,
    labels = n_grid
  ) +
  labs(
    x = "Sample size (n)",
    y = expression(Delta)
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = c(0.83, 0.80),
    legend.background = element_rect(color = "black", fill = "white"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

P_shd

ggsave(
  filename = file.path(out_dir, "theorem_valid_order_SHD.pdf"),
  plot = P_shd,
  device = cairo_pdf,
  width = 7.0,
  height = 4.8,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "theorem_valid_order_SHD.png"),
  plot = P_shd,
  width = 7.0,
  height = 4.8,
  dpi = 300
)

############################################################
## 8. Save processed summaries
############################################################
save(
  res_df,
  sep_df,
  summary_sep,
  rec_df,
  summary_rec,
  file = file.path(out_dir, "theorem_numeric_processed_results.RData")
)

cat("\nAll done.\n", flush = TRUE)