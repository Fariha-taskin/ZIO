## ---------------------------------------------------------------------
## 0. Working directory and source files
## ---------------------------------------------------------------------
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

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

p <- 10
K <- 1
n <- 1000
operator <- "adj"  # "r2r", "rts", "adj"
operator <- tolower(operator)
if (!operator %in% c("r2r", "rts", "adj")) stop("Choose r2r, rts, or adj")

helper_file <- switch(operator,
                      r2r = "helper_zio.R",
                      rts = "helper_zio_rts.R",
                      adj = "helper_zio_adj.R")
source(file.path(helper_dir, helper_file))
source(file.path(helper_dir, "data_generator_for_given_order.R"))
source(file.path(helper_dir, "metric.R"))

base_dir <- file.path(project_dir, "output", "simulation", "operator_comparison")
res_name <- if (operator == "r2r") {
  "final_simulation_results_fixed_data"
} else {
  paste0("final_simulation_results_fixed_data_", operator)
}
result_dir <- file.path(base_dir, res_name)

# ADJ uses set.seed(seed + 30) and saves that shifted seed in filenames.
reference_seed <- if (operator == "adj") 31L else 1L
reference_file <- file.path(
  result_dir,
  sprintf("simulation_result_K_%d_P_%d", K, p),
  sprintf("p%d_K%d_seed%d.RData", p, K, reference_seed)
)
if (!file.exists(reference_file)) stop("Missing reference result: ", reference_file)

env <- new.env(parent = globalenv())
load(reference_file, envir = env)
out <- env$ZIPBN_DAG$out

c0 <- 1
criterion <- "BIC"
dist <- "poisson"
direction <- "both"
max_steps <- 1000L

curr_order <- c(1:p)

init_fit <- zip_evaluate_order(
  dat_list = out$dat,
  topo = curr_order,
  c0 = c0,
  direction = direction,
  criterion = criterion,
  dist = dist,
  trace = FALSE,
  max_steps = max_steps
)

true_score <- init_fit$total_bic


## Folder containing checkpoint files
log_path <- file.path(result_dir, sprintf("log_result_K_%d_P_%d", K, p))
if (!dir.exists(log_path)) stop("Missing checkpoint directory: ", log_path)

## Settings
seeds <- 1:30

## Store BIC trajectory for each seed
BIC_score <- vector("list", length(seeds))
names(BIC_score) <- paste0("seed", seeds)

for (s in seeds) {
  stored_seed <- if (operator == "adj") s + 30L else s
  
  files_s <- list.files(
    path = log_path,
    pattern = paste0(
      "^ZIPBN_checkpoint_P_", p,
      "_K_", K,
      "_seed", stored_seed,
      "_iter[0-9]+\\.RData$"
    ),
    full.names = TRUE
  )
  
  if (length(files_s) == 0) {
    warning("No checkpoint found for seed ", stored_seed)
    next
  }
  
  ## Extract iteration number
  iter_s <- as.integer(
    sub(".*_iter([0-9]+)\\.RData$", "\\1", basename(files_s))
  )
  
  ## Choose latest checkpoint, e.g. iter500
  last_file <- files_s[which.max(iter_s)]
  
  ## Load checkpoint object
  checkpoint_env <- new.env(parent = globalenv())
  load(last_file, envir = checkpoint_env)
  checkpoint <- checkpoint_env$checkpoint
  if (is.null(checkpoint)) stop("Missing checkpoint object in ", last_file)
  
  ## Save full BIC trajectory for that seed
  BIC_score[[paste0("seed", s)]] <- unlist(checkpoint$BIC_score)
}

# ## Check
# length(BIC_score)
# length(BIC_score$seed1)
# head(BIC_score$seed1)
# tail(BIC_score$seed1)
# 

## Remove empty seeds, if any
BIC_score_clean <- BIC_score[sapply(BIC_score, length) > 0]

if (!length(BIC_score_clean)) stop("No checkpoint trajectories found.")

# Handle partial runs with different checkpoint lengths.
traj_list <- BIC_score_clean
df <- dplyr::bind_rows(lapply(seq_along(traj_list), function(i) {
  data.frame(
    run = paste0("run", i),
    iter = seq_along(traj_list[[i]]),
    logpost = as.numeric(traj_list[[i]])
  )
}))

# Clean any non-finite values
df <- df %>% filter(is.finite(logpost))

## ---- Target score and tolerance --------------------------------------------
Lstar <- sum(true_score)

tol <- max(1e-8, 1e-6 * abs(Lstar))

## ---- Distance from target --------------------------------------------------
## BIC is minimized, but for visualization we plot closeness to the target.
## gap = 0 means the trajectory reaches the target score.
df <- df %>%
  mutate(
    gap = abs(logpost - Lstar),
    near_star = gap <= tol,
    new_delta = 1 / (gap + 1)
  )

## ---- Find first hit per run ------------------------------------------------
hits <- df %>%
  group_by(run) %>%
  arrange(iter, .by_group = TRUE) %>%
  filter(near_star) %>%
  slice_min(iter, n = 1, with_ties = FALSE) %>%
  transmute(run, first_hit_iter = iter, hit_logpost = logpost) %>%
  ungroup()

## ---- Join hit info to df ---------------------------------------------------
df <- df %>%
  left_join(hits, by = "run") %>%
  mutate(ever_hit = !is.na(first_hit_iter))

## ---- Hit points ------------------------------------------------------------
hit_points2 <- df %>%
  filter(ever_hit, iter == first_hit_iter) %>%
  mutate(y_plot = 1) %>%
  distinct(run, iter, y_plot)

## ---- Better y-axis range ----------------------------------------------------
y_min <- min(df$new_delta[df$new_delta > 0], na.rm = TRUE)

## Round down to nearest power of 10
y_lower <- 10^floor(log10(y_min))

## Do not let it go too small unless needed
y_lower <- max(y_lower, 1e-6)

## y-axis breaks
y_breaks <- c(1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
y_breaks <- y_breaks[y_breaks >= y_lower]

y_labels <- lapply(y_breaks, function(x) {
  if (x == 1) {
    "0"
  } else {
    bquote(-10^.(as.integer(-log10(x))))
  }
})

## ---- Plot ------------------------------------------------------------------
P <- ggplot(df, aes(iter, new_delta, group = run)) +
  geom_line(color = "grey70", linewidth = 0.35, alpha = 0.9) +
  geom_point(
    data = hit_points2,
    aes(iter, y_plot, group = run),
    shape = 4,
    stroke = 1.0,
    size = 2.2,
    color = "red"
  ) +
  geom_hline(yintercept = 1, linewidth = 0.8, color = "blue") +
  scale_y_log10(
    limits = c(y_lower, 1.05),
    breaks = y_breaks,
    labels = y_labels
  ) +
  scale_x_continuous(
    limits = c(1, max(df$iter)),
    breaks = unique(c(1, floor(max(df$iter) / 2), max(df$iter)))
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 35) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.grid = element_blank(),
    plot.tag = element_text(size = 35),
    plot.tag.position = c(0.85, 0.15),
    axis.text.x = element_text(size = 35),
    axis.text.y = element_text(size = 35),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(4, "pt")
  )

P



figure_dir <- file.path(base_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(figure_dir, paste0("Traj_plot_", operator, ".pdf")),
  plot = P,
  device = cairo_pdf,
  width = 8.5,
  height = 10,
  dpi = 300
)
