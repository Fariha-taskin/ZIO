# Simulation 4.1: performance table for all p and K settings
project_dir <- "ZIO"
if (!dir.exists(file.path(project_dir, "code", "helper"))) {
  if (dir.exists(file.path("code", "helper"))) {
    project_dir <- "."
  } else {
    stop("Cannot locate ZIO/code/helper. Run from ZIO/ or its parent directory.")
  }
}
base_dir <- file.path(project_dir, "output", "simulation", "4.1")

# ==============================================================================
# Summarize simulation results for all P and K settings
# ==============================================================================
P_values <- c(5, 10, 20, 40)
K_values <- c(1, 5, 10, 20)
seeds <- 1:30

# Preserve the orientation adjustment in your original code.
transpose_true <- TRUE

metrics <- c("Delta", "tau", "TPR", "FDR")

if (!dir.exists(file.path(base_dir, "final_simulation_results"))) {
  stop("Simulation results directory does not exist: ",
       file.path(base_dir, "final_simulation_results"))
}

# ==============================================================================
# 1. Helper functions
# ==============================================================================

finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

mean_se <- function(x) {
  x <- x[is.finite(x)]
  
  c(
    mean = if (length(x)) mean(x) else NA_real_,
    se = if (length(x) > 1L) {
      sd(x) / sqrt(length(x))
    } else {
      NA_real_
    },
    n = length(x)
  )
}

get_graph_list <- function(x, k, field, path) {
  # Allow a single matrix when K = 1.
  if (is.matrix(x)) x <- list(x)
  
  if (!is.list(x) || length(x) != k) {
    stop(
      "Expected ", k, " adjacency matrices in ", field,
      ". Check the saved object structure in: ", path
    )
  }
  
  x
}

validate_graph <- function(x, p, path) {
  if (!is.matrix(x) ||
      !identical(dim(x), c(as.integer(p), as.integer(p))) ||
      anyNA(x) ||
      !all(x %in% c(0, 1)) ||
      any(diag(x) != 0)) {
    stop("Invalid binary adjacency matrix in: ", path)
  }
}

format_cell <- function(mean, se) {
  if (!is.finite(mean)) return("NA")
  
  if (!is.finite(se)) {
    return(sprintf("%.3f (NA)", mean))
  }
  
  sprintf("%.3f (%.3f)", mean, se)
}

# ==============================================================================
# 2. Read one (P, K) setting
# ==============================================================================

summarize_setting <- function(p, k) {
  
  result_dir <- file.path(
    base_dir,
    "final_simulation_results",
    sprintf("simulation_result_K_%d_P_%d", k, p)
  )
  
  seed_results <- matrix(
    NA_real_,
    nrow = length(seeds),
    ncol = length(metrics),
    dimnames = list(as.character(seeds), metrics)
  )
  
  loaded <- integer(0)
  missing <- integer(0)
  
  for (i in seq_along(seeds)) {
    seed <- seeds[i]
    
    path <- file.path(
      result_dir,
      sprintf("p%d_K%d_seed%d.RData", p, k, seed)
    )
    
    if (!file.exists(path)) {
      missing <- c(missing, seed)
      next
    }
    
    env <- new.env()
    load(path, envir = env)
    
    fit <- env$ZIPBN_DAG
    
    if (is.null(fit)) {
      stop("ZIPBN_DAG was not found in: ", path)
    }
    
    estimated_graphs <- get_graph_list(
      fit$Adjacency_matrix, k, "Adjacency_matrix", path
    )
    
    true_graphs <- get_graph_list(
      fit$out$E_true, k, "out$E_true", path
    )
    
    # Calculate metrics separately for every dataset.
    graph_metrics <- vapply(
      seq_len(k),
      function(subject) {
        
        estimated <- estimated_graphs[[subject]]
        truth <- true_graphs[[subject]]
        
        if (transpose_true) truth <- t(truth)
        
        validate_graph(estimated, p, path)
        validate_graph(truth, p, path)
        
        TP <- sum(estimated == 1 & truth == 1)
        FP <- sum(estimated == 1 & truth == 0)
        FN <- sum(estimated == 0 & truth == 1)
        
        c(
          Delta = FP + FN,
          TPR = if (TP + FN > 0) TP / (TP + FN) else NA_real_,
          FDR = if (TP + FP > 0) FP / (TP + FP) else NA_real_
        )
      },
      FUN.VALUE = c(Delta = 0, TPR = 0, FDR = 0)
    )
    
    # Average across the K datasets within this replicate.
    seed_results[i, c("Delta", "TPR", "FDR")] <-
      apply(graph_metrics, 1, finite_mean)
    
    # Joint ordering has one Kendall's tau per replicate.
    tau <- unlist(fit$tau[[1]], use.names = FALSE)
    
    if (!is.numeric(tau) ||
        length(tau) != 1L ||
        !is.finite(tau) ||
        abs(tau) > 1) {
      stop("Expected one finite Kendall's tau in: ", path)
    }
    
    seed_results[i, "tau"] <- tau
    loaded <- c(loaded, seed)
  }
  
  cat(sprintf(
    "P = %d, K = %d: loaded %d/%d replicates.\n",
    p, k, length(loaded), length(seeds)
  ))
  
  if (length(missing)) {
    cat("  Missing seeds:", paste(missing, collapse = ", "), "\n")
  }
  
  # Standard errors are calculated across replicates.
  summaries <- vapply(
    metrics,
    function(metric) mean_se(seed_results[, metric]),
    FUN.VALUE = c(mean = 0, se = 0, n = 0)
  )
  
  data.frame(
    K = k,
    P = p,
    Metric = metrics,
    Mean = unname(summaries["mean", ]),
    SE = unname(summaries["se", ]),
    N = unname(summaries["n", ]),
    row.names = NULL
  )
}

# ==============================================================================
# 3. Summarize all 16 settings
# ==============================================================================

settings <- expand.grid(
  K = K_values,
  P = P_values,
  KEEP.OUT.ATTRS = FALSE
)

summaries <- lapply(seq_len(nrow(settings)), function(i) {
  summarize_setting(
    p = settings$P[i],
    k = settings$K[i]
  )
})

summary_data <- do.call(rbind, summaries)

# Number of available replicates for each setting.
replicate_counts <- summary_data[
  summary_data$Metric == "Delta",
  c("K", "P", "N")
]

replicate_counts <- replicate_counts[
  order(replicate_counts$K, replicate_counts$P),
]

cat("\nAvailable replicates:\n")
print(replicate_counts, row.names = FALSE)

# ==============================================================================
# 4. Create the complete formatted table
# ==============================================================================

table_data <- data.frame(
  K = rep(K_values, each = length(metrics)),
  Metric = rep(metrics, times = length(K_values)),
  stringsAsFactors = FALSE
)

for (p in P_values) {
  
  table_data[[paste0("p=", p)]] <- vapply(
    seq_len(nrow(table_data)),
    function(i) {
      
      row <- summary_data[
        summary_data$K == table_data$K[i] &
          summary_data$P == p &
          summary_data$Metric == table_data$Metric[i],
      ]
      
      format_cell(row$Mean, row$SE)
    },
    character(1)
  )
}

cat("\nMean (standard error):\n")
print(table_data, row.names = FALSE)

# ==============================================================================
# 5. Generate the LaTeX table
# ==============================================================================

metric_labels <- c(
  Delta = "$\\Delta$",
  tau = "$\\tau$",
  TPR = "TPR",
  FDR = "FDR"
)

latex_rows <- character(0)

for (k in K_values) {
  
  for (j in seq_along(metrics)) {
    
    metric <- metrics[j]
    
    row <- table_data[
      table_data$K == k & table_data$Metric == metric,
    ]
    
    cells <- unlist(
      row[paste0("p=", P_values)],
      use.names = FALSE
    )
    
    k_label <- if (j == 1L) {
      sprintf("\\multirow{4}{*}{%d}", k)
    } else {
      ""
    }
    
    latex_rows <- c(
      latex_rows,
      paste0(
        paste(
          c(k_label, metric_labels[metric], cells),
          collapse = " & "
        ),
        " \\\\"
      )
    )
  }
  
  if (k != tail(K_values, 1)) {
    latex_rows <- c(latex_rows, "\\midrule")
  }
}

caption <- paste(
  "Graph and ordering recovery for varying $K$ and $p$.",
  "Entries report means (standard errors) across available",
  "simulation replicates.",
  "Graph metrics are averaged across datasets within each replicate;",
  "undefined rates are excluded from their respective averages."
)

latex_table <- paste(
  "\\begin{table}[!h]",
  "\\centering",
  "\\captionsetup{font={normalsize,stretch=1}}",
  paste0("\\caption{", caption, "}"),
  "\\label{tab:performance_multiple}",
  "\\begingroup",
  "\\small",
  "\\renewcommand{\\arraystretch}{1.3}",
  "\\begin{tabular}{clcccc}",
  "\\toprule",
  "$K$ & Metric & \\multicolumn{4}{c}{$p$} \\\\",
  "\\cmidrule(lr){3-6}",
  "& & 5 & 10 & 20 & 40 \\\\",
  "\\midrule",
  paste(latex_rows, collapse = "\n"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\endgroup",
  "\\end{table}",
  sep = "\n"
)

cat("\n", latex_table, "\n", sep = "")

# ==============================================================================
# 6. Save numerical summaries and LaTeX
# ==============================================================================

write.csv(
  summary_data,
  file.path(base_dir, "simulation_metrics_all_P_K.csv"),
  row.names = FALSE
)

write.csv(
  table_data,
  file.path(base_dir, "simulation_table_all_P_K.csv"),
  row.names = FALSE
)

writeLines(
  latex_table,
  file.path(base_dir, "simulation_table_all_P_K.tex")
)