# ----------------------------------------------------------------------
# 1. Project paths (no machine-specific drive letter)
# ----------------------------------------------------------------------
project_dir <- "ZIO"
if (!dir.exists(file.path(project_dir, "code", "helper"))) {
  if (dir.exists(file.path("code", "helper"))) {
    project_dir <- "."
  } else {
    stop("Cannot find the ZIO repository. Run from ZIO/ or its parent.")
  }
}
base_dir <- file.path(project_dir, "output", "simulation", "zero_proportions")

# ==============================================================================
# TPR and FDR under different zero proportions, with p = 20
# ==============================================================================

p <- 20
n_values <- c(250, 500, 1000, 2000)
K_values <- c(1, 5)
zeros <- c(25, 50, 75)
seeds <- 1:30

transpose_true <- TRUE
metrics <- c("TPR", "FDR")

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
    se = if (length(x) > 1L) sd(x) / sqrt(length(x)) else NA_real_,
    n = length(x)
  )
}

get_graph_list <- function(x, k, field, path) {
  if (is.matrix(x)) x <- list(x)
  
  if (!is.list(x) || length(x) != k) {
    stop(
      "Expected ", k, " matrices in ", field,
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
# 2. Read one (n, K, zero proportion) setting
# ==============================================================================

summarize_setting <- function(n, k, zero) {
  
  result_dir <- file.path(
    base_dir,
    "final_simulation_results",
    sprintf(
      "simulation_result_n_%d_K_%d_P_%d_zeros_%d%%",
      n, k, p, zero
    )
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
    
    # Calculate metrics for every dataset in this replicate.
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
          TPR = if (TP + FN > 0) TP / (TP + FN) else NA_real_,
          FDR = if (TP + FP > 0) FP / (TP + FP) else NA_real_
        )
      },
      FUN.VALUE = c(TPR = 0, FDR = 0)
    )
    
    # Average across all K datasets within the replicate.
    seed_results[i, ] <- apply(graph_metrics, 1, finite_mean)
    
    loaded <- c(loaded, seed)
  }
  
  # Calculate means and SEs across independent replicates.
  summaries <- vapply(
    metrics,
    function(metric) mean_se(seed_results[, metric]),
    FUN.VALUE = c(mean = 0, se = 0, n = 0)
  )
  
  cat(sprintf(
    "n = %d, K = %d, zeros = %d%%: loaded %d/%d replicates; TPR n = %d, FDR n = %d.\n",
    n, k, zero, length(loaded), length(seeds),
    summaries["n", "TPR"], summaries["n", "FDR"]
  ))
  
  if (length(missing)) {
    cat("  Missing seeds:", paste(missing, collapse = ", "), "\n")
  }
  
  data.frame(
    Zero = zero,
    K = k,
    n = n,
    Metric = metrics,
    Mean = unname(summaries["mean", ]),
    SE = unname(summaries["se", ]),
    N = unname(summaries["n", ]),
    row.names = NULL
  )
}

# ==============================================================================
# 3. Summarize all 24 settings
# ==============================================================================

settings <- expand.grid(
  n = n_values,
  K = K_values,
  Zero = zeros,
  KEEP.OUT.ATTRS = FALSE
)

summaries <- lapply(seq_len(nrow(settings)), function(i) {
  summarize_setting(
    n = settings$n[i],
    k = settings$K[i],
    zero = settings$Zero[i]
  )
})

summary_data <- do.call(rbind, summaries)

# ==============================================================================
# 4. Create the formatted table
# ==============================================================================

table_data <- do.call(
  rbind,
  lapply(zeros, function(zero) {
    data.frame(
      Zero = zero,
      K = rep(K_values, each = length(metrics)),
      Metric = rep(metrics, times = length(K_values)),
      stringsAsFactors = FALSE
    )
  })
)

rownames(table_data) <- NULL

for (sample_size in n_values) {
  
  table_data[[paste0("n=", sample_size)]] <- vapply(
    seq_len(nrow(table_data)),
    function(i) {
      
      row <- summary_data[
        summary_data$Zero == table_data$Zero[i] &
          summary_data$K == table_data$K[i] &
          summary_data$n == sample_size &
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
# 5. Generate the complete LaTeX table
# ==============================================================================

latex_rows <- character(0)
rows_per_zero <- length(K_values) * length(metrics)

for (zero in zeros) {
  
  for (k in K_values) {
    
    for (j in seq_along(metrics)) {
      
      metric <- metrics[j]
      
      row <- table_data[
        table_data$Zero == zero &
          table_data$K == k &
          table_data$Metric == metric,
      ]
      
      cells <- unlist(
        row[paste0("n=", n_values)],
        use.names = FALSE
      )
      
      zero_label <- if (k == K_values[1] && j == 1L) {
        paste0(
          "\\multirow{", rows_per_zero, "}{*}{$",
          zero, "\\%$}"
        )
      } else {
        ""
      }
      
      k_label <- if (j == 1L) {
        paste0("\\multirow{", length(metrics), "}{*}{", k, "}")
      } else {
        ""
      }
      
      latex_rows <- c(
        latex_rows,
        paste0(
          paste(
            c(zero_label, k_label, metric, cells),
            collapse = " & "
          ),
          " \\\\"
        )
      )
    }
  }
  
  if (zero != tail(zeros, 1)) {
    latex_rows <- c(latex_rows, "\\midrule")
  }
}

caption <- paste(
  sprintf(
    "Performance of the proposed method under different zero proportions, with $p=%d$.",
    p
  ),
  "Entries report means (standard errors) across available simulation replicates.",
  "Metrics are averaged across datasets within each replicate;",
  "undefined rates are excluded from their respective averages."
)

latex_table <- paste(
  "\\begin{table}[!t]",
  "\\centering",
  "\\captionsetup{font={normalsize,stretch=1}}",
  paste0("\\caption{", caption, "}"),
  "\\label{tab:zero_inflation_results}",
  "\\begingroup",
  "\\renewcommand{\\arraystretch}{1.3}",
  "\\begin{tabular}{ccccccc}",
  "\\toprule",
  "& & & \\multicolumn{4}{c}{$n$} \\\\",
  "\\cmidrule(lr){4-7}",
  paste0(
    "Zero proportion & $K$ & Metric & ",
    paste(n_values, collapse = " & "),
    " \\\\"
  ),
  "\\midrule",
  paste(latex_rows, collapse = "\n"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\endgroup",
  "\\end{table}",
  sep = "\n"
)

cat("\n", latex_table, "\n", sep = "")

# ----------------------------------------------------------------------
# 6. Save the table and its numerical summaries
# ----------------------------------------------------------------------
write.csv(
  summary_data,
  file.path(base_dir, "zero_proportion_metrics.csv"),
  row.names = FALSE
)
write.csv(
  table_data,
  file.path(base_dir, "zero_proportion_table.csv"),
  row.names = FALSE
)
writeLines(
  latex_table,
  file.path(base_dir, "zero_proportion_table.tex")
)
