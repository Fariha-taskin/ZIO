# ==============================================================================
# 1. Settings
# ==============================================================================

# ------------------------------------------------------------
# Portable project paths: run from ZIO/ or its parent directory
# ------------------------------------------------------------
project_dir <- "ZIO"
if (!dir.exists(file.path(project_dir, "code", "helper"))) {
  if (dir.exists(file.path("code", "helper"))) {
    project_dir <- "."
  } else {
    stop("Cannot find ZIO/code/helper. Run from ZIO/ or its parent.")
  }
}
helper_dir <- file.path(project_dir, "code", "helper")

# Set to "negbin" to summarize runs using ZINB-generated data.
input_model <- "poisson"
c0 <- 1L
base_dir <- file.path(project_dir, "output", "simulation", "model_misspecification",
                      paste0("Model_misspecification_", input_model, "_data"))

n <- 1000
p <- 25
K <- 1
seeds <- 1:30

fit_models <- c(ZIP = "poisson", ZINB = "negbin")
metrics <- c("TPR", "FDR", "SHD", "BIC")
digits <- c(TPR = 3, FDR = 3, SHD = 2, BIC = 2)

# Preserve the orientation adjustment in your original script.
transpose_true <- TRUE

if (K != 1) stop("This script summarizes K = 1 results only.")
if (!dir.exists(base_dir)) stop("Result directory does not exist: ", base_dir)

# ==============================================================================
# 2. Read results and calculate metrics
# ==============================================================================

read_distribution <- function(distribution) {
  result_dir <- file.path(
    base_dir,
    sprintf("simulation_result_n_%d_K_%d_P_%d_dist_%s",
            n, K, p, distribution)
  )
  
  rows <- list()
  missing_seeds <- integer(0)
  
  for (seed in seeds) {
    path <- file.path(
      result_dir, sprintf("p%s_K%s_seed%s.RData", p, K, seed)
    )
    
    if (!file.exists(path)) {
      missing_seeds <- c(missing_seeds, seed)
      next
    }
    
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    fit <- env$ZIO_DAG
    
    estimated <- fit$Adjacency_matrix[[1]]
    truth <- fit$out$E_true[[1]]
    if (transpose_true) truth <- t(truth)
    
    # Incomplete or nonbinary graphs should not silently affect the summary.
    if (!is.matrix(estimated) || !is.matrix(truth) ||
        !all(dim(estimated) == c(p, p)) ||
        !all(dim(truth) == c(p, p)) ||
        anyNA(estimated) || anyNA(truth) ||
        !all(estimated %in% c(0, 1)) ||
        !all(truth %in% c(0, 1)) ||
        any(diag(estimated) != 0) || any(diag(truth) != 0)) {
      stop("Invalid adjacency matrices in: ", path)
    }
    
    bic <- unlist(fit$BIC_score, use.names = FALSE)
    if (!is.numeric(bic) || length(bic) != 1L || !is.finite(bic)) {
      stop("Expected one finite BIC_score in: ", path)
    }
    
    TP <- sum(estimated == 1 & truth == 1)
    FP <- sum(estimated == 1 & truth == 0)
    true_edges <- sum(truth)
    estimated_edges <- sum(estimated)
    
    rows[[length(rows) + 1L]] <- data.frame(
      seed = seed,
      TPR = if (true_edges > 0) TP / true_edges else NA_real_,
      FDR = if (estimated_edges > 0) FP / estimated_edges else NA_real_,
      SHD = sum(abs(estimated - truth)),
      BIC = bic
    )
  }
  
  if (!length(rows)) {
    stop("No result files were found in: ", result_dir)
  }
  
  result <- do.call(rbind, rows)
  cat(sprintf("\n%s: loaded %d/%d seeds.\n",
              distribution, nrow(result), length(seeds)))
  if (length(missing_seeds)) {
    cat("Missing seeds:", paste(missing_seeds, collapse = ", "), "\n")
  }
  
  # FDR is undefined for an empty estimated graph, matching the original script.
  cat("Finite values used:",
      paste(sprintf("%s=%d", metrics,
                    vapply(result[metrics], function(x) sum(is.finite(x)),
                           integer(1))),
            collapse = ", "), "\n")
  
  result
}

mean_se <- function(x) {
  x <- x[is.finite(x)]
  c(
    mean = if (length(x)) mean(x) else NA_real_,
    se = if (length(x) > 1L) sd(x) / sqrt(length(x)) else NA_real_
  )
}

format_cell <- function(summary, digits) {
  if (!is.finite(summary["mean"])) return("NA")
  if (!is.finite(summary["se"])) {
    return(sprintf("%.*f (NA)", digits, summary["mean"]))
  }
  sprintf("%.*f (%.*f)", digits, summary["mean"], digits, summary["se"])
}

# ==============================================================================
# 3. Summarize both distributions
# ==============================================================================

results <- lapply(fit_models, read_distribution)
summaries <- lapply(results, function(x) {
  vapply(x[metrics], mean_se, c(mean = 0, se = 0))
})

table_data <- data.frame(
  Fitted_model = names(fit_models),
  stringsAsFactors = FALSE
)

for (metric in metrics) {
  table_data[[metric]] <- vapply(
    summaries,
    function(x) format_cell(x[, metric], digits[metric]),
    character(1)
  )
}

cat("\nMean (standard error):\n")
print(table_data, row.names = FALSE)

# ==============================================================================
# 4. Print and save the complete LaTeX table
# ==============================================================================

latex_rows <- apply(table_data, 1, function(row) {
  paste0(paste(row, collapse = " & "), " \\\\")
})

replication_text <- paste(
  sprintf("%s: %d simulations", names(results),
          vapply(results, nrow, integer(1))),
  collapse = "; "
)

caption <- sprintf(
  paste0(
    "Performance of ZIP and ZINB fitted models on %s-generated data with ",
    "$n=%d$, $p=%d$, and $K=%d$. ",
    "Entries report means (standard errors) across available repetitions (%s). ",
    "Undefined metrics are excluded from their respective summaries. ",
    "SHD is the adjacency-matrix Hamming distance, counting a reversed edge twice. ",
    "BIC scores are reported as stored in the result files."
  ),
  input_model, n, p, K, replication_text
)

latex_table <- paste(
  "\\begin{table}[!h]",
  "\\centering",
  "\\captionsetup{font={normalsize,stretch=1}}",
  paste0("\\caption{", caption, "}"),
  "\\label{tab:model_misspecification}",
  "\\begingroup",
  "\\renewcommand{\\arraystretch}{1.3}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "Fitted model & TPR & FDR & SHD & BIC score \\\\",
  "\\midrule",
  paste(latex_rows, collapse = "\n\\midrule\n"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\endgroup",
  "\\end{table}",
  sep = "\n"
)

cat("\n", latex_table, "\n", sep = "")

# Export one table and the numeric values for reproducibility.
output_file <- file.path(base_dir, "model_misspecification_table.tex")
writeLines(latex_table, output_file)
write.csv(table_data, file.path(base_dir, "model_misspecification_table.csv"), row.names = FALSE)
cat("\nSaved table to:", output_file, "\n")
