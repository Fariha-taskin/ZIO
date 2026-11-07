## ---------------------------------------------------------------------
## 0. Working directory and source files
## ---------------------------------------------------------------------
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

source(file.path(helper_dir, "data_generator_from_order_zip_zeros.R"))
source(file.path(helper_dir, "order_sampling_for_fixed_tau.R"))

seed <- 1
set.seed(seed)

p <- 20
K <- 5
n <- 1000
zeros <- "50%"

## Choose target Kendall tau
tau_target <- 1   # change this as needed

out <- generate_data_from_order(
  n = n,
  p = p,
  K = K,
  topo = 1:p,
  zeros = zeros
)

X <- out$dat
G <- if (!is.null(out$E_true)) out$E_true else out$graph

## ------------------------------------------------------------
## 1) Sample one ordering per dataset k
## ------------------------------------------------------------
perm_bank <- lapply(1:K, function(k) {
  sample_perm_fixed_tau(
    n = p,
    tau = tau_target,
    baseline = 1:p
  )
})

## ------------------------------------------------------------
## 2) Inverse permutations
## ------------------------------------------------------------
inv_perm <- lapply(1:K, function(k) {
  order(perm_bank[[k]]$perm)
})

## ------------------------------------------------------------
## 3) Apply permutation to each DAG and data matrix
## ------------------------------------------------------------
G_perm <- vector("list", K)
X_perm <- vector("list", K)

for (k in 1:K) {
  idx <- inv_perm[[k]]
  
  G_perm[[k]] <- G[[k]][idx, idx]
  X_perm[[k]] <- X[[k]][, idx]
}

## ------------------------------------------------------------
## 4) Store sampled permutations
## ------------------------------------------------------------
perm <- vector("list", K)

for (k in 1:K) {
  perm[[k]] <- perm_bank[[k]]$perm
}

perm_mat <- do.call(cbind, perm)
colnames(perm_mat) <- paste0("sigma", 1:K)

perm_mat

cor_mat_kendall <- cor(perm_mat, method = "kendall")
mean_upper <- if (K > 1L) mean(cor_mat_kendall[upper.tri(cor_mat_kendall)]) else NA_real_

cat("For tau =", tau_target, "upper_triangle_ave =", round(mean_upper, 4))


# Save the generated result (single .RData file).
output_dir <- file.path(project_dir, "output", "simulation", "fixed_tau", "U_mapping")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
U_result <- list(
  seed = seed, n = n, p = p, K = K, zeros = zeros,
  tau_target = tau_target, mean_upper_tau = mean_upper,
  permutations = perm_mat, pairwise_tau = cor_mat_kendall
)
result_file <- file.path(
  output_dir, sprintf("U_p%d_K%d_tau%s_seed%d.RData", p, K, tau_target, seed)
)
save(U_result, file = result_file)
cat("\nSaved:", result_file, "\n")
