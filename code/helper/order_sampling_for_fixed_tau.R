############################################################
# Numerically stable log-sum-exp
# Computes log(sum(exp(x))) without overflow/underflow
############################################################
logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(-Inf)   # if all entries are -Inf
  m + log(sum(exp(x - m)))
}

############################################################
# Build DP table for counting inversion sequences
# dp[i+1, s+1] = log(number of ways to assign inversions
#                    to first i elements summing to s)
############################################################
build_log_dp <- function(n, K) {
  dp <- matrix(-Inf, nrow = n + 1, ncol = K + 1)
  dp[1, 1] <- 0  # Base case: 0 elements, 0 inversions
  
  for (i in 1:n) {           # number of elements placed so far
    for (s in 0:K) {         # current inversion sum
      tmax <- min(i - 1, s)  # a_i can be at most i-1 and ≤ remaining s
      vals <- dp[i, (s - (0:tmax)) + 1]
      dp[i + 1, s + 1] <- logsumexp(vals)  # sum counts in log-scale
    }
  }
  dp
}

############################################################
# Sample an inversion sequence a = (a1,...,an)
# such that sum(a) = K using the DP table probabilities
############################################################
sample_a_given_K <- function(dp, n, K) {
  a <- integer(n)
  s <- K  # remaining inversion budget
  
  for (i in n:1) {                 # work backward
    tmax <- min(i - 1, s)           # allowed values for a[i]
    logw <- dp[i, (s - (0:tmax)) + 1]  # log weights of choices
    m <- max(logw)
    w <- exp(logw - m)              # convert to probabilities safely
    t <- sample.int(length(w), size = 1, prob = w) - 1
    a[i] <- t
    s <- s - t                      # update remaining inversions
  }
  a
}

############################################################
# Convert inversion sequence to permutation
# This is the Lehmer code → permutation transformation
############################################################
a_to_perm <- function(a) {
  n <- length(a)
  perm <- integer(0)
  
  for (i in 1:n) {
    val <- n - i + 1       # insert values n, n-1, ..., 1
    pos <- a[i] + 1        # insertion position
    if (pos == 1) {
      perm <- c(val, perm)
    } else if (pos == length(perm) + 1) {
      perm <- c(perm, val)
    } else {
      perm <- c(perm[1:(pos - 1)], val, perm[pos:length(perm)])
    }
  }
  perm
}

############################################################
# Sample permutation at fixed inversion count K
# Returns baseline reordered by that permutation
############################################################
sample_perm_fixed_K <- function(n , K, baseline = 1:n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(baseline) == n)
  
  M <- n * (n - 1) / 2     # maximum possible inversions
  if (K < 0 || K > M) stop(sprintf("K must be in [0, %d] for n=%d.", M, n))
  
  dp <- build_log_dp(n, K)
  a  <- sample_a_given_K(dp, n, K)
  p  <- a_to_perm(a)
  baseline[p]              # apply permutation to baseline ordering
}

############################################################
# Convert between Kendall tau and inversion count K
# tau = 1 - 4K / (n(n-1))
############################################################
tau_to_K <- function(tau, n, rule = c("nearest","floor","ceil")) {
  rule <- match.arg(rule)
  K_raw <- (1 - tau) * n * (n - 1) / 4
  
  K <- switch(rule,
              nearest = as.integer(round(K_raw)),
              floor   = as.integer(floor(K_raw)),
              ceil    = as.integer(ceiling(K_raw)))
  
  M <- n * (n - 1) / 2
  K <- max(0L, min(M, K))
  K
}

K_to_tau <- function(K, n) {
  1 - 4 * K / (n * (n - 1))
}

############################################################
# Sample permutation with (approximately) fixed Kendall tau
############################################################
sample_perm_fixed_tau <- function(n, tau, baseline = 1:n, seed = NULL, rule="nearest") {
  K <- tau_to_K(tau, n, rule = rule)
  perm <- sample_perm_fixed_K(n = n, K = K, baseline = baseline, seed = seed)
  list(perm = perm, K = K, tau_target = tau, tau_achieved = K_to_tau(K, n))
}

############################################################
# Generate permutations for tau = 0, 0.1, ..., 1
# R replicates per tau
############################################################
# p <- 5
# tau_grid <- seq(0, 1, by = 0.1)
# R <- 5
# baseline <- 1:p
# 
# perm_bank <- lapply(tau_grid, function(tau) {
#   reps <- lapply(1:R, function(r) {
#     sample_perm_fixed_tau(
#       n = p,
#       tau = tau,
#       baseline = baseline,
#       seed = 10000 + 100*r + round(100*tau)  # ensures reproducibility
#     )
#   })
#   reps
# })
# names(perm_bank) <- sprintf("tau_%.1f", tau_grid)
# 
# ############################################################
# # Check the *actual* Kendall tau achieved (first replicate)
# ############################################################
# achieved <- sapply(perm_bank, function(reps) reps[[1]]$tau_achieved)
# achieved
# 
# saveRDS(perm_bank, "order_sampling_for_fixed_tau.rds")
# 
# cor.test(1:20, perm_bank[["tau_0.0"]][[1]][["perm"]], method = "kendall")$estimate
