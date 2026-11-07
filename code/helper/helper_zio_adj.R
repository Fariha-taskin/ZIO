suppressPackageStartupMessages(library(parallel, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(pscl, warn.conflicts = FALSE))

# -----------------------------------------------------------------------------
# Adjacent-transposition move in POSITION space.
# -----------------------------------------------------------------------------

adjacent_transposition <- function(ordering, i){
  # Swap the elements at adjacent positions i and i+1
  
  p <- length(ordering)
  
  if(i < 1 || i >= p){
    stop("i must be between 1 and p-1")
  }
  
  # make a copy
  ordering_next <- ordering
  
  # swap adjacent positions
  ordering_next[c(i, i+1)] <- ordering_next[c(i+1, i)]
  
  return(ordering_next)
}


# -----------------------------------------------------------------------------
# Worker-count helpers.
# - On schedulers, prefer the core count assigned by the job environment.
# - Otherwise, fall back to detectCores(logical = FALSE) when possible.
# -----------------------------------------------------------------------------
zip_available_cores <- function(logical = FALSE) {
  env_vars <- c("SLURM_CPUS_PER_TASK", "PBS_NP", "NSLOTS")
  for (nm in env_vars) {
    val <- suppressWarnings(as.integer(Sys.getenv(nm, unset = "")))
    if (!is.na(val) && val > 0L) {
      return(val)
    }
  }
  
  cores <- suppressWarnings(parallel::detectCores(logical = logical))
  if (is.na(cores) || cores < 1L) {
    cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
  }
  if (is.na(cores) || cores < 1L) {
    cores <- 1L
  }
  as.integer(cores)
}

zip_choose_workers <- function(total_tasks,
                               requested_workers = NULL,
                               leave_free = 4L,
                               logical = FALSE) {
  total_tasks <- as.integer(total_tasks)
  if (is.na(total_tasks) || total_tasks < 1L) {
    return(1L)
  }
  
  if (!is.null(requested_workers) && !is.na(requested_workers)) {
    return(max(1L, min(as.integer(requested_workers), total_tasks)))
  }
  
  available <- zip_available_cores(logical = logical)
  usable <- max(1L, available - as.integer(leave_free))
  max(1L, min(usable, total_tasks))
}

# -----------------------------------------------------------------------------
# Data normalization.
# -----------------------------------------------------------------------------
zip_normalize_data_list <- function(dat_list) {
  dat_list <- lapply(dat_list, function(X) {
    X <- as.data.frame(X)
    if (is.null(colnames(X))) {
      colnames(X) <- paste0("V", seq_len(ncol(X)))
    }
    X
  })
  
  ref_names <- colnames(dat_list[[1L]])
  for (k in seq_along(dat_list)) {
    if (!identical(colnames(dat_list[[k]]), ref_names)) {
      stop("All datasets must use the same column names / node labels.")
    }
  }
  dat_list
}

# -----------------------------------------------------------------------------
# ZIP model helpers.
# -----------------------------------------------------------------------------
zip_k_params <- function(fit, dist) {
  k <- length(coef(fit, model = "count")) + length(coef(fit, model = "zero"))
  if (identical(dist, "negbin")) {
    k <- k + 1L
  }
  as.integer(k)
}

zip_build_formula <- function(response, parents) {
  rhs <- if (length(parents)) paste(parents, collapse = " + ") else "1"
  stats::as.formula(paste0(response, " ~ ", rhs, " | ", rhs))
}

zip_fit_safe <- function(formula_obj, data, dist) {
  out <- try(pscl::zeroinfl(formula_obj, data = data, dist = dist), silent = TRUE)
  if (inherits(out, "try-error")) {
    return(NULL)
  }
  out
}

zip_get_ic <- function(fit, c0, dist, criterion = c("BIC", "AIC")) {
  criterion <- match.arg(criterion)
  loglik <- as.numeric(logLik(fit))
  k <- zip_k_params(fit, dist)
  
  if (identical(criterion, "AIC")) {
    penalty <- 2 * k
  } else {
    n_obs <- NROW(model.frame(fit))
    penalty <- k * (c0 * log(n_obs))
  }
  
  -2 * loglik + penalty
}

# -----------------------------------------------------------------------------
# Fit one node for one dataset under a fixed ordering position.
# This is the fundamental independent unit of work for parallelization.
# -----------------------------------------------------------------------------
zip_stepwise_fit_node <- function(X,
                                  topo,
                                  pos,
                                  c0,
                                  direction = c("both", "forward", "backward"),
                                  criterion = c("BIC", "AIC"),
                                  dist = c("poisson", "negbin"),
                                  trace = FALSE,
                                  max_steps = 1000L) {
  direction <- match.arg(direction)
  criterion <- match.arg(criterion)
  dist <- match.arg(dist)
  
  X <- as.data.frame(X)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("V", seq_len(ncol(X)))
  }
  vars <- colnames(X)
  
  pos <- as.integer(pos)
  if (pos < 1L || pos > length(topo)) {
    stop("'pos' must lie in 1:length(topo).")
  }
  
  response <- vars[topo[pos]]
  parents_all <- if (pos > 1L) vars[topo[seq_len(pos - 1L)]] else character(0)
  dat <- X[, c(response, parents_all), drop = FALSE]
  
  selected <- if (identical(direction, "backward")) parents_all else character(0)
  
  formula_best <- zip_build_formula(response, selected)
  fit_best <- zip_fit_safe(formula_best, dat, dist = dist)
  if (is.null(fit_best)) {
    stop(sprintf("Initial zeroinfl() failed for response '%s' at position %d.", response, pos))
  }
  best_ic <- zip_get_ic(fit_best, c0 = c0, dist = dist, criterion = criterion)
  
  if (length(parents_all) == 0L) {
    return(list(
      pos = pos,
      response = response,
      selected = character(0),
      IC = best_ic
    ))
  }
  
  if (isTRUE(trace)) {
    cat(sprintf("\n[pos=%d] start: %s | %s = %.4f\n",
                pos, deparse(formula_best), criterion, best_ic))
  }
  
  step <- 0L
  repeat {
    step <- step + 1L
    if (step > max_steps) {
      break
    }
    
    candidates <- list()
    outside <- setdiff(parents_all, selected)
    
    if (direction %in% c("forward", "both") && length(outside)) {
      for (v in outside) {
        candidates[[paste0("add:", v)]] <- c(selected, v)
      }
    }
    if (direction %in% c("backward", "both") && length(selected)) {
      for (v in selected) {
        candidates[[paste0("drop:", v)]] <- setdiff(selected, v)
      }
    }
    if (!length(candidates)) {
      break
    }
    
    tried <- lapply(names(candidates), function(lbl) {
      parent_set <- candidates[[lbl]]
      form <- zip_build_formula(response, parent_set)
      fit <- zip_fit_safe(form, dat, dist = dist)
      if (is.null(fit)) {
        return(NULL)
      }
      list(
        label = lbl,
        selected = parent_set,
        IC = zip_get_ic(fit, c0 = c0, dist = dist, criterion = criterion)
      )
    })
    tried <- Filter(Negate(is.null), tried)
    if (!length(tried)) {
      break
    }
    
    IC_values <- vapply(tried, `[[`, numeric(1L), "IC")
    best_move <- tried[[which.min(IC_values)]]
    
    if (isTRUE(trace)) {
      cat(sprintf("  best move: %-12s | %s = %.4f (current %.4f)\n",
                  best_move$label, criterion, best_move$IC, best_ic))
    }
    
    if (best_move$IC + 1e-8 < best_ic) {
      selected <- best_move$selected
      best_ic <- best_move$IC
    } else {
      break
    }
  }
  
  list(
    pos = pos,
    response = response,
    selected = selected,
    IC = best_ic
  )
}

# -----------------------------------------------------------------------------
# Global task queue.
# - Each task is one (dataset k, ordering position pos) fit.
# - Tasks are ordered by decreasing pos so the heavier regressions start first.
# -----------------------------------------------------------------------------
zip_make_tasks <- function(K, positions_by_k) {
  tasks <- vector("list", sum(lengths(positions_by_k)))
  cursor <- 1L
  
  for (k in seq_len(K)) {
    pos_vec <- sort(unique(as.integer(positions_by_k[[k]])), decreasing = TRUE)
    for (pos in pos_vec) {
      tasks[[cursor]] <- list(k = k, pos = pos)
      cursor <- cursor + 1L
    }
  }
  
  if (length(tasks) == 0L) {
    return(tasks)
  }
  
  ord <- order(vapply(tasks, function(x) x$pos, integer(1L)), decreasing = TRUE)
  tasks[ord]
}

# -----------------------------------------------------------------------------
# Cluster helpers.
# - Data are broadcast once to workers as ZIO_DAT_LIST.
# - This avoids re-serializing all K datasets at every hill-climbing step.
# -----------------------------------------------------------------------------
zip_make_cluster <- function(workers, dat_list) {
  workers <- as.integer(workers)
  if (workers < 1L) {
    stop("'workers' must be at least 1.")
  }
  
  dat_list <- zip_normalize_data_list(dat_list)
  
  cl <- parallel::makePSOCKcluster(workers)
  
  parallel::clusterEvalQ(cl, {
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1"
    )
    suppressPackageStartupMessages(library(pscl, warn.conflicts = FALSE))
    NULL
  })
  
  parallel::clusterExport(
    cl,
    varlist = c(
      "zip_k_params",
      "zip_build_formula",
      "zip_fit_safe",
      "zip_get_ic",
      "zip_stepwise_fit_node"
    ),
    envir = globalenv()
  )
  
  tmp_env <- list2env(list(ZIO_DAT_LIST = dat_list), parent = emptyenv())
  parallel::clusterExport(cl, varlist = "ZIO_DAT_LIST", envir = tmp_env)
  
  cl
}

zip_stop_cluster <- function(cl) {
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
  }
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Internal result assembly.
# -----------------------------------------------------------------------------
zip_init_state <- function(K, p, vars, base_IC = NULL, base_adj = NULL, base_selected = NULL) {
  if (is.null(base_IC)) {
    IC_mat <- matrix(NA_real_, nrow = K, ncol = p)
  } else {
    IC_mat <- base_IC
  }
  
  if (is.null(base_adj)) {
    adj_list <- replicate(
      K,
      matrix(0L, nrow = p, ncol = p, dimnames = list(vars, vars)),
      simplify = FALSE
    )
  } else {
    adj_list <- base_adj
    for (k in seq_len(K)) {
      storage.mode(adj_list[[k]]) <- "integer"
      dimnames(adj_list[[k]]) <- list(vars, vars)
    }
  }
  
  if (is.null(base_selected)) {
    selected_list <- replicate(K, vector("list", p), simplify = FALSE)
  } else {
    selected_list <- base_selected
  }
  
  list(IC_mat = IC_mat, adjacency = adj_list, selected = selected_list)
}

zip_evaluate_order <- function(dat_list,
                               topo,
                               c0,
                               cl = NULL,
                               positions_by_k = NULL,
                               base_IC = NULL,
                               base_adj = NULL,
                               base_selected = NULL,
                               direction = c("both", "forward", "backward"),
                               criterion = c("BIC", "AIC"),
                               dist = c("poisson", "negbin"),
                               trace = FALSE,
                               max_steps = 1000L) {
  direction <- match.arg(direction)
  criterion <- match.arg(criterion)
  dist <- match.arg(dist)
  
  dat_list <- zip_normalize_data_list(dat_list)
  K <- length(dat_list)
  p <- ncol(dat_list[[1L]])
  vars <- colnames(dat_list[[1L]])
  
  if (is.null(positions_by_k)) {
    positions_by_k <- replicate(K, seq_len(p), simplify = FALSE)
  }
  
  tasks <- zip_make_tasks(K = K, positions_by_k = positions_by_k)
  state <- zip_init_state(K = K, p = p, vars = vars,
                          base_IC = base_IC,
                          base_adj = base_adj,
                          base_selected = base_selected)
  
  if (length(tasks) == 0L) {
    total_bic <- sum(state$IC_mat, na.rm = FALSE)
    return(list(
      topo = topo,
      IC_mat = state$IC_mat,
      selected = state$selected,
      adjacency = state$adjacency,
      total_bic = total_bic,
      n_tasks = 0L
    ))
  }
  
  if (is.null(cl)) {
    task_results <- lapply(tasks, function(task) {
      fit <- zip_stepwise_fit_node(
        X = dat_list[[task$k]],
        topo = topo,
        pos = task$pos,
        c0 = c0,
        direction = direction,
        criterion = criterion,
        dist = dist,
        trace = trace,
        max_steps = max_steps
      )
      fit$k <- task$k
      fit
    })
  } else {
    task_results <- parallel::parLapplyLB(
      cl,
      tasks,
      fun = function(task, topo, c0, direction, criterion, dist, trace, max_steps) {
        fit <- zip_stepwise_fit_node(
          X = ZIO_DAT_LIST[[task$k]],
          topo = topo,
          pos = task$pos,
          c0 = c0,
          direction = direction,
          criterion = criterion,
          dist = dist,
          trace = trace,
          max_steps = max_steps
        )
        fit$k <- task$k
        fit
      },
      topo = topo,
      c0 = c0,
      direction = direction,
      criterion = criterion,
      dist = dist,
      trace = trace,
      max_steps = max_steps
    )
  }
  
  for (res in task_results) {
    k <- res$k
    pos <- res$pos
    response <- res$response
    chosen <- res$selected
    
    state$IC_mat[k, pos] <- res$IC
    state$selected[[k]][[pos]] <- chosen
    
    col_idx <- match(response, vars)
    if (is.na(col_idx)) {
      stop("Response variable not found while assembling adjacency.")
    }
    
    state$adjacency[[k]][, col_idx] <- 0L
    if (length(chosen)) {
      row_idx <- match(chosen, vars)
      if (anyNA(row_idx)) {
        stop("Selected parent not found while assembling adjacency.")
      }
      state$adjacency[[k]][row_idx, col_idx] <- 1L
    }
  }
  
  if (anyNA(state$IC_mat)) {
    stop("IC_mat contains NA values after evaluation. Check failed node fits.")
  }
  
  list(
    topo = topo,
    IC_mat = state$IC_mat,
    selected = state$selected,
    adjacency = state$adjacency,
    total_bic = sum(state$IC_mat),
    n_tasks = length(tasks)
  )
}

zip_affected_positions <- function(pos_i, pos_j) {
  seq.int(min(pos_i, pos_j), max(pos_i, pos_j))
}

zip_evaluate_move <- function(dat_list,
                              curr_topo,
                              pos_i,
                              pos_j = NULL,
                              curr_IC,
                              curr_adj,
                              curr_selected = NULL,
                              c0,
                              cl = NULL,
                              direction = c("both", "forward", "backward"),
                              criterion = c("BIC", "AIC"),
                              dist = c("poisson", "negbin"),
                              trace = FALSE,
                              max_steps = 1000L) {
  direction <- match.arg(direction)
  criterion <- match.arg(criterion)
  dist <- match.arg(dist)

  p <- length(curr_topo)

  ## pos_i is the adjacent-swap position.
  ## Swap positions pos_i and pos_i + 1.
  pos_i <- as.integer(pos_i)

  if (pos_i < 1L || pos_i >= p) {
    stop("For adjacent transposition, pos_i must be between 1 and p-1.")
  }

  pos_j <- pos_i + 1L

  next_topo <- adjacent_transposition(curr_topo, i = pos_i)

  affected <- zip_affected_positions(pos_i, pos_j)

  positions_by_k <- replicate(length(dat_list), affected, simplify = FALSE)

  fit <- zip_evaluate_order(
    dat_list = dat_list,
    topo = next_topo,
    c0 = c0,
    cl = cl,
    positions_by_k = positions_by_k,
    base_IC = curr_IC,
    base_adj = curr_adj,
    base_selected = curr_selected,
    direction = direction,
    criterion = criterion,
    dist = dist,
    trace = trace,
    max_steps = max_steps
  )

  fit$move_positions <- c(pos_i, pos_j)
  fit$affected_positions <- affected
  fit
}
# -----------------------------------------------------------------------------
# Backward-compatible single-dataset wrappers.
# These now run sequentially and are intentionally free of nested parallelism.
# -----------------------------------------------------------------------------
zip_stepwise_chain <- function(X,
                               topo,
                               c0,
                               direction = c("both", "forward", "backward"),
                               criterion = c("BIC", "AIC"),
                               dist = c("poisson", "negbin"),
                               trace = FALSE,
                               max_steps = 1000L,
                               parallel = FALSE,
                               workers = NULL,
                               n_pair = NULL,
                               keep_pairs_together = TRUE) {
  direction <- match.arg(direction)
  criterion <- match.arg(criterion)
  dist <- match.arg(dist)
  
  fit <- zip_evaluate_order(
    dat_list = list(X),
    topo = topo,
    c0 = c0,
    cl = NULL,
    direction = direction,
    criterion = criterion,
    dist = dist,
    trace = trace,
    max_steps = max_steps
  )
  
  list(
    models = vector("list", length(topo)),
    selected = fit$selected[[1L]],
    adjacency = fit$adjacency[[1L]],
    IC = fit$IC_mat[1L, ],
    criterion = criterion,
    dist = dist,
    parallel = FALSE,
    workers_used = 1L
  )
}

zip_stepwise_chain_next <- function(X,
                                    topo,
                                    i,
                                    j,
                                    IC_vec,
                                    c0,
                                    adj,
                                    direction = c("both", "forward", "backward"),
                                    criterion = c("BIC", "AIC"),
                                    dist = c("poisson", "negbin"),
                                    trace = FALSE,
                                    max_steps = 1000L,
                                    parallel = FALSE,
                                    workers = NULL,
                                    n_pair = NULL,
                                    keep_pairs_together = TRUE) {
  direction <- match.arg(direction)
  criterion <- match.arg(criterion)
  dist <- match.arg(dist)
  
  X <- as.data.frame(X)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("V", seq_len(ncol(X)))
  }
  p <- ncol(X)
  vars <- colnames(X)
  
  if (is.null(dim(adj)) || any(dim(adj) != c(p, p))) {
    adj <- matrix(0L, nrow = p, ncol = p, dimnames = list(vars, vars))
  } else {
    storage.mode(adj) <- "integer"
    dimnames(adj) <- list(vars, vars)
  }
  
  IC_mat <- matrix(as.numeric(IC_vec), nrow = 1L)
  fit <- zip_evaluate_move(
    dat_list = list(X),
    curr_topo = topo,
    pos_i = i,
    pos_j = j,
    curr_IC = IC_mat,
    curr_adj = list(adj),
    curr_selected = list(vector("list", p)),
    c0 = c0,
    cl = NULL,
    direction = direction,
    criterion = criterion,
    dist = dist,
    trace = trace,
    max_steps = max_steps
  )
  
  list(
    models = vector("list", length(topo)),
    selected = fit$selected[[1L]],
    adjacency = fit$adjacency[[1L]],
    IC = fit$IC_mat[1L, ],
    criterion = criterion,
    dist = dist
  )
}
