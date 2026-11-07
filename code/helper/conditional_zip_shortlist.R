## Conditional ZIP-score shortlist override for helper_zio.R.
##
## Source this file immediately after helper_zio.R.  It replaces only
## zip_stepwise_fit_node(); the order search, task queue, adjacency assembly,
## and information criterion remain unchanged.
##
## Configuration (environment variables):
##   ZIO_PARENT_SELECTION = shortlist (default) or full
##   ZIO_SHORTLIST_B      = number of fully fitted additions (default 10)
##   ZIO_DMAX             = maximum selected parents (default 5; Inf allowed)
##
## In shortlist mode, every forward round:
##   1. scores all remaining predecessors using a vectorized two-component
##      ZIP score conditional on the currently selected parents;
##   2. fully fits only the top b additions;
##   3. fully fits every possible backward deletion;
##   4. accepts the move with the largest improvement in the original BIC/AIC.
##
## All fitted supports use the original unpenalized pscl::zeroinfl likelihood.

zip_stepwise_fit_node <- function(
    X,
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

  selection_mode <- tolower(Sys.getenv(
    "ZIO_PARENT_SELECTION",
    unset = "shortlist"
  ))
  if (!selection_mode %in% c("shortlist", "full")) {
    stop("ZIO_PARENT_SELECTION must be 'shortlist' or 'full'.")
  }
  if (selection_mode == "shortlist" && dist != "poisson") {
    stop("Conditional ZIP-score shortlist currently supports dist='poisson'.")
  }

  top_b <- suppressWarnings(as.integer(Sys.getenv(
    "ZIO_SHORTLIST_B",
    unset = "10"
  )))
  if (is.na(top_b) || top_b < 1L) {
    stop("ZIO_SHORTLIST_B must be a positive integer.")
  }

  dmax_text <- Sys.getenv("ZIO_DMAX", unset = "5")
  if (tolower(dmax_text) %in% c("inf", "infinity")) {
    dmax_requested <- Inf
  } else {
    dmax_requested <- suppressWarnings(as.numeric(dmax_text))
  }
  if (length(dmax_requested) != 1L || is.na(dmax_requested) ||
      dmax_requested < 0 ||
      (!is.infinite(dmax_requested) &&
       dmax_requested != floor(dmax_requested))) {
    stop("ZIO_DMAX must be a nonnegative integer or Inf.")
  }

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
  parents_all <- if (pos > 1L) {
    vars[topo[seq_len(pos - 1L)]]
  } else {
    character(0)
  }
  dat <- X[, c(response, parents_all), drop = FALSE]

  ## Full mode reproduces the original candidate breadth and has no dmax cap.
  dmax <- if (selection_mode == "full" || is.infinite(dmax_requested)) {
    length(parents_all)
  } else {
    min(as.integer(dmax_requested), length(parents_all))
  }

  cache <- new.env(hash = TRUE, parent = emptyenv())
  full_zip_fits <- 0L
  cache_hits <- 0L
  scored_candidates <- 0L
  shortlisted_candidates <- 0L

  normalize_support <- function(parent_set) {
    parents_all[parents_all %in% unique(as.character(parent_set))]
  }

  support_key <- function(parent_set) {
    parent_set <- normalize_support(parent_set)
    if (!length(parent_set)) "<empty>" else paste(parent_set, collapse = "\r")
  }

  mapped_start <- function(old_fit, new_parent_set) {
    if (is.null(old_fit)) return(NULL)
    count_old <- try(stats::coef(old_fit, model = "count"), silent = TRUE)
    zero_old <- try(stats::coef(old_fit, model = "zero"), silent = TRUE)
    if (inherits(count_old, "try-error") ||
        inherits(zero_old, "try-error")) {
      return(NULL)
    }

    count_new <- c("(Intercept)" = unname(count_old["(Intercept)"]))
    zero_new <- c("(Intercept)" = unname(zero_old["(Intercept)"]))
    for (parent in new_parent_set) {
      count_new[parent] <- if (parent %in% names(count_old)) {
        count_old[parent]
      } else {
        0
      }
      zero_new[parent] <- if (parent %in% names(zero_old)) {
        zero_old[parent]
      } else {
        0
      }
    }
    if (any(!is.finite(c(count_new, zero_new)))) return(NULL)
    list(count = unname(count_new), zero = unname(zero_new))
  }

  fit_support <- function(parent_set, start_fit = NULL) {
    parent_set <- normalize_support(parent_set)
    key <- support_key(parent_set)
    if (exists(key, envir = cache, inherits = FALSE)) {
      cache_hits <<- cache_hits + 1L
      return(get(key, envir = cache, inherits = FALSE))
    }

    formula_obj <- zip_build_formula(response, parent_set)
    start <- if (selection_mode == "shortlist") {
      mapped_start(start_fit, parent_set)
    } else {
      NULL
    }

    call_fit <- function(start_value) {
      full_zip_fits <<- full_zip_fits + 1L
      out <- try(
        pscl::zeroinfl(
          formula_obj,
          data = dat,
          dist = dist,
          control = pscl::zeroinfl.control(
            maxit = 10000L,
            start = start_value
          )
        ),
        silent = TRUE
      )
      if (inherits(out, "try-error")) NULL else out
    }

    fit <- call_fit(start)
    if (is.null(fit) && !is.null(start)) {
      fit <- call_fit(NULL)
    }

    result <- if (is.null(fit)) {
      list(selected = parent_set, IC = Inf, fit = NULL)
    } else {
      IC <- try(
        zip_get_ic(
          fit,
          c0 = c0,
          dist = dist,
          criterion = criterion
        ),
        silent = TRUE
      )
      if (inherits(IC, "try-error") || !is.finite(IC)) {
        list(selected = parent_set, IC = Inf, fit = fit)
      } else {
        list(selected = parent_set, IC = as.numeric(IC), fit = fit)
      }
    }
    assign(key, result, envir = cache)
    result
  }

  conditional_scores <- function(current, candidates, selected) {
    candidates <- as.character(candidates)
    if (!length(candidates)) return(numeric(0))
    fit <- current$fit
    if (is.null(fit)) {
      return(stats::setNames(rep(-Inf, length(candidates)), candidates))
    }

    y <- fit$y
    if (is.null(y)) {
      y <- stats::model.response(stats::model.frame(fit))
    }
    y <- as.numeric(y)
    lambda <- as.numeric(stats::predict(fit, type = "count"))
    pi_zero <- as.numeric(stats::predict(fit, type = "zero"))
    lambda <- pmin(pmax(lambda, 1e-10), 1e10)
    pi_zero <- pmin(pmax(pi_zero, 1e-8), 1 - 1e-8)

    ## Posterior probability that an observed zero is structural:
    ## pi / {pi + (1-pi) exp(-lambda)}.
    posterior_structural <- numeric(length(y))
    zero_idx <- y == 0
    posterior_structural[zero_idx] <- stats::plogis(
      stats::qlogis(pi_zero[zero_idx]) + lambda[zero_idx]
    )

    count_residual <- (1 - posterior_structural) * (y - lambda)
    zero_residual <- posterior_structural - pi_zero
    X_candidate <- as.matrix(X[, candidates, drop = FALSE])

    ## Residualize candidates against the intercept and selected parents so
    ## that ranking is conditional rather than marginal.
    if (length(selected)) {
      nuisance <- cbind(1, as.matrix(X[, selected, drop = FALSE]))
      X_candidate <- qr.resid(qr(nuisance), X_candidate)
    } else {
      X_candidate <- scale(
        X_candidate,
        center = TRUE,
        scale = FALSE
      )
    }

    scale_sq <- colSums(X_candidate^2)
    valid <- is.finite(scale_sq) & scale_sq > 1e-12
    score <- rep(-Inf, length(candidates))
    if (any(valid)) {
      Xv <- X_candidate[, valid, drop = FALSE]
      U_count <- as.numeric(crossprod(Xv, count_residual))
      U_zero <- as.numeric(crossprod(Xv, zero_residual))
      I_count <- as.numeric(crossprod(
        Xv^2,
        (1 - posterior_structural) * lambda
      ))
      I_zero <- as.numeric(crossprod(
        Xv^2,
        pi_zero * (1 - pi_zero)
      ))
      score[valid] <-
        U_count^2 / pmax(I_count, 1e-10) +
        U_zero^2 / pmax(I_zero, 1e-10)
    }
    names(score) <- candidates
    score
  }

  selected <- if (identical(direction, "backward")) {
    parents_all
  } else {
    character(0)
  }
  current <- fit_support(selected)
  if (!is.finite(current$IC)) {
    stop(sprintf(
      "Initial zeroinfl() failed for response '%s' at position %d.",
      response,
      pos
    ))
  }

  if (!length(parents_all)) {
    return(list(
      pos = pos,
      response = response,
      selected = character(0),
      IC = current$IC,
      selection_method = selection_mode,
      full_zip_fits = full_zip_fits,
      cache_hits = cache_hits,
      scored_candidates = 0L,
      shortlisted_candidates = 0L
    ))
  }

  if (isTRUE(trace)) {
    cat(sprintf(
      "\n[pos=%d] start | mode=%s | b=%d | dmax=%s | %s=%.4f\n",
      pos,
      selection_mode,
      top_b,
      if (selection_mode == "full") "Inf" else as.character(dmax),
      criterion,
      current$IC
    ))
  }

  for (step in seq_len(max_steps)) {
    tried <- list()
    outside <- setdiff(parents_all, selected)

    if (direction %in% c("forward", "both") &&
        length(outside) &&
        length(selected) < dmax) {
      additions <- outside
      if (selection_mode == "shortlist" && length(outside) > top_b) {
        score <- conditional_scores(current, outside, selected)
        scored_candidates <- scored_candidates + length(outside)
        ranked <- order(-score, names(score), na.last = TRUE)
        additions <- names(score)[head(ranked, top_b)]
      } else if (selection_mode == "shortlist") {
        scored_candidates <- scored_candidates + length(outside)
      }
      shortlisted_candidates <-
        shortlisted_candidates + length(additions)

      for (candidate in additions) {
        result <- fit_support(
          c(selected, candidate),
          start_fit = current$fit
        )
        if (is.finite(result$IC)) {
          result$label <- paste0("add:", candidate)
          tried[[length(tried) + 1L]] <- result
        }
      }
    }

    ## Backward deletions remain exhaustive.
    if (direction %in% c("backward", "both") && length(selected)) {
      for (candidate in selected) {
        result <- fit_support(
          setdiff(selected, candidate),
          start_fit = current$fit
        )
        if (is.finite(result$IC)) {
          result$label <- paste0("drop:", candidate)
          tried[[length(tried) + 1L]] <- result
        }
      }
    }

    if (!length(tried)) break
    IC_values <- vapply(tried, `[[`, numeric(1L), "IC")
    best_move <- tried[[which.min(IC_values)]]

    if (isTRUE(trace)) {
      cat(sprintf(
        "  best move: %-12s | %s=%.4f (current %.4f)\n",
        best_move$label,
        criterion,
        best_move$IC,
        current$IC
      ))
    }

    if (best_move$IC + 1e-8 < current$IC) {
      current <- best_move
      selected <- current$selected
    } else {
      break
    }
  }

  list(
    pos = pos,
    response = response,
    selected = selected,
    IC = current$IC,
    selection_method = selection_mode,
    shortlist_b = top_b,
    dmax = if (selection_mode == "full") Inf else dmax,
    full_zip_fits = full_zip_fits,
    cache_hits = cache_hits,
    scored_candidates = scored_candidates,
    shortlisted_candidates = shortlisted_candidates
  )
}

message(sprintf(
  paste0(
    "Loaded conditional ZIP-score parent selection: mode=%s, b=%s, ",
    "dmax=%s"
  ),
  Sys.getenv("ZIO_PARENT_SELECTION", unset = "shortlist"),
  Sys.getenv("ZIO_SHORTLIST_B", unset = "10"),
  Sys.getenv("ZIO_DMAX", unset = "5")
))
