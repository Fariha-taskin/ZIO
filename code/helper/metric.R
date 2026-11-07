############### Edge-level metrics ###############

check_edge <- function(adj_pre, adj_act) {
  sum(abs(adj_pre - adj_act))
}

TP <- function(adj_pre, adj_act) {
  true_edges <- which(adj_act == 1)
  pred_edges <- which(adj_pre == 1)
  length(intersect(true_edges, pred_edges))
}

FP <- function(adj_pre, adj_act) {
  true_non_edges <- which(adj_act == 0)
  pred_edges <- which(adj_pre == 1)
  length(intersect(true_non_edges, pred_edges))
}

FN <- function(adj_pre, adj_act) {
  true_edges <- which(adj_act == 1)
  pred_non_edges <- which(adj_pre == 0)
  length(intersect(true_edges, pred_non_edges))
}

FPR <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  true_non_edges <- which(adj_act == 0)
  pred_edges <- which(adj_pre == 1)
  length(intersect(true_non_edges, pred_edges)) / length(true_non_edges)
}

TPR <- function(adj_pre, adj_act) {
  tp <- TP(adj_pre, adj_act)
  fn <- FN(adj_pre, adj_act)
  if ((tp + fn) == 0) return(0)
  tp / (tp + fn)
}

FDR <- function(adj_pre, adj_act) {
  tp <- TP(adj_pre, adj_act)
  fp <- FP(adj_pre, adj_act)
  if ((tp + fp) == 0) return(0)
  fp / (tp + fp)
}

############### Undirected mismatch count ###############
check_adj_l2 <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  sum((adj_pre - adj_act)^2) / 2
}
