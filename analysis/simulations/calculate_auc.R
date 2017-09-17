library(scater)
library(readr)
library(pROC)
library(aargh)
library(dplyr)
library(stringr)



sim_dir <- file.path("data", "simulations")

parse_split <- function(s, path) {
  N <- G <- p <- rep <- alg <- NULL
  for(i in seq_along(s)) {
    if(s[i] == "N") N <- as.numeric(s[i+1])
    if(s[i] == "G") G <- as.numeric(s[i+1])
    if(s[i] == "p") p <- as.numeric(s[i+1])
    if(s[i] == "rep") rep <- as.numeric(s[i+1])
    if(s[i] == "alg") {
      ss <- s[i+1]
      alg <- str_split(ss, fixed("."))[[1]][1]
    }
  }
  data.frame(N, G, p, rep, alg, path)
}

parse_sceset_path <- function(df_small) {
  sceset_path <- file.path(sim_dir, "scesets",
                           paste0("sceset_N_",
                                  df_small$N,
                                  "_G_",
                                  df_small$G,
                                  "_p_",
                                  df_small$p,
                                  "_rep_",
                                  df_small$rep,
                                  ".rds"))
  return(sceset_path)
}

calculate_rocs <- function(df_small) {
  qvals <- read_csv(df_small$path)$qval
  sce <- readRDS(df_small$sceset_path)
  is_interaction <- 1 * (fData(sce)$is_interaction)
  
  roc_obj <- roc(is_interaction, qvals)
  return(roc_obj$auc[1])
}

calculate_auc <- function(output_file = "output.csv") {

  sim_dir <- file.path("data", "simulations")
  
  all_qvals <- dir(file.path(sim_dir, "qvals"))
  
  split <- str_split(all_qvals, "_")
  
  df_split <- bind_rows(Map(parse_split, split, file.path(sim_dir, "qvals", all_qvals)))
  
  sceset_paths <- sapply(seq_len(nrow(df_split)), function(i) parse_sceset_path(df_split[i,]))
  
  df_split$sceset_path <- sceset_paths
  
  rocs <- sapply(seq_len(nrow(df_split)), function(i) {
    print(df_split[i,1:5])
    calculate_rocs(df_split[i,])
  })
  
  df_split$rocs <- roc
  
  write_csv(df_split, output_file)
  
}

aargh(calculate_auc)
