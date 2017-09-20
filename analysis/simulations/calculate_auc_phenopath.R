library(scater)
library(readr)
library(AUC)
library(aargh)
library(dplyr)
library(stringr)

source("analysis/simulations/parse_funcs.R")

sim_dir <- file.path("data", "simulations")

calculate_rocs_pp <- function(df_small) {
  df <- read_csv(df_small$path)
  quasi_p_val <- dnorm(0, df$m_beta, sqrt(df$s_beta))
  sce <- readRDS(df_small$sceset_path)
  is_interaction <- 1 * (fData(sce)$is_interaction)
  
  roc_obj <- AUC::roc(1 - quasi_p_val, factor(is_interaction))
  
  return(AUC::auc(roc_obj))
}

calculate_auc_pp <- function(output_file = "output.csv") {

  sim_dir <- file.path("data", "simulations")
  
  all_pp <- dir(file.path(sim_dir, "phenopath_fdata"))
  
  split <- str_split(all_pp, "_")
  
  df_split <- bind_rows(Map(parse_split, split, file.path(sim_dir, "phenopath_fdata", all_pp)))
  
  sceset_paths <- sapply(seq_len(nrow(df_split)), function(i) parse_sceset_path(df_split[i,]))
  
  df_split$sceset_path <- sceset_paths
  
  rocs <- sapply(seq_len(nrow(df_split)), function(i) {
    print(df_split[i,1:5])
    calculate_rocs_pp(df_split[i,])
  })
  
  df_split$auc <- rocs
  
  write_csv(df_split, output_file)
  
}

aargh(calculate_auc_pp)
