
library(stringr)
library(readr)
library(aargh)
library(dplyr)

source("analysis/simulations/parse_funcs.R")

sim_dir <- file.path("data", "simulations")

correlate_pseudotimes <- function(df_small) {
  pseudotimes <- read_csv(df_small$path)$pst
  true_pseudotimes <- readRDS(df_small$sceset_path)$pst
  abs_kendall <- abs(cor(pseudotimes, true_pseudotimes, method = "kendall", use = "na"))
  return(abs_kendall)
}


parse_pseudotimes <- function(output_file = "output.csv") {
  
  
  all_pseudotimes <- dir(file.path(sim_dir, "pseudotimes"))
  
  split <- str_split(all_pseudotimes, "_")
  s <- split[1]
  
  df_split <- bind_rows(Map(parse_split, split, file.path(sim_dir, "pseudotimes", all_pseudotimes)))
  
  sceset_paths <- sapply(seq_len(nrow(df_split)), function(i) parse_sceset_path(df_split[i,]))
  
  df_split$sceset_path <- sceset_paths
  
  kendall_correlations <- sapply(seq_len(nrow(df_split)), function(i) {
    correlate_pseudotimes(df_split[i,])
    })
  
  df_split$kendall_correlation <- kendall_correlations
  
  write_csv(df_split, output_file)
}

aargh(parse_pseudotimes)


