
library(stringr)
library(readr)
library(aargh)
library(dplyr)

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
    print(i)
    correlate_pseudotimes(df_split[i,])
    })
  
  df_split$kendall_correlation <- kendall_correlations
  
  write_csv(df_split, output_file)
}

aargh(parse_pseudotimes)


