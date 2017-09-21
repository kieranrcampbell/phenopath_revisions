library(stringr)

parse_split <- function(s, path) {
  N <- G <- p <- rep <- alg <- noise <- NULL
  for(i in seq_along(s)) {
    if(s[i] == "N") N <- as.numeric(s[i+1])
    if(s[i] == "G") G <- as.numeric(s[i+1])
    if(s[i] == "p") p <- as.numeric(s[i+1])
    if(s[i] == "rep") rep <- as.numeric(s[i+1])
    if(s[i] == "noise") {
      noise <- s[i+1]
    }
    if(s[i] == "alg") {
      ss <- s[i+1]
      alg <- str_split(ss, fixed("."))[[1]][1]
    }
  }
  data.frame(N, G, p, rep, alg, path, noise)
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
                                  "_noise_",
                                  df_small$noise,
                                  ".rds"))
  return(sceset_path)
}
