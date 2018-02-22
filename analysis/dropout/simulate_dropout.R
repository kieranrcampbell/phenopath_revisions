library(SummarizedExperiment)
library(scater)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(polyester)
library(Biostrings)
library(readr)

source("analysis/simulations/simulate_funcs.R")

set.seed(12345L)

get_dropout_relationship <- function() {
  sce_rds <- "data/dropout/tec_sceset_qc.rds"
  sce <- readRDS(sce_rds)
  cc <- assayData(sce)[['counts']]
  cc <- cc[rowMeans(cc) > 0, ] # genes that are expressed at all
  
  pdrop_df <- data.frame(p_dropout = rowMeans(cc == 0), base_expression = rowMeans(cc))
  
  fit <- glm(p_dropout ~ base_expression, data = pdrop_df, family = "binomial")
  coefs <- coef(fit)
  
  pdrop <- function(x) 1 / (1 + exp(-(coefs[1] + coefs[2] * x)))  
  pdrop
}


tidy_sim_data <- function(mat, pst, x, sim_str) {
  df <- as_data_frame(mat)
  names(df) <- paste0(sim_str, "_gene_", seq_len(ncol(mat)))
  df <- mutate(df, pst, x, sim_str)
  gather(df, gene, expression, -pst, -x, -sim_str)
}

#' Simulate the mean function for sigmoidal interactions
#' @param N Number of cells
#' @param G Number of genes
#' @param prop_interaction Proportion of genes exhibiting interaction, in [0,1]
simulate_mean_function <- function(N, G, prop_interaction){
  pst <- runif(N)
  x <- sample(c(0, 1), N, replace = TRUE)
  
  G_pst <- round(G * (1 - prop_interaction))
  G_int <- G - G_pst
  G_class <- rep(round(G_int / 4), 3)
  G_class <- c(G_class, G_int - sum(G_class))
  
  k_lower <- 2
  k_upper <- 20
  sample_k <- function() runif(1, k_lower, k_upper)
  
  t0_lower <- 0.2
  t0_upper <- 0.8
  sample_t0 <- function() runif(1, t0_lower, t0_upper)
  
  mu0_lower <- 3
  mu0_upper <- 6
  sample_mu0 <- function() runif(1, mu0_lower, mu0_upper)

  
  y_pst_only <- replicate(G_pst, f_pseudotime_sample(pst, sample_k(), sample_t0(), sample_mu0()))
  df_pst_only <- tidy_sim_data(y_pst_only, pst, x, "pseudotime_only")
  
  y_interaction_1 <- replicate(G_class[1], fx_1(pst, x, sample_k(), sample_t0(), sample_mu0()))
  df_interaction_1 <- tidy_sim_data(y_interaction_1, pst, x, "int_type_1")
  
  y_interaction_2 <- replicate(G_class[2], fx_2(pst, x, sample_k(), sample_t0(), sample_mu0()))
  df_interaction_2 <- tidy_sim_data(y_interaction_2, pst, x, "int_type_2")
  
  y_interaction_3 <- replicate(G_class[3], fx_3(pst, x, sample_k(), sample_t0(), sample_mu0()))
  df_interaction_3 <- tidy_sim_data(y_interaction_3, pst, x, "int_type_3")
  
  y_interaction_4 <- replicate(G_class[4], fx_4(pst, x, sample_k(), sample_t0(), sample_mu0()))
  df_interaction_4 <- tidy_sim_data(y_interaction_4, pst, x, "int_type_4")
  
  bind_rows(
    df_pst_only,
    df_interaction_1,
    df_interaction_2,
    df_interaction_3,
    df_interaction_4
  )
  
}

mean_to_wide <- function(df) {
  ydf <- select(df, -sim_str) %>% 
    spread(gene, expression)
  ymat <- select(ydf, -pst, -x) %>% as.matrix() %>% t()
  colnames(ymat) <- paste0("cell_", seq_len(ncol(ymat)))
  
  is_interaction <- grepl("int", rownames(ymat))
  
  list(
    exprsmat = ymat,
    pst = ydf$pst,
    x = ydf$x,
    is_interaction = is_interaction
  )
}

#' Simulate the mean function for sigmoidal interactions
#' @param N Number of cells
#' @param G Number of genes
#' @param prop_interaction Proportion of genes exhibiting interaction, in [0,1]
#' @param replication The particular replication for saving files
#' @param ad Additional dropout (proportion of dropout to add to a dataset)
simulate_counts <- function(N, G, prop_interaction, replication, noise, ad) {
  
  # The signature string for this run, to be used in output files
  # So we can always parse what comes back
  sig_str <- paste0("ad_", ad,  
                    "_rep_", replication)
  
  df <- simulate_mean_function(N, G, prop_interaction)
  sim <- mean_to_wide(df)
  
  
  pos_gex <- 2^sim$exprsmat - 1 # Make expression positive
  count_mat <- sapply(seq_len(nrow(pos_gex)), function(i) {
    x <- pos_gex[i,]
    if(noise == "low") {
      return(NB(x, x / 3 + 1))
    } else {
      return(NB(x, 1))
    }
  })
  
  count_mat <- t(count_mat) - 1
  
  rownames(count_mat) <- rownames(sim$exprsmat)
  colnames(count_mat) <- colnames(sim$exprsmat)
  
  # Time to add in dropout -------
  
  pdrop <- get_dropout_relationship()
  
  # how many counts are non-zero?
  is_zero <- count_mat == 0
  original_prop_zero <- mean(is_zero)
  non_zero_counts <- sum(!is_zero)
  
  # of these we want to set proportion ad to 0
  num_to_zero <- ad * non_zero_counts
  
  potential_zero <- which(!is_zero)
  
  # probabilities of being set to zero - 
  p_zero <- pdrop(count_mat[potential_zero])
  
  to_zero <- sample(potential_zero, num_to_zero, prob = p_zero)
  
  count_mat[to_zero] <- 0
  
  new_prop_zero <- mean(count_mat == 0)
  
  if(any(count_mat < 0)) stop("Negative counts!")
  
  output_sce_file <- file.path("data", "dropout", "scesets", 
                           paste0("sceset_", sig_str, ".rds"))

  pdata <- data.frame(pst = sim$pst, x = sim$x,
                      original_prop_zero = original_prop_zero,
                      new_prop_zero = new_prop_zero,
                      proportion_set_to_zero = ad)
  
  rownames(pdata) <- colnames(count_mat)
  fdata <- data.frame(is_interaction = sim$is_interaction)
  rownames(fdata) <- rownames(count_mat)
  
  bioc_version <- tools:::.BioC_version_associated_with_R_version() 
  sce <- NULL
  
  if(bioc_version > 3.5) {
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(assays = list(counts = count_mat), 
              colData = pdata,
              rowData = fdata)
    sce <- normaliseExprs(sce)
    } else {
    sce <- newSCESet(countData = count_mat, 
                     phenoData = new("AnnotatedDataFrame", pdata),
                     featureData = new("AnnotatedDataFrame", fdata))
  }
  
  saveRDS(sce, output_sce_file)
}

## Make sure these line up with the snakefile

N <- 200
G <- 500

# prop_interactions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
pi <- 0.2
reps <- 40
# noises <- c("low", "high")
noise <- "high"
additional_dropout <- c(0.05, 0.1, 0.2, 0.5, 0.8, 0.9)

for(r in seq_len(reps)) {
  for(ad in additional_dropout) {
    simulate_counts(N, G, pi, r, noise, ad)
  }
}




