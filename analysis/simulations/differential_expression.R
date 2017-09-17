library(scater)
library(limma)
library(edgeR)
library(readr)
library(aargh)
library(magrittr)
library(dplyr)

dex_analysis <- function(input_sceset = "sce.rds",
                         pseudotime_file = "pseudotime.csv",
                         output_file = "qvals.csv") {
  
  sce <- readRDS(input_sceset)
  pseudotime_df <- read_csv(pseudotime_file)
  pseudotime <- pseudotime_df$pst
  
  dmat <- dplyr::select(pData(sce), x) %>% 
    dplyr::mutate(pseudotime)
  
  design <- model.matrix(~ x + pseudotime + x:pseudotime, data = dmat)
  
  

  
  dge <- DGEList(counts = counts(sce))
  
  v <- voom(dge, design, plot = FALSE)
  
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  interaction_pval <- fit$p.value[,"x:pseudotime"]
  interaction_qval <- p.adjust(interaction_pval, method = "BH")
  
  output_data_frame <- data_frame(qval = interaction_qval)
  write_csv(output_data_frame, output_file)
}

aargh(dex_analysis)

