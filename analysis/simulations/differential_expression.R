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
  
  interaction_qval <- NULL
  sce <- readRDS(input_sceset)
  pseudotime_df <- read_csv(pseudotime_file)
  
  if(all(is.na(pseudotime_df$pst))) {
    interaction_qval <- rep(NA, nrow(sce))
  } else {
    pseudotime <- scale(pseudotime_df$pst)[,1]
    cells_non_na <- which(!is.na(pseudotime))
    
    dmat <- dplyr::select(pData(sce), x) %>% 
      dplyr::mutate(pseudotime)
    
    design <- model.matrix(~ x + pseudotime + x:pseudotime, data = dmat[cells_non_na, ])
    
    
  
    
    dge <- DGEList(counts = counts(sce[, cells_non_na]))
    
    v <- voom(dge, design, plot = FALSE)
    
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    
    interaction_pval <- fit$p.value[,"x:pseudotime"]
    interaction_qval <- p.adjust(interaction_pval, method = "BH")
  }

  output_data_frame <- data_frame(qval = interaction_qval)
  write_csv(output_data_frame, output_file)
}

aargh(dex_analysis)

