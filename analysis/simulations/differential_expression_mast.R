library(scater)
library(DESeq2)
library(readr)
library(aargh)
library(magrittr)
library(MAST)
library(dplyr)

dex_analysis_mast <- function(input_sceset = "sce.rds",
                         pseudotime_file = "pseudotime.csv",
                         output_file = "qvals.csv") {
  
  interaction_qval <- NULL
  sce <- readRDS(input_sceset)
  pseudotime_df <- read_csv(pseudotime_file)
  
  if(all(is.na(pseudotime_df$pst))) {
    interaction_qval <- rep(NA, nrow(sce))
  } else {
    pst <- pseudotime_df$pst
    pseudotime <- (pst - min(pst, na.rm = TRUE)) # / (max(pst) - min(pst)) # Scale to [0,range(pst))
    
    cells_non_na <- which(!is.na(pseudotime))
    
    pdata <- data.frame(x = factor(sce$x), pseudotime)
    rownames(pdata) <- colnames(sce)
    
    sca <- FromMatrix(exprs(sce), pdata)
    fit <- zlm(~ x + pseudotime + x:pseudotime, sca[, cells_non_na])
  
    lrt <- lrTest(fit, "x:pseudotime")
    p_value <- lrt[,3,3]
    interaction_qval <- p.adjust(p_value, method = "BH")
  }
  output_data_frame <- data_frame(qval = interaction_qval)
  write_csv(output_data_frame, output_file)
}

aargh(dex_analysis_mast)

