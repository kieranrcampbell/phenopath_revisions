library(scater)
library(DESeq2)
library(readr)
library(aargh)
library(magrittr)
library(dplyr)

dex_analysis_deseq2 <- function(input_sceset = "sce.rds",
                         pseudotime_file = "pseudotime.csv",
                         output_file = "qvals.csv") {
  
  sce <- readRDS(input_sceset)
  pseudotime_df <- read_csv(pseudotime_file)
  pseudotime <- scale(pseudotime_df$pst)[,1]
  
  cells_non_na <- which(!is.na(pseudotime))
  
  sce$pseudotime <- pseudotime
  sce$x <- factor(sce$x)
  
  count_mat <- counts(sce[, cells_non_na])
  coldata = pData(sce[, cells_non_na])
  
  interaction_qval <- tryCatch({
    dds <- DESeqDataSetFromMatrix(countData = count_mat,
                                  colData = coldata, 
                                  design = ~ pseudotime + x + pseudotime:x)
    dds <- DESeq(dds)
    res <- results(dds)
    interaction_qval <- res$padj
    interaction_qval
  },
  error = function(cond) {
    return(rep(NA, nrow(sce)))
  })
  
  output_data_frame <- data_frame(qval = interaction_qval)
  write_csv(output_data_frame, output_file)
}

aargh(dex_analysis_deseq2)

