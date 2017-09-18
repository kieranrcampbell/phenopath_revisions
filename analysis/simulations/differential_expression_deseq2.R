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
  pseudotime <- pseudotime_df$pst
  
  cells_non_na <- which(!is.na(pseudotime))
  
  sce$pseudotime <- pseudotime
  sce$x <- factor(sce$x)
  
  count_mat <- counts(sce)
  coldata = pData(sce)
  
  dds <- DESeqDataSetFromMatrix(countData = count_mat,
                                colData = coldata, 
                                design = ~ pseudotime + x + pseudotime:x)
  dds <- DESeq(dds)
  res <- results(dds)
  
  interaction_qval <- res$padj
  
  output_data_frame <- data_frame(qval = interaction_qval)
  write_csv(output_data_frame, output_file)
}

aargh(dex_analysis_deseq2)

