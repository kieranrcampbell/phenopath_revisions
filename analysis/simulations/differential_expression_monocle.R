library(scater)
library(readr)
library(aargh)
library(magrittr)
library(monocle)
library(dplyr)

dex_analysis_monocle <- function(input_sceset = "sce.rds",
                         pseudotime_file = "pseudotime.csv",
                         output_file = "qvals.csv",
                         random = 0) {
  
  sce <- readRDS(input_sceset)
  pseudotime_df <- read_csv(pseudotime_file)
  pseudotime <- scale(pseudotime_df$pst)[,1]
  
  cells_non_na <- which(!is.na(pseudotime))
  
  sce$pseudotime <- pseudotime
  sce <- sce[, cells_non_na]
  
  cds <- newCellDataSet(counts(sce), new("AnnotatedDataFrame", pData(sce)))
  cds$covariate <- cds$x
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  de <- differentialGeneTest(cds, 
                             fullModelFormulaStr = "~ covariate*sm.ns(pseudotime, df = 3)",
                             reducedModelFormulaStr = "~ covariate + sm.ns(pseudotime, df = 3)")
  

  interaction_qval <- de$qval
  
  output_data_frame <- data_frame(qval = interaction_qval)
  write_csv(output_data_frame, output_file)
}

aargh(dex_analysis_monocle)

