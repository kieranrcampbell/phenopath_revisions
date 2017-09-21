library(scater)
library(phenopath)
library(monocle)
library(readr)
library(dplyr)
library(aargh)
library(matrixStats)

fit_phenopath <- function(exprs_mat, x) {
  fit <- phenopath(t(exprs_mat), x)
  return(fit)
}

fit_monocle2 <- function(exprs_mat) {
  cds <- newCellDataSet(exprs_mat)
  sizeFactors(cds) <- rep(1, ncol(cds))
  cds <- setOrderingFilter(cds, rownames(exprs_mat))
  cds <- reduceDimension(cds, norm_method = "none")
  cds <- orderCells(cds)
  cds$Pseudotime
}

fit_hvg_pseudotime <- function(input_sceset = "input.rds",
                               algorithm = "phenopath",
                               hvg = 100,
                               dataset = "coad",
                               output_csv = "output.csv") {
  sce <- readRDS(input_sceset)
  gene_vars <- rowVars(exprs(sce))
  hvg_threshold <- sort(gene_vars, decreasing = TRUE)[hvg]
  is_hvg <- gene_vars >= hvg_threshold
  
  sce_hvg <- sce[is_hvg, ]
  
  pseudotime <- switch(algorithm,
                       phenopath = trajectory(fit_phenopath(exprs(sce_hvg), sce_hvg$x)),
                       monocle = fit_monocle2(exprs(sce_hvg)))
  
  sample_name = paste0("sample_", seq_len(ncol(sce_hvg)))
  
  output_dataframe <- data_frame(sample = sample_name, 
                                 pseudotime = pseudotime,
                                 algorithm = algorithm,
                                 hvg = hvg,
                                 dataset = dataset)
  write_csv(output_dataframe, output_csv)
}

aargh(fit_hvg_pseudotime)