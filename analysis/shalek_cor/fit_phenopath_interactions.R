library(scater)
library(phenopath)
library(monocle)
library(readr)
library(dplyr)
library(aargh)
library(matrixStats)
library(TSCAN)
library(dpt)

fit_phenopath_interactions <- function(exprs_mat, x, pst_init = 1) {
  fit <- phenopath(t(exprs_mat), x, z_init = pst_init,
                   a_beta = 1e-2, b_beta = 1e-2)
  return(significant_interactions(fit))
}

fit_monocle2 <- function(exprs_mat) {
  cds <- newCellDataSet(exprs_mat)
  sizeFactors(cds) <- rep(1, ncol(cds))
  cds <- setOrderingFilter(cds, rownames(exprs_mat))
  cds <- reduceDimension(cds, norm_method = "none")
  cds <- orderCells(cds)
  cds$Pseudotime
}

fit_shalek_interactions <- function(input_sceset = "input.rds",
                                algorithm = "phenopath",
                                hvg = "all",
                               output_csv = "output.csv") {
  sce <- readRDS(input_sceset)
  is_hvg <- rep(TRUE, nrow(sce))
  if(hvg != "all") {
    hvg <- as.numeric(hvg)
    gene_vars <- rowVars(exprs(sce))
    hvg_threshold <- sort(gene_vars, decreasing = TRUE)[hvg]
    is_hvg <- gene_vars >= hvg_threshold
  }

  pst_init <- scale(as.numeric(gsub("h", "", sce$time)))[,1]  
  sce_hvg <- sce[is_hvg, ]
  
  interactions <- switch(algorithm,
                       phenopath_init_time = fit_phenopath_interactions(exprs(sce_hvg), sce_hvg$x, pst_init),
                       phenopath_init_pc1 = fit_phenopath_interactions(exprs(sce_hvg), sce_hvg$x),
                       phenopath_init_monocle = fit_phenopath_interactions(exprs(sce_hvg), sce_hvg$x, fit_monocle2(exprs(sce_hvg)))
                      )
  
  gene_name = paste0("gene_", seq_len(nrow(sce_hvg)))
  
  output_dataframe <- data_frame(gene = gene_name, 
                                 interactions = interactions,
                                 algorithm = algorithm,
                                 hvg = hvg)
  
  write_csv(output_dataframe, output_csv)
}

aargh(fit_shalek_interactions)