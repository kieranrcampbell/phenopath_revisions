library(scater)
library(phenopath)
library(monocle)
library(readr)
library(dplyr)
library(aargh)
library(matrixStats)
library(TSCAN)
library(dpt)

fit_phenopath <- function(exprs_mat, x, pst_init) {
  fit <- phenopath(t(exprs_mat), x, z_init = pst_init)
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


fit_dpt <- function(exprs_mat) {
  pt <- Transitions(t(exprs_mat))
  DPT <- dpt(pt, branching = FALSE)
  DPT$DPT
}

fit_tscan <- function(exprs_mat) {
  tscan_pst <- rep(NA, ncol(exprs_mat))
  cl_data <- exprmclust(exprs_mat)
  tscan_order <- TSCANorder(cl_data, orderonly = FALSE)
  tscan_cell_inds <- match(tscan_order$sample_name, colnames(exprs_mat))
  tscan_pst[tscan_cell_inds] <- tscan_order$Pseudotime
  tscan_pst
}

fit_shalek_pseudotime <- function(input_sceset = "input.rds",
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
  
  pseudotime <- switch(algorithm,
                       phenopath_init_time = trajectory(fit_phenopath(exprs(sce_hvg), sce_hvg$x, pst_init)),
                       phenopath_init_pc1 = trajectory(fit_phenopath(exprs(sce_hvg), sce_hvg$x)),
                       monocle2 = fit_monocle2(exprs(sce_hvg)),
                       dpt = fit_dpt(exprs(sce_hvg)),
                       tscan = fit_tscan(exprs(sce_hvg)))
  
  sample_name = paste0("sample_", seq_len(ncol(sce_hvg)))
  
  output_dataframe <- data_frame(sample = sample_name, 
                                 pseudotime = pseudotime,
                                 algorithm = algorithm,
                                 hvg = hvg)
  
  write_csv(output_dataframe, output_csv)
}

aargh(fit_shalek_pseudotime)