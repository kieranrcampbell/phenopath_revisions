
library(phenopath)
library(dpt)
library(monocle)
library(TSCAN)
library(aargh)

fit_phenopath <- function(exprs_mat, x) {
  fit <- phenopath(t(exprs_mat), x)
  return(fit)
}

fit_dpt <- function(exprs_mat) {
  pt <- Transitions(t(exprs_mat))
  DPT <- dpt(pt, branching = FALSE)
  DPT$DPT
}

fit_monocle2 <- function(exprs_mat) {
  cds <- newCellDataSet(exprs_mat)
  sizeFactors(cds) <- rep(1, ncol(cds))
  cds <- setOrderingFilter(cds, rownames(exprs_mat))
  cds <- reduceDimension(cds, norm_method = "none")
  cds <- orderCells(cds)
  cds$Pseudotime
}

fit_tscan <- function(exprs_mat) {
  tscan_pst <- rep(NA, ncol(exprs_mat))
  cl_data <- exprmclust(exprs_mat)
  tscan_order <- TSCANorder(cl_data, orderonly = FALSE)
  tscan_cell_inds <- match(tscan_order$sample_name, colnames(exprs_mat))
  tscan_pst[tscan_cell_inds] <- tscan_order$Pseudotime
  tscan_pst
}


pseudotime_inference <- function(algorithm = "phenopath",
                                 input_file = "sceset.rds",
                                 output_file = "myfile.csv",
                                 phenopath_fdata_file = "pdata.csv") {
  
  stopifnot(algorithm %in% c("phenopath", "dpt", "tscan", "monocle2"))
  
  sce <- readRDS(input_file)
  exprs_mat <- exprs(sce)
  
  pst <- NULL
  
  if(algorithm == "phenopath") {
    fit <- fit_phenopath(exprs_mat, sce$x)
    pst <- trajectory(fit)

    fdata_df <- data.frame(gene = featureNames(sce), 
                           m_beta = as.vector(fit$m_beta),
                           s_beta = as.vector(fit$s_beta))
    write.csv(fdata_df, phenopath_fdata_file)
  } else {
    pst <- switch(algorithm,
                  phenopath = fit_phenopath(exprs_mat, sce$x),
                  dpt = fit_dpt(exprs_mat),
                  monocle2 = fit_monocle2(exprs_mat),
                  tscan = fit_tscan(exprs_mat))  
  }
  output_df <- data.frame(pst = pst)
  write.csv(output_df, output_file)
}

aargh(pseudotime_inference)
