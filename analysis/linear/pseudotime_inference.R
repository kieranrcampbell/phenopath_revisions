
library(scater)
library(dpt)
library(monocle)
library(TSCAN)
library(aargh)
library(dplyr)
library(readr)

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

fit_pca <- function(exprs_mat) {
  pca <- prcomp(t(exprs_mat))
  return(pca$x[,1])
}

fit_linear_model <- function(pst_df, sce) {
  get_pval <- function(pst, y) {
    s <- summary(lm(y ~ pst, data = pst_df, na.omit = TRUE))
    s$coefficients[2,4]
  }
  
  apply(exprs(sce), 1, function(y) {
    apply(pst_df[,2:4], 2, get_pval, y)
  })
}


pseudotime_inference <- function(input_file = "sceset.rds",
                                 output_file = "myfile.csv",
                                 output_linear_model = "linmod.csv",
                                 dataset = "trapnell",
                                 hvg = 500) {
  
  sce <- readRDS(input_file)
  row_vars <- matrixStats::rowVars(exprs(sce))
  high_var <- row_vars >= sort(row_vars, decreasing = TRUE)[hvg]
  exprs_mat <- exprs(sce[high_var, ])
  
  pst_df <- data_frame(
    pca = fit_pca(exprs_mat),
    dpt = fit_dpt(exprs_mat),
    monocle2 = fit_monocle2(exprs_mat),
    tscan = fit_tscan(exprs_mat),
    dataset = dataset
  )
  
  lin_mod <- t(fit_linear_model(pst_df, sce[high_var, ]))
  lin_mod <- apply(lin_mod, 2, p.adjust, method = "BH")

  write_csv(pst_df, output_file)
}

aargh(pseudotime_inference)
