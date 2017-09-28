
library(aargh)
library(phenopath)
library(monocle)
library(scater)
library(tidyverse)
library(matrixStats)

fit_monocle2 <- function(exprs_mat) {
  cds <- newCellDataSet(exprs_mat)
  sizeFactors(cds) <- rep(1, ncol(cds))
  cds <- setOrderingFilter(cds, rownames(exprs_mat))
  cds <- reduceDimension(cds, norm_method = "none")
  cds <- orderCells(cds)
  cds$Pseudotime
}

scale_vec <- function(x) (x - mean(x)) / sd(x)

init_and_hypers <- function(input_sceset = "sce.rds",
                            output_csv = "output.csv",
                            control = FALSE,
                            z_init = "pc1",
                            elbo_tol = 1e-3,
                            tau_alpha = 1,
                            ab_beta_ratio = 1,
                            hvg = 1000) {
  sce <- readRDS(input_sceset)
  
  gene_vars <- rowVars(exprs(sce))
  var_cutoff <- sort(gene_vars, decreasing = TRUE)[hvg]
  sce_hvg <- sce[gene_vars >= var_cutoff, ]
  
  time_numeric <- as.numeric(gsub("h", "", sce_hvg$time))
  
  pst_init <- switch(z_init,
                     pc1 = 1,
                     monocle = scale_vec(fit_monocle2(exprs(sce_hvg))),
                     time = scale_vec(time_numeric))
  
  a_beta = ab_beta_ratio^2 / 10
  b_beta = ab_beta_ratio / 10
  
  fit <- output_df <- NULL
  if(control) {
    fit <- phenopath(sce_hvg, sce_hvg$x)
    output_df <- data_frame(
      pseudotime = trajectory(fit)
    )
  } else {
    fit <- phenopath(sce_hvg, sce_hvg$x,
                     elbo_tol = elbo_tol,
                     tau_alpha = tau_alpha,
                     z_init = pst_init,
                     a_beta = a_beta,
                     b_beta = b_beta)
    
    output_df <- data_frame(
      pseudotime = trajectory(fit),
      elbo_tol = elbo_tol,
      tau_alpha = tau_alpha,
      z_init = z_init,
      ab_beta_ratio = ab_beta_ratio
    )
  }
  
  write_csv(output_df, output_csv)
  
}

aargh(init_and_hypers)






