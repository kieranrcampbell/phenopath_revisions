library(scater)
library(phenopath)
library(readr)
library(dplyr)
library(aargh)
library(matrixStats)
library(TSCAN)
library(dpt)
library(reticulate)

fit_phenopath <- function(exprs_mat, x, pst_init = 1) {
  fit <- phenopath(t(exprs_mat), x, z_init = pst_init)
  return(fit)
}

fit_monocle2 <- function(exprs_mat) {
  library(monocle)
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


fit_wishbone <- function(exprs_mat, sce) {
  use_python("/apps/well/python/3.4.3/bin/python3")
  wishbone <- import('wishbone')
  csv_file <- tempfile()
  write.csv(t(exprs_mat), csv_file, row.names = TRUE)
  scdata <- wishbone$wb$SCData$from_csv(csv_file,
                                        data_type = 'sc-seq',
                                        normalize = FALSE)
  scdata$run_pca()
  # scdata$run_tsne()
  n <- ncol(exprs_mat)
  scdata$run_diffusion_map(knn = as.integer(n/4))
  wb <- wishbone$wb$Wishbone(scdata)
  
  root_cell <- colnames(sce)[which.min(sample(pData(sce)$pst))]
  
  nwp <- as.integer(round(n/5))
  kk <- as.integer(round(n/5))
  
  wishbone_pst = tryCatch(
    {
      wb$run_wishbone(
        start_cell = root_cell,
        branch = FALSE,
        k = kk,
        num_waypoints = nwp,
        components_list = c(1,2)
      )
      wishbone_pst <- wb$trajectory$as_matrix()
      wishbone_pst
    }, error = function(e) {
      message(e)
      rep(NA, n)
    })
  # plot(sce$pst, wishbone_pst)
  wishbone_pst
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
                       phenopath_init_monocle = trajectory(fit_phenopath(exprs(sce_hvg), sce_hvg$x, fit_monocle2(exprs(sce_hvg)))),
                       monocle2 = fit_monocle2(exprs(sce_hvg)),
                       dpt = fit_dpt(exprs(sce_hvg)),
                       tscan = fit_tscan(exprs(sce_hvg)),
                       wishbone = fit_wishbone(exprs(sce_hvg)))
  
  sample_name = paste0("sample_", seq_len(ncol(sce_hvg)))
  
  output_dataframe <- data_frame(sample = sample_name, 
                                 pseudotime = pseudotime,
                                 algorithm = algorithm,
                                 hvg = hvg)
  
  write_csv(output_dataframe, output_csv)
}

aargh(fit_shalek_pseudotime)