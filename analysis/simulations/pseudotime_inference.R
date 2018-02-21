
library(phenopath)
library(aargh)
library(reticulate)

fit_phenopath <- function(exprs_mat, x) {
  fit <- phenopath(t(exprs_mat), x)
  return(fit)
}

fit_dpt <- function(exprs_mat) {
  library(dpt)
  pt <- Transitions(t(exprs_mat))
  DPT <- dpt(pt, branching = FALSE)
  DPT$DPT
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

fit_tscan <- function(exprs_mat) {
  library(TSCAN)
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
  
  root_cell <- colnames(sce)[which.min((pData(sce)$pst))]
  
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


pseudotime_inference <- function(algorithm = "phenopath",
                                 input_file = "sceset.rds",
                                 output_file = "myfile.csv",
                                 phenopath_fdata_file = "pdata.csv") {
  
  stopifnot(algorithm %in% c("phenopath", "dpt", "tscan", "monocle2", "wishbone"))
  
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
                  tscan = fit_tscan(exprs_mat),
                  wishbone = fit_wishbone(exprs_mat, sce))  
  }
  output_df <- data.frame(pst = pst)
  write.csv(output_df, output_file)
}

aargh(pseudotime_inference)
