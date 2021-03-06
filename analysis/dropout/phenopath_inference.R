library(phenopath)
library(aargh)
library(scater)

fit_phenopath <- function(exprs_mat, x) {
  fit <- phenopath(t(exprs_mat), x)
  return(fit)
}


pseudotime_inference <- function(input_file = "sceset.rds",
                                 output_file = "myfile.csv") {
  
  sce <- readRDS(input_file)
  exprs_mat <- exprs(sce)
  
  non_zero_counts <- rowMeans(exprs_mat) > 0
  exprs_mat <- exprs_mat[non_zero_counts,]
  
  fit <- fit_phenopath(exprs_mat, sce$x)
  pst <- trajectory(fit)

  true_pseudotimes <- sce$pst
  abs_kendall <- abs(cor(pst, true_pseudotimes, method = "kendall", use = "na"))  
  
  original_prop_zero <- sce$original_prop_zero[1]
  new_prop_zero <- sce$new_prop_zero[1]
  proportion_set_to_zero <- sce$proportion_set_to_zero[1]
  
  output_df <- data.frame(input_sceset = input_file,
                          kendall_correlation = abs_kendall,
                          original_prop_zero,
                          new_prop_zero,
                          proportion_set_to_zero)
  
  write.csv(output_df, output_file)
  

}

aargh(pseudotime_inference)
