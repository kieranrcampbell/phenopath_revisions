
# Run CAVI for CLVM and save result

# Format:
# Rscript clvm [input sceset] [output file] [pc to initialise] [optional covariate name]

library(phenopath)
library(scater)

args <- commandArgs(trailingOnly = TRUE)

input_rds <- args[1]
output_file <- args[2]

pc_initialise <- as.numeric(args[3])

elbo_tol <- as.numeric(args[4])

cov_name <- NULL

if(length(args) > 4) {
  cov_name <- args[5]
} else {
  cov_name <- 'x'
}

sce <- readRDS(input_rds)
sce <- updateSCESet(sce)

y <- t(assay(sce, "logcounts"))
x <- cbind(colData(sce)[[ cov_name ]])

G <- nrow(sce)
true_tau <- rep(1, G)

pcavi <- phenopath(y, x, z_init = 1, maxiter = 800, 
              elbo_tol = elbo_tol, thin = 4)

saveRDS(pcavi, file = output_file)



