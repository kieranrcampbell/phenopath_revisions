
library(MultiAssayExperiment)
library(scater)
library(biomaRt)

scale_vec <- function(x) (x - mean(x)) / sd(x)


mae <- readRDS("data/paper-scesets/GSE48968-GPL13112.rds")

cts <- assays(experiments(mae)[["gene"]])[["count_lstpm"]]
tpms <- assays(experiments(mae)[["gene"]])[["TPM"]]
phn <- colData(mae)

sce <- newSCESet(countData = cts, 
                  phenoData = new("AnnotatedDataFrame", data = as.data.frame(phn)))
tpm(sce) <- tpms
exprs(sce) <- log2(tpm(sce) + 1)

is_lps_pam <- grepl("LPS|PAM", sce$description)

sce <- sce[, is_lps_pam]

split <- strsplit(as.character(sce$description), "_", fixed = TRUE)
stimulant <- sapply(split, `[`, 1)
time <- sapply(split, `[`, 2)

sce$stimulant <- stimulant
sce$time <- time

sce <- calculateQCMetrics(sce)

plotPhenoData(sce, aes(x = total_counts, y = total_features))
sce$to_keep <- sce$total_counts > 5e5 & sce$total_features > 5e3
plotPhenoData(sce, aes(x = total_counts, y = total_features, color = to_keep))
sce_qc <- sce[, sce$to_keep]

ensembl_gene_ids <- sapply(strsplit(featureNames(sce_qc), ".", fixed = TRUE), `[`, 1)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
bm <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
            filters = "ensembl_gene_id",
            values = ensembl_gene_ids,
            mart = mart)

fData(sce_qc)$mgi_symbol <- rep(NA, nrow(sce_qc))

mm2 <- match(bm$ensembl_gene_id, ensembl_gene_ids)
fData(sce_qc)$mgi_symbol[mm2] <- bm$mgi_symbol

Lyz1_index <- grep("Lyz1", fData(sce_qc)$mgi_symbol)
SerpinB6b_index <- grep("SerpinB6b", fData(sce_qc)$mgi_symbol, ignore.case = TRUE)

Lyz1 <- exprs(sce_qc)[Lyz1_index,]
Serpinb6b <- exprs(sce_qc)[SerpinB6b_index,]

Serpinb6b_threshold <- 2.5
Lyz1_threshold <- 0

to_keep <- Lyz1 > Lyz1_threshold & Serpinb6b < Serpinb6b_threshold

qplot(Lyz1, Serpinb6b, color = to_keep) +
  geom_vline(xintercept = Lyz1_threshold, linetype = 2) +
  geom_hline(yintercept = Serpinb6b_threshold, linetype = 2)

sce_qc2 <- sce_qc[, to_keep]

sce_qc2$x <- sce_qc2$stimulant == "LPS"

saveRDS(sce_qc2, "data/paper-scesets/sce_shalek_qc.rds")

