
library(MultiAssayExperiment)
library(scater)
library(biomaRt)
library(cowplot)

theme_set(theme_cowplot(font_size = 11))

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

plotPhenoData(sce, aes(x = total_counts, y = total_features, color = to_keep)) +
  labs(subtitle = "QC pass: total_features > 5000 and total_counts > 50000")
ggsave("figs/supp_shalek_qc.png", width = 5, height = 3)

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
  geom_hline(yintercept = Serpinb6b_threshold, linetype = 2) +
  scale_color_brewer(palette = "Dark2") +
  labs(subtitle = "Non-cluster-disrupted: Serpinb6b > 2.5 and Lyz1 > 0")

ggsave("figs/supp_shalek_qc_2.png", width = 5, height = 3)

sce_qc2 <- sce_qc[, to_keep]

sce_qc2$x <- 1 * (sce_qc2$stimulant == "LPS")

qplot(sce_qc2$total_features) + xlab("Total features")
hist_plot <- last_plot()
plotQC(sce_qc2, type = 'find', var = 'total_features', ntop = 2e3)
pc_plot <- last_plot()

plot_grid(hist_plot, pc_plot + xlab("Total features"), labels = "AUTO")
ggsave("figs/supp_shalek_qc_3.png", width = 7, height = 3)

m <- model.matrix(~ sce_qc2$total_features)
sce_qc2 <- normaliseExprs(sce_qc2, design = m)
exprs(sce_qc2) <- norm_exprs(sce_qc2)

saveRDS(sce_qc2, "data/paper-scesets/sce_shalek_qc.rds")

