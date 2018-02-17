library(scater)
library(dplyr)
library(magrittr)
library(biomaRt)

sce <- readRDS("data/paper-scesets/sce_shalek_clvm.rds")

sce <- plotPCA(sce, colour_by = "time", ntop = 1e3, return_SCESet = TRUE)

cluster <- 1 * (redDim(sce)[,1] < 0 & redDim(sce)[,2] < 0)

sce$cluster <- cluster

plotPCA(sce, colour_by = "cluster", ntop = 1e3)

lm_fits <- apply(exprs(sce), 1, function(y) {
  fit <- lm(y ~ sce$cluster)
  s <- summary(fit)
  coef(s)[2, c("Estimate", "Pr(>|t|)")]
})

df <- as_data_frame(t(lm_fits))
names(df) <- c("estimate", "p_val")

df <- mutate(df, q_val = p.adjust(p_val, method = "BH"),
             gene = stringr::str_to_title(featureNames(sce)))

ggplot(df, aes(x = estimate, y = -log10(q_val))) +
  geom_point() +
  labs(x = "Fold change", y = "-log10(q-value)")

arrange(df, q_val)

plotPCA(sce, colour_by = "LYZ1", ntop = 1e3)

## Need to add ensembl mus id
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
bm <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
            filters = "mgi_symbol",
            values = stringr::str_to_title(df$gene),
            mart = mart)

bm$gene <- bm$mgi_symbol

df <- inner_join(df, bm, by = "gene")

library(goseq)

#' Run all steps for a goseq analysis
#'
#' @param genes The genes "differentially expressed" or in the condition
#' @param all_genes The background gene set
#' @param genome The genome of genes and all_genes, e.g. mm10 or hg19
#' @param test_cats The test categories, any combination of GO:BP, GO:CC, or GO:MF
#' @return The output of a call to \code{Goseq}
go_analysis <- function(genes, all_genes, genome = "mm10", test_cats = "GO:BP") {
  lgl <- 1 * (all_genes %in% genes)
  names(lgl) <- all_genes
  pwf <- nullp(lgl, genome, "ensGene")
  go <- goseq(pwf, genome, "ensGene", test.cats = test_cats)
  return(go)
}

genes_up <- filter(df, estimate > 0, q_val < 0.05) %>% .$ensembl_gene_id
genes_down <- filter(df, estimate < 0, q_val < 0.05) %>% .$ensembl_gene_id

go_up <- go_analysis(genes_up, df$ensembl_gene_id)
go_down <- go_analysis(genes_down, df$ensembl_gene_id)


library(gplots)
library(viridis)
library(RColorBrewer)

top_genes <- filter(df, q_val < 0.05) %>% 
  arrange(q_val) %>% 
  head(n = 80) %>% 
  .$gene


small_exprs <- exprs(sce[stringr::str_to_upper(top_genes), ])

csc <- brewer.pal(3, "Set1")[sce$cluster + 1]

png("~/Desktop/heatmap.png", width = 1000, height = 800)
gplots::heatmap.2(small_exprs, trace = "none",
                  scale = "none", col = "viridis",
                  ColSideColors = csc)
dev.off()

