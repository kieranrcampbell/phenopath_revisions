library(scater)

scale_vec <- function(x) (x - mean(x)) / sd(x)

sce <- readRDS("data/paper-scesets/sce_shalek.rds")

sce <- sce[, sce$stimulant != "PIC"]
sce$stimulant <- droplevels(sce$stimulant)


plotPhenoData(sce, aes(x = total_counts, y = total_features))

# Remove cells with < 1e4 counts

sce <- sce[, sce$total_counts > 1e4]



Lyz1 <- exprs(sce)["LYZ1",]
Serpinb6b <- exprs(sce)["SERPINB6B",]

Serpinb6b_threshold <- 4
Lyz1_threshold <- 5

to_keep <- Lyz1 > Lyz1_threshold & Serpinb6b < Serpinb6b_threshold

qplot(Lyz1, Serpinb6b, color = to_keep) +
  geom_vline(xintercept = Lyz1_threshold, linetype = 2) +
  geom_hline(yintercept = Serpinb6b_threshold, linetype = 2)

sce <- sce[, to_keep]

m <- model.matrix(~ sce$pct_dropout)
sc <- normaliseExprs(sce, design = m)
exprs(sc) <- norm_exprs(sc)

sc$x <- scale_vec( sc$stimulant == "LPS" )

saveRDS(sc, "data/paper-scesets/sce_shalek_qc.rds")
