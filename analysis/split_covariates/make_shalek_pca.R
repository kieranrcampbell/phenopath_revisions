library(scater)
library(cowplot)

sce <- readRDS("data/paper-scesets/sce_shalek_qc.rds")

pca_time <- plotPCA(sce, colour_by = "time", ncomponents = 3)
pca_stimulant <- plotPCA(sce, colour_by = "stimulant", ncomponents = 3)

plot_grid(pca_time, pca_stimulant, ncol = 1, labels = "AUTO")

ggsave("figs/supp_shalek_pca.png", width = 5, height = 8)