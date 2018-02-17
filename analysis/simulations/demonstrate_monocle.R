# library(scater)
# library(ggplot2)
library(AUC)
library(dplyr)
# library(cowplot)
library(monocle)
library(readr)
library(ggplot2)

set.seed(246L)

theme_set(theme_cowplot(font_size = 11))

sce <- readRDS("data/simulations/scesets/sceset_N_200_G_500_p_0.05_rep_1.rds")
pseudotime_df <- read_csv("data/simulations/pseudotimes/pseudofit_N_200_G_500_p_0.05_rep_1_alg_monocle2.csv")
pseudotime <- scale(pseudotime_df$pst)[,1]


N <- ncol(sce)



sce$pseudotime <- scale(pseudotime)[,1]

cds <- newCellDataSet(counts(sce), new("AnnotatedDataFrame", pData(sce)))



cds$covariate <- cds$x

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

monocle_auc <- function(z) {
  cds$z <- z
  de <- differentialGeneTest(cds, 
                                  fullModelFormulaStr = "~ covariate*sm.ns(z, df = 3)",
                                  reducedModelFormulaStr = "~ covariate + sm.ns(z, df = 3)")
  
  plot(-log10(de$pval))
  
  roc_obj <- roc(1 - de$qval, factor(fData(sce)$is_interaction))
  
  a <- auc(roc_obj)
  
  auc_text <- frame_data(
    ~ x, ~ y, ~ auc,
    0.25, 0.75, paste0("AUC = ", signif(a, 2))
  )
  
  qplot(roc_obj$fpr, roc_obj$tpr, geom = 'line') +
    stat_function(fun = function(x) x, linetype = 2) +
    geom_text(data = auc_text, aes(x = x, y = y, label = auc)) +
    labs(x = "FPR", y = "TPR")
}

z <- rnorm(N)
monocle_auc(sce$pst)
monocle_auc(z)
monocle_auc(z - min(z))

z2 <- rnorm(N, sce$x, 0.1)
monocle_auc(z2)


a <- auc(roc_obj)


de_true <- differentialGeneTest(cds, 
                           fullModelFormulaStr = "~ x*sm.ns(pst, df = 3)",
                           reducedModelFormulaStr = "~ x + sm.ns(pst, df = 3)")

plot(-log10(de_true$pval))


de_ns <- differentialGeneTest(cds, 
                              fullModelFormulaStr = "~ Branch*sm.ns(Pseudotime, df = 3)",
                              reducedModelFormulaStr = "~ Branch + sm.ns(Pseudotime, df = 3)")

plot(-log10(de_ns$pval))

cds$absolute_crap <- cds$x

de_ns <- differentialGeneTest(cds, 
                              fullModelFormulaStr = "~ absolute_crap*sm.ns(pst, df = 3)",
                              reducedModelFormulaStr = "~ absolute_crap + sm.ns(pst, df = 3)")

plot(-log10(de_ns$pval))

de_inferred <- differentialGeneTest(cds, 
                              fullModelFormulaStr = "~ Branch*sm.ns(pseudotime, df = 3)",
                              reducedModelFormulaStr = "~ Branch + sm.ns(pseudotime, df = 3)")

plot(-log10(de_inferred$pval))

roc_obj <- roc(1 - de_inferred$qval, factor(fData(sce)$is_interaction))

a <- auc(roc_obj)

de_true_lin <- differentialGeneTest(cds, 
                                fullModelFormulaStr = "~ x*pst",
                                reducedModelFormulaStr = "~ x + pst")

plot(-log10(de_true_lin$pval))

pvals <- apply(exprs(cds), 1, function(y) {
  fit <- lm(y ~ cds$pst * cds$x)
  coef(summary(fit))[4,'Pr(>|t|)']
})
plot(-log10(pvals))


pval_plot <- qplot(fData(sce)$is_interaction, -log10(lrt[,3,3]), geom = "violin") +
  labs(x = "Interaction simulated",
       y = "-log10(p-value)")

roc_obj <- roc(1 - lrt[,3,3], factor(fData(sce)$is_interaction))

a <- auc(roc_obj)

auc_text <- frame_data(
  ~ x, ~ y, ~ auc,
  0.25, 0.75, paste0("AUC = ", signif(a, 2))
)

auc_plot <- qplot(roc_obj$fpr, roc_obj$tpr, geom = 'line') +
  stat_function(fun = function(x) x, linetype = 2) +
  geom_text(data = auc_text, aes(x = x, y = y, label = auc)) +
  labs(x = "FPR", y = "TPR")

plot_grid(pst_plot, pval_plot, auc_plot, nrow = 1, labels = 'AUTO')

ggsave("figs/mast.png", width = 7, height = 2.5)

