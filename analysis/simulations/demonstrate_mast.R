library(scater)
library(ggplot2)
library(AUC)
library(dplyr)
library(cowplot)
library(MAST)

set.seed(246L)

theme_set(theme_cowplot(font_size = 11))

sce <- readRDS("data/simulations/scesets/sceset_N_200_G_500_p_0.5_rep_8.rds")

N <- ncol(sce)

pseudotime <- rnorm(N)

make_pst_plot <- function(pseudotime, labels = c("A", "B", "C")) {
  pst_plot <- qplot(factor(sce$x), pseudotime, geom = 'boxplot') +
    labs(x = "x", y = "Pseudotime")
  
  
  pdata <- data.frame(x = sce$x, pseudotime)
  rownames(pdata) <- colnames(sce)
  
  ## Remove ----
  sce$pseudotime <- pseudotime
  lm_pvals <- apply(exprs(sce), 1, function(y) {
    fit <- lm(y ~ sce$pseudotime * sce$x)
    coef(summary(fit))[4,'Pr(>|t|)']
  })
  
  ## End remove ---
  
  sca <- FromMatrix(exprs(sce), pdata)
  fit <- zlm(~ x + pseudotime + x:pseudotime, sca)
  
  lrt <- lrTest(fit, "x:pseudotime")
  
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
  
  plot_grid(pst_plot, pval_plot, auc_plot, nrow = 1, labels = labels)
}

top_grid <- make_pst_plot(rnorm(N))
bottom_grid <- make_pst_plot(rnorm(N, sce$x, 0.5), labels = c("D", "E", "F"))

plot_grid(top_grid, bottom_grid, ncol = 1)

# ggsave("figs/mast.png", width = 7, height = 5)

