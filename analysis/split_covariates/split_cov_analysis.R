library(scater)
library(cowplot)
library(monocle)
library(phenopath)
library(tidyverse)
library(magrittr)

theme_set(theme_cowplot(font_size = 11))

fit_monocle2 <- function(sce) {
  exprs_mat <- exprs(sce)
  cds <- newCellDataSet(exprs_mat)
  sizeFactors(cds) <- rep(1, ncol(cds))
  cds <- setOrderingFilter(cds, rownames(exprs_mat))
  cds <- reduceDimension(cds, norm_method = "none")
  cds <- orderCells(cds)
  cds
}

get_R2 <- function(pseudotime, time) {
  time <- as.numeric(gsub("h", "", time))
  fit <- lm(pseudotime ~ time)
  s <- summary(fit)
  s$r.squared
}

sce <- readRDS("data/paper-scesets/sce_shalek_qc.rds")

pca_time <- plotPCA(sce, colour_by = "time", ncomponents = 3)
pca_stimulant <- plotPCA(sce, colour_by = "stimulant", ncomponents = 3)

plot_grid(pca_time, pca_stimulant, ncol = 1, labels = "AUTO")

ggsave("figs/supp_shalek_pca.png", width = 5, height = 8)

row_vars <- matrixStats::rowVars(exprs(sce))
high_var <- row_vars >= sort(row_vars, decreasing = TRUE)[5000]

sce_hvg <- sce[high_var, ]

sce_lps <- sce_hvg[, sce$stimulant == "LPS"]
sce_pam <- sce_hvg[, sce$stimulant == "PAM"]

sces <- list(
  lps = sce_lps,
  pam = sce_pam
)

cds_sets <- lapply(sces, fit_monocle2)
pseudotimes <- lapply(cds_sets, function(cds) cds$Pseudotime)

get_R2(pseudotimes[[1]], sces[[1]]$time)
get_R2(pseudotimes[[2]], sces[[2]]$time)

pp_fit <- phenopath(sce_hvg, sce_hvg$x)

get_R2(trajectory(pp_fit), sce_hvg$time)

# cds_sets_counts <- lapply(names(cds_sets), function(stimulant) {
#   cds <- cds_sets[[stimulant]]
#   exprs(cds) <- counts(sces[[stimulant]])
#   cds
# })



de <- lapply(cds_sets_counts, function(cds) {
  differentialGeneTest(cds, fullModelFormulaStr = "~ sm.ns(Pseudotime, df = 3)")
})

log10_qvals <- -log10(de[[2]]$qval)

log10_qvals[is.infinite(log10_qvals)] <- max(log10_qvals[is.finite(log10_qvals)]) + 50
plot(log10_qvals)

de_linear <- lapply(cds_sets_counts, function(cds) {
  cds$random <- runif(ncol(cds))
  differentialGeneTest(cds, fullModelFormulaStr = "~ random",
                       reducedModelFormulaStr = "~ 1")
})


y <- exprs(cds_sets[[1]])[1,]
pseudotime <- pseudotimes[[1]]

de_cds <- function(cds) {
  apply(exprs(cds), 1, function(y) {
    fit <- vglm(y ~ sm.ns(cds$Pseudotime, df = 3), family = gaussianff())
    fit_null <- vglm(y ~ 1, family = gaussianff())
    lrt <- lrtest(fit, fit_null)
    p_value <- lrt@Body$`Pr(>Chisq)`[2]
    p_value
  })
}

for(i in 1:2) {
  cds_sets[[i]]$time <- sces[[i]]$time
  cds <- cds_sets[[i]]
  time_numeric <- as.numeric(gsub("h", "", cds$time))
  if(cor(cds$Pseudotime, time_numeric) < 0) {
    cds_sets[[i]]$Pseudotime <- max(cds$Pseudotime) - cds$Pseudotime
  }
}

de_cds_2 <- function(cds) {
  apply(exprs(cds), 1, function(y) {
    fit <- lm(y ~ cds$Pseudotime)
    s <- summary(fit)
    p_value <- coef(s)[2, c("Estimate", "Pr(>|t|)")]
    t(p_value)
  })
}

des <- lapply(cds_sets, de_cds_2)

beta_diff <- des[[1]][1,] - des[[2]][1,]
pvals <- lapply(des, function(d) d[2,])
qvals <- lapply(pvals, p.adjust, method = "BH")
log10qvals <- lapply(qvals, function(q) -log10(q))

qval_diff <- log10qvals[[1]] - log10qvals[[2]]

plot(beta_diff, qval_diff)


sig <- significant_interactions(pp_fit, 3)
qplot(beta_diff, pp_fit$m_beta[1,]) +
  labs(x = expression(beta[LPS]-beta[PAM]~"from linear model"),
       y = expression("PhenoPath"~beta))

ggsave("figs/shalek-linear-model-fit.png", width = 4, height = 3)

which_min <- which.min(beta_diff)

plot(cds_sets[[1]]$Pseudotime, exprs(cds_sets[[1]])[which_min,])
plot(cds_sets[[2]]$Pseudotime, exprs(cds_sets[[2]])[which_min,])

inds <- which_min # match(genes, fData(sce_hvg)$mgi_symbol)

df <- t(exprs(sce_hvg))[, inds, drop=FALSE] 
colnames(df) <- fData(sce_hvg)$mgi_symbol[which_min]

df <- as_data_frame(df) %>% 
  dplyr::mutate(capture_time = sce_hvg$time, 
                x = as.character(sce_hvg$stimulant),
                phenopath = trajectory(pp_fit)) %>% 
  gather(gene, expression, -capture_time, -x, -phenopath)

# df$gene <- factor(df$gene, levels = genes)

df_g <- group_by(df, capture_time, x, gene) %>% 
  summarise(median_expression = median(expression))

ggplot(df, aes(x = capture_time, y = expression, fill = x, color = x)) + 
  geom_jitter(shape = 21, color = 'grey30', alpha = 0.3, width = 0.1) +
  facet_grid(gene ~ x, scales = "free") + 
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "none", strip.text = element_text(size = 10, face = "bold"),
        strip.text.y = element_text(size = 8),
        strip.background = element_rect(colour = "black", fill = "grey95", 
                                        size = 0.5, linetype = 1)) +
  theme(axis.text = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ylab(expression("Normalised expression")) + 
  xlab("Capture time") +
  geom_point(data = df_g, aes(y = median_expression), 
             shape = 21, color = 'black', size = 2) +
  geom_line(data = df_g, aes(y = median_expression, group = 1))

ggplot(df, aes(x = phenopath, y = expression, color = x)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ x)


tnf <- grep("^Tnf$", fData(sce_hvg)$mgi_symbol)
plot(cds_sets[[1]]$Pseudotime, exprs(cds_sets[[1]])[tnf,])
plot(cds_sets[[2]]$Pseudotime, exprs(cds_sets[[2]])[tnf,])



