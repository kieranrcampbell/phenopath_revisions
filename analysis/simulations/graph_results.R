library(readr)
library(ggplot2)
library(cowplot)
library(dplyr)

theme_set(theme_cowplot(font_size = 11))

allcor_file <- "data/simulations/all_pseudotime_correlations.csv"

df_split <- read_csv(allcor_file)

ggplot(df_split, aes(x = as.factor(100 * p), 
                     y = kendall_correlation, fill = alg)) +
  geom_boxplot() + # geom_line() +
  labs(y = "Absolute Kendall correlation\nto true pseudotime",
       x = "% genes covariate interaction",
       subtitle = "Correlation to true pseudotime across 40 replications for each condition") 


dfg <- group_by(df_split, N, G, p, alg, noise) %>% 
  summarise(mean_correlation = median(kendall_correlation),
            lower = quantile(kendall_correlation, 0.25),
            upper = quantile(kendall_correlation, 0.75))

alg_from <- c("dpt", "monocle2", "phenopath", "tscan")
alg_to <- c("DPT", "Monocle 2", "PhenoPath", "TSCAN")

dfg$alg <- plyr::mapvalues(dfg$alg, from = alg_from, to = alg_to)

dfg$N_str <- paste0(dfg$N, " cells")

ggplot(dfg, aes(x = as.factor(100 * p), y = mean_correlation,
                color = alg, group = alg)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  geom_line(size = 1.5) + 
  geom_point(aes(fill = alg), shape = 21, color = 'black', size = 2) +
  scale_colour_brewer(palette = "Set2", name = "Algorithm") +
  scale_fill_brewer(palette = "Set2", name = "Algorithm") +
  labs(y = "Median absolute Kendall-tau\ncorrelation to true pseudotime",
                     x = "% genes covariate interaction") + #,
                     # subtitle = "Correlation to true pseudotime across 40 replications for each condition")+
  ylim(c(-0.01, 1)) + 
  facet_grid(noise ~ N_str) +
  theme(strip.background = element_rect(fill = "white"))

ggsave("figs/pseudotime_correlation.png", width = 8, height = 3)

ggplot(filter(dfg, N == 200, noise == "low"), aes(x = as.factor(100 * p), y = mean_correlation,
                color = alg, group = alg)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  geom_line(size = 1.5) + 
  geom_point(aes(fill = alg), shape = 21, color = 'black', size = 2) +
  scale_colour_brewer(palette = "Set2", name = "Algorithm") +
  scale_fill_brewer(palette = "Set2", name = "Algorithm") +
  labs(y = "Median absolute Kendall-tau\ncorrelation to true pseudotime",
       x = "% genes covariate interaction") + #,
  # subtitle = "Correlation to true pseudotime across 40 replications for each condition")+
  ylim(c(-0.01, 1)) + 
  # facet_wrap(~ N_str) +
  theme(strip.background = element_rect(fill = "white"))

pst_plot <- last_plot()


# AUC ---------------------------------------------------------------------

aucs <- read_csv("data/simulations/roc.csv")
aucs <- mutate(aucs, DE_alg = "limma")

aucs_phenopath <- read_csv("data/simulations/roc_phenopath.csv")
aucs_phenopath <- mutate(aucs_phenopath, DE_alg = "phenopath")

aucs_deseq <- read_csv("data/simulations/roc_deseq.csv")
aucs_deseq <- mutate(aucs_deseq, DE_alg = "deseq2")

aucs_mast <- read_csv("data/simulations/roc_mast.csv")
aucs_mast <- mutate(aucs_mast, DE_alg = "mast")

aucs_monocle <- read_csv("data/simulations/roc_monocle.csv")
aucs_monocle <- mutate(aucs_monocle, DE_alg = "monocle")


aucs <- bind_rows(aucs, aucs_phenopath, aucs_deseq, aucs_mast, aucs_monocle)#, aucs_randmonocle)

aucs <- dplyr::filter(aucs, auc != -1) # This is when the fit didn't work

dfg <- group_by(aucs, N, G, p, alg, DE_alg, noise) %>% 
  summarise(mean_auc = median(auc),
            lower = quantile(auc, 0.25),
            upper = quantile(auc, 0.75))

alg_from <- c("dpt", "monocle2", "phenopath", "tscan")
alg_to <- c("DPT", "Monocle 2", "PhenoPath", "TSCAN")

dfg$alg <- plyr::mapvalues(dfg$alg, from = alg_from, to = alg_to)

dfg$N_str <- paste0(dfg$N, " cells")

lower <- min(dfg$mean_auc) - 0.01
upper <- 1.05

graph_auc <- function(dfg, de_alg, de_str, noi) {
  filter(dfg, DE_alg %in% c(de_alg, "phenopath"), noise == noi) %>% 
    ggplot(aes(x = as.factor(100 * p), y = mean_auc, color = alg, group = alg)) + 
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
      geom_line(size = 1.5) + 
      geom_point(aes(fill = alg), shape = 21, color = 'black', size = 2) +
      scale_colour_brewer(palette = "Set2", name = "Algorithm") +
      scale_fill_brewer(palette = "Set2", name = "Algorithm") +
      labs(y = "Median AUC",
           x = "% genes covariate interaction",
           subtitle = paste("DE with", de_str)) +
      facet_wrap(~ N_str) +
      ylim(lower, upper)
}


auc_plots <- lapply(c("low", "high"), function(noise) {
  voom_plot <- graph_auc(dfg, "limma", "Limma Voom", noise)
  deseq_plot <- graph_auc(dfg, "deseq2", "DESeq2", noise)
  mast_plot <- graph_auc(dfg, "mast", "MAST", noise)
  monocle_plot <- graph_auc(dfg, "monocle", "Monocle 2", noise)
  plot_grid(voom_plot, deseq_plot, mast_plot, monocle_plot, ncol = 1)
})

names(auc_plots) <- c('low', 'high')


ggsave("figs/supp_auc_low.png", auc_plots[['low']], width = 8, height = 11)
ggsave("figs/supp_auc_high.png", auc_plots[['high']], width = 8, height = 11)




voom_paper_plot <- 
  filter(dfg, DE_alg %in% c("limma", "phenopath"), noise == "low", N == 200) %>% 
  ggplot(aes(x = as.factor(100 * p), y = mean_auc,
                                            color = alg, group = alg)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  geom_line(size = 1.5) + 
  geom_point(aes(fill = alg), shape = 21, color = 'black', size = 2) +
  scale_colour_brewer(palette = "Set2", name = "Algorithm") +
  scale_fill_brewer(palette = "Set2", name = "Algorithm") +
  labs(y = "Median AUC",
       x = "% genes covariate interaction") +
  # facet_wrap(~ N_str) +
  # ylim(lower, upper) +
  theme(legend.position = "none")


# Paper figure ------------------------------------------------------------

sim_example <- readRDS("figs/simulation_example.rds")

bottom_grid <- plot_grid(pst_plot, voom_paper_plot, rel_widths = c(4,3),
                         labels = c("C", "D"), scale = 0.95)

plot_grid(plot_grid(sim_example, scale = 0.95), NULL, bottom_grid, ncol = 1, labels = c("A", "B", NULL),
          rel_heights = c(1.2,2, 2))

ggsave("figs/sim.png", width = 10, height = 7)
