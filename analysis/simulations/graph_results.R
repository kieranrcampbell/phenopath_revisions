library(readr)
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot(font_size = 11))

allcor_file <- "data/simulations/all_pseudotime_correlations.csv"

df_split <- read_csv(allcor_file)

ggplot(df_split, aes(x = as.factor(100 * p), 
                     y = kendall_correlation, fill = alg)) +
  geom_boxplot() + # geom_line() +
  labs(y = "Absolute Kendall correlation\nto true pseudotime",
       x = "% genes covariate interaction",
       subtitle = "Correlation to true pseudotime across 40 replications for each condition") 


dfg <- group_by(df_split, N, G, p, alg) %>% 
  summarise(mean_correlation = mean(kendall_correlation),
            lower = quantile(kendall_correlation, 0.25),
            upper = quantile(kendall_correlation, 0.75))

alg_from <- c("dpt", "monocle2", "phenopath", "tscan")
alg_to <- c("DPT", "Monocle 2", "PhenoPath", "TSCAN")

dfg$alg <- plyr::mapvalues(dfg$alg, from = alg_from, to = alg_to)

dfg$N_str <- paste0(dfg$N, " cells")

ggplot(dfg, aes(x = as.factor(100 * p), y = mean_correlation,
                color = alg, group = alg)) + 
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  geom_line(size = 1.5) + 
  geom_point(aes(fill = alg), shape = 21, color = 'black', size = 2) +
  scale_colour_brewer(palette = "Set2", name = "Algorithm") +
  scale_fill_brewer(palette = "Set2", name = "Algorithm") +
  labs(y = "Absolute Kendall-tau\ncorrelation to true pseudotime",
                     x = "% genes covariate interaction",
                     subtitle = "Correlation to true pseudotime across 40\nreplications for each condition")+
  ylim(c(-0.01, 1)) + 
  facet_wrap(~ N_str)

ggsave("figs/pseudotime_correlation.png", width = 8, height = 3)


# AUC ---------------------------------------------------------------------

aucs <- read_csv("data/simulations/roc.csv")
aucs_phenopath <- read_csv("data/simulations/roc_phenopath.csv")
aucs <- bind_rows(aucs, aucs_phenopath)


dfg <- group_by(aucs, N, G, p, alg) %>% 
  summarise(mean_auc = mean(auc),
            lower = quantile(auc, 0.25),
            upper = quantile(auc, 0.75))

alg_from <- c("dpt", "monocle2", "phenopath", "tscan")
alg_to <- c("DPT", "Monocle 2", "PhenoPath", "TSCAN")

dfg$alg <- plyr::mapvalues(dfg$alg, from = alg_from, to = alg_to)

dfg$N_str <- paste0(dfg$N, " cells")

ggplot(dfg, aes(x = as.factor(100 * p), y = mean_auc,
                color = alg, group = alg)) + 
  # geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  geom_line(size = 1.5) + 
  geom_point(aes(fill = alg), shape = 21, color = 'black', size = 2) +
  scale_colour_brewer(palette = "Set2", name = "Algorithm") +
  scale_fill_brewer(palette = "Set2", name = "Algorithm") +
  labs(y = "Mean AUC",
       x = "% genes covariate interaction") +
  facet_wrap(~ N_str)

ggsave("figs/auc.png", width = 8, height = 3)
