library(ggplot2)
library(readr)

aucs <- read_csv("data/simulations/roc.csv")

ggplot(aucs, aes(x = as.factor(100 * p), y = auc, fill = alg)) +
  geom_boxplot()



dfg <- group_by(aucs, N, G, p, alg) %>% 
  summarise(mean_auc = mean(auc),
            lower = quantile(auc, 0.25),
            upper = quantile(auc, 0.75))

alg_from <- c("dpt", "monocle2", "phenopath", "tscan")
alg_to <- c("DPT", "Monocle 2", "PhenoPath", "TSCAN")

dfg$alg <- plyr::mapvalues(dfg$alg, from = alg_from, to = alg_to)

ggplot(dfg, aes(x = as.factor(100 * p), y = mean_auc,
                color = alg, group = alg)) + 
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  geom_line(size = 1.5) + 
  geom_point(aes(fill = alg), shape = 21, color = 'black', size = 2) +
  scale_colour_brewer(palette = "Set2", name = "Algorithm") +
  scale_fill_brewer(palette = "Set2", name = "Algorithm") +
  labs(y = "Mean AUC",
       x = "% genes covariate interaction")+
  ylim(c(-0.01, 1))

