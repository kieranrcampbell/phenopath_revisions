library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(RColorBrewer)

theme_set(theme_cowplot(font_size = 11))

files <- dir("data/linpst", full.names = TRUE)

dfs <- lapply(files, read_csv)
df <- bind_rows(dfs)

dfg <- group_by(df, dataset) %>% 
  summarise(dpt = cor(dpt, pca),
            monocle2 = cor(monocle2, pca),
            tscan = cor(tscan, pca, use = "na"))
dfg <- gather(dfg, algorithm, correlation, -dataset) %>% 
  mutate(correlation = abs(correlation))

alg_from <- c("dpt", "monocle2", "tscan")
alg_to <- c("DPT ", "Monocle 2 ", "TSCAN ")

dataset_from <- c("chu", "dulken", "hsc", "trapnell")
dataset_to <- c("Chu et al.", "Dulken et al.", "Zhou et al.", "Trapnell et al.")

dfg$algorithm <- plyr::mapvalues(dfg$algorithm, from = alg_from, to = alg_to)
dfg$dataset <- plyr::mapvalues(dfg$dataset, from = dataset_from, to = dataset_to)

cols <- brewer.pal(4, "Set2")[c(1,2,4)]

ggplot(dfg, aes(x = dataset, y = correlation, fill = algorithm)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  labs(y = "Correlation\nto PC1") +
  # scale_fill_manual(values = cols, name = "Algorithm") +
  scale_fill_brewer(palette = "Accent", name = "Algorithm") +
  xlab("Dataset") +
  theme(legend.position = "none",
        axis.text = element_text(size = 10))

cor_to_pc1 <- last_plot()

# ggsave("figs/vs-pc1.png", width = 4.5, height = 3)


# Proportion genes linear -------------------------------------------------


files <- dir("data/lincoef", full.names = TRUE)

dfs <- lapply(files, read_csv)
df <- bind_rows(dfs)
dfg <- gather(df, algorithm, qval, -dataset)
dfg <- group_by(dfg, dataset, algorithm) %>% 
  summarise(mean_signif = mean(qval < 0.05))

dfg$algorithm <- plyr::mapvalues(dfg$algorithm, from = alg_from, to = alg_to)
dfg$dataset <- plyr::mapvalues(dfg$dataset, from = dataset_from, to = dataset_to)

ggplot(dfg, aes(x = dataset, y = mean_signif, fill = algorithm)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  labs(y = "Proportion of genes\nsignificant linear trend") +
  # scale_fill_manual(values = cols, name = "Algorithm") +
  scale_fill_brewer(palette = "Accent", name = "Algorithm") +
  xlab("Dataset") +
  theme(legend.position = "top",
        axis.text = element_text(size = 10))

prop_lin_trend <- last_plot()

plot_grid(cor_to_pc1, prop_lin_trend, nrow = 2, rel_heights = c(3,4),
          labels = "AUTO")

ggsave("figs/linear.png", width = 5, height = 5)
