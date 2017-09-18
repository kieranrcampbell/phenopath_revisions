library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(RColorBrewer)

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
alg_to <- c("DPT", "Monocle 2", "TSCAN")

dataset_from <- c("chu", "dulken", "hsc", "trapnell")
dataset_to <- c("Chu et al.", "Dulken et al.", "Zhou et al.", "Trapnell et al.")

dfg$algorithm <- plyr::mapvalues(dfg$algorithm, from = alg_from, to = alg_to)
dfg$dataset <- plyr::mapvalues(dfg$dataset, from = dataset_from, to = dataset_to)

cols <- brewer.pal(4, "Set2")[c(1,2,4)]

ggplot(dfg, aes(x = dataset, y = correlation, fill = algorithm)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  labs(y = "Correlation to PC1") +
  # scale_fill_manual(values = cols, name = "Algorithm") +
  scale_fill_brewer(palette = "Accent", name = "Algorithm") +
  xlab("Dataset") +
  theme(legend.position = "top",
        axis.text = element_text(size = 10))

ggsave("figs/vs-pc1.png", width = 4.5, height = 3)
