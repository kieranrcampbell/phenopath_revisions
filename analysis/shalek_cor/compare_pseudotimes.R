

library(tidyverse)
library(cowplot)

theme_set(theme_cowplot(font_size = 11))

hvg_dir <- "data/hvg"

all_files <- dir(hvg_dir, full.names = TRUE)

df_list <- lapply(all_files, read_csv)

df <- bind_rows(df_list)

df_ref <- filter(df, hvg == 4000) %>% select(-hvg)
df_comp <- filter(df, hvg != 4000)

df_all <- inner_join(df_ref, df_comp,
                     by = c("sample", "algorithm", "dataset"),
                    suffix = c("_all", "_4000"))

df_cor <- group_by(df_all, algorithm, dataset, hvg) %>% 
  summarise(cor_to_4000 = cor(pseudotime_all, pseudotime_4000))

alg_from <- c("monocle", "phenopath")
alg_to <- c("Monocle 2", "PhenoPath")
dataset_from <- c("coad", "brca", "shalek")
dataset_to <- c("COAD", "BRCA", "Shalek et al.")

df_cor$algorithm <- plyr::mapvalues(df_cor$algorithm, from = alg_from, to = alg_to)
df_cor$dataset <- plyr::mapvalues(df_cor$dataset, from = dataset_from, to = dataset_to)

df_cor$hvg <- factor(df_cor$hvg)

ggplot(df_cor, aes(x = hvg, y = abs(cor_to_4000), color = algorithm, group = algorithm)) + 
  geom_point() + facet_wrap(~ dataset) +
  geom_line() +
  labs(subtitle = "Correlation to 4000 highly-variable-gene (HVG) pseudotime",
       x = "Number of HVGs", y = "Correlation") +
  scale_color_brewer(palette = "Set1",
                     name = "Algorithm") +
  theme(strip.background = element_rect(colour = "black", fill = "grey95", 
                                        linetype = "solid", size = 1),
        legend.position = "right",
        axis.text.x = element_text(size = 8))

ggsave("figs/hvg.png", width = 6, height = 2.5)
