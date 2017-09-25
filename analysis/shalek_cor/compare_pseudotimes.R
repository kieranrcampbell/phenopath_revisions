

library(tidyverse)
library(cowplot)

theme_set(theme_cowplot(font_size = 11))

hvg_dir <- "data/shalek_cor"

sce <- readRDS("data/paper-scesets/sce_shalek_qc.rds")

time_df <- select(pData(sce), time) %>% 
  as_data_frame() %>% 
  mutate(sample = paste0("sample_", seq_len(ncol(sce))))

all_files <- dir(hvg_dir, full.names = TRUE)
all_files <- all_files[grepl("pseudotime", all_files)]


df_list <- lapply(all_files, read_csv)
df_list <- lapply(df_list, function(d) { d$hvg <- as.character(d$hvg); d$x <- sce$x; d} )

df <- bind_rows(df_list)

df <- inner_join(df, time_df, by = "sample")

df_maxmin <- group_by(df, algorithm, hvg) %>% 
  summarise(max_pst = max(pseudotime, na.rm = TRUE), min_pst = min(pseudotime, na.rm = TRUE))

df <- inner_join(df, df_maxmin, by = c("algorithm", "hvg"))
df_norm <- mutate(df, pseudotime_norm = (pseudotime - min_pst) / (max_pst - min_pst))

plt1 <- ggplot(df_norm, aes(x = algorithm, y = pseudotime_norm, fill = time)) +
  geom_boxplot() +
  facet_wrap(~ hvg)


# R2s ---------------------------------------------------------------------

get_R2 <- function(pseudotime, time) {
  time <- as.numeric(gsub("h", "", time))
  fit <- lm(pseudotime ~ time)
  s <- summary(fit)
  s$r.squared
  #abs(cor(time, pseudotime, use = "na"))
}

df_R2 <- group_by(df_norm, algorithm, hvg) %>% 
  summarise(R2_to_time = get_R2(pseudotime, time))

# df_R2$hvg <- plyr::mapvalues(df_R2$hvg, from = "all", to = as.character(nrow(sce)))
df_R2$hvg <- factor(as.numeric(df_R2$hvg))

# df_R2 <- filter(df_R2, !grepl("init_monocle", algorithm) & !grepl("init_time", algorithm))

alg_from <- c("dpt", "monocle2", "phenopath_init_pc1", "tscan")
alg_to <- c("DPT", "Monocle 2", "PhenoPath", "TSCAN")
df_R2$algorithm <- plyr::mapvalues(df_R2$algorithm, from = alg_from, to = alg_to)

ggplot(df_R2, aes(x = hvg, y = R2_to_time, group = algorithm, color = algorithm)) +
  geom_point() + geom_line() +
  labs(y = "R2 to true time", x = "Number of highly variable genes") +
  scale_color_brewer(palette = "Set1", name = "Algorithm")

# 
# all_files <- dir(hvg_dir, full.names = TRUE)
# all_files <- all_files[grepl("interaction", all_files)]
# 
# 
# df_list <- lapply(all_files, read_csv)
# df_list <- lapply(df_list, function(d) { d$hvg <- as.character(d$hvg);  d} )
# 
# df_int <- bind_rows(df_list)
# 
# df_int_g <- group_by(df_int, algorithm, hvg) %>% 
#   summarise(pct_interactions = 100 * mean(interactions),
#             num_interactions = sum(interactions))
# 
# df_int_g$hvg <- plyr::mapvalues(df_int_g$hvg, from = "all", to = as.character(nrow(sce)))
# df_int_g$hvg <- factor(as.numeric(df_int_g$hvg))
# 
# ggplot(df_int_g, aes(x = hvg, y = pct_interactions, color = algorithm, group = algorithm)) +
#   geom_point() + geom_line()
# 
# ggplot(df_int_g, aes(x = hvg, y = num_interactions, color = algorithm, group = algorithm)) +
#   geom_point() + geom_line()


