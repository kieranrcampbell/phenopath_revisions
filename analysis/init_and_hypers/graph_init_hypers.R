

library(tidyverse)
library(cowplot)

theme_set(theme_cowplot(font_size = 11))

hvg_dir <- "data/init_and_hypers"

control_df <- read_csv("data/init_and_hypers/control.csv")

all_files <- dir(hvg_dir, full.names = TRUE)
all_files <- all_files[!grepl("control", all_files)]

df_list <- lapply(all_files, read_csv)

df <- bind_rows(df_list)

cor_df <- group_by(df, elbo_tol, tau_alpha, z_init, ab_beta_ratio) %>% 
  summarise(cor_to_control = cor(pseudotime, control_df$pseudotime))

yl <- "Absolute correlation\nto default pseudotime"

## Sort the x axes the way we want them

cor_df$elbo_tol_fct <- factor(cor_df$elbo_tol,
                              levels = sort(unique(cor_df$elbo_tol), decreasing = TRUE))

cor_df$z_init_chr <- plyr::mapvalues(cor_df$z_init,
                                     from = c("monocle", "pc1", "time"),
                                     to = c("Monocle 2", "PC1+noise", "Capture times\n(scaled)"))

## The graphs

plt1 <- ggplot(cor_df, aes(x = elbo_tol_fct, y = abs(cor_to_control))) +
  geom_boxplot() +
  labs(x = "ELBO change for convergence (%)", y = yl)

plt2 <- ggplot(cor_df, aes(x = factor(tau_alpha), y = abs(cor_to_control))) +
  geom_boxplot() +
  labs(x = expression(tau[alpha]), y = yl)

plt3 <- ggplot(cor_df, aes(x = z_init_chr, y = abs(cor_to_control))) +
  geom_boxplot() +
  labs(x = "Initialisation for z", y = yl)

plt4 <- ggplot(cor_df, aes(x = factor(ab_beta_ratio), y = abs(cor_to_control))) +
  geom_boxplot() +
  labs(x = expression(a[beta]~"/"~b[beta]), y = yl)

plot_grid(plt1, plt2, plt3, plt4, labels = 'AUTO')

ggsave("figs/supp_robustness_to_init_hyper.png", width = 8, height = 6)
