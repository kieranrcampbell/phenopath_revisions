library(readr)
library(ggplot2)

allcor_file <- "data/simulations/all_pseudotime_correlations.csv"

df_split <- read_csv(allcor_file)

ggplot(df_split, aes(x = as.factor(100 * p), 
                     y = kendall_correlation, fill = alg)) +
  geom_boxplot() + # geom_line() +
  labs(y = "Absolute Kendall correlation\nto true pseudotime",
       x = "% genes covariate interaction") 
