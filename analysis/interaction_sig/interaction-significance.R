library(tidyverse)

brca <- read_csv("data/interactions/brca_interactions.csv")

filter(brca, is_sig) %>% 
  arrange(desc(abs(beta))) %>% 
  head(n = 20) %>% 
  select(gene_symbol) %>% 
  write_csv("data/interactions/brca_top10.csv")

coad <- read_csv("data/interactions/coad_interactions.csv")

filter(coad, is_sig) %>% 
  arrange(desc(abs(beta_msi))) %>% 
  head(n = 20) %>% 
  select(gene_symbol) %>% 
  write_csv("data/interactions/coad_top.csv")

