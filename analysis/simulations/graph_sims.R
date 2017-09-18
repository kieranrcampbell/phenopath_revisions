library(dplyr)
library(tidyr)
library(ggplot2)

source("analysis/simulations/simulate_funcs.R")

set.seed(1234L)

N <- 1000
pst <- runif(N)
x <- sample(c(0,1), N, size = N)
k <- 10

y <- f_pseudotime(pst, k)

y1 <- fx_1(pst, x, k)
y2 <- fx_2(pst, x, k)
y3 <- fx_3(pst, x, k)
y4 <- fx_4(pst, x, k)

df <- data_frame(y, y1, y2, y3, y4)

i_names <- c("No interaction",
             "Interaction type 1", 
               "Interaction type 2",
               "Interaction type 3",
               "Interaction type 4")

names(df) <- i_names

df <- mutate(df, pst, x) %>% 
  gather(interaction_type, expression, -pst, -x) %>% 
  mutate(x = factor(x))
df$interaction_type <- factor(df$interaction_type, levels = i_names)

ggplot(df, aes(x = pst, y = expression, color = x, group = x)) +
  geom_line(size = 1.5, alpha = 0.8) + 
  facet_wrap(~ interaction_type, nrow = 1) +
  scale_colour_brewer(palette = "Set1", name = "Covariate") +
  theme(strip.background = element_rect(fill = "white"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 10)) +
  labs(y = "Simulated expression") +
  xlab(expression(Pseudotime ~ symbol('\256')))

saveRDS(last_plot(), "figs/simulation_example.rds")







  