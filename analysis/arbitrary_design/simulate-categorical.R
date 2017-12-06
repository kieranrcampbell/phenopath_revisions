library(tidyverse)
library(viridis)
library(phenopath)

scale_vec <- function(x) (x - mean(x)) / sd(x)


N <- 200
G <- 80

set.seed(5)

pst <- rnorm(N)

x <- sample(1:3, N, replace = TRUE)
x_mat <- as.matrix(model.matrix(~ 0 + factor(x)))

alpha <- matrix(rnorm(G * 3), ncol = 3)

# Make sure beta has one 1 in each row
b <- sample(1:3, G, replace = TRUE)
beta <- as.matrix(model.matrix(~ 0 + factor(b)))

beta_sign <- matrix(sample(c(-1, 1), G * 3, replace = TRUE), ncol = 3)
beta <- beta * beta_sign

lambda <- rnorm(G)

Y <- x_mat %*% t(alpha) + x_mat %*% t(beta) * pst + pst %*% t(lambda)

Y <- apply(Y, 2, function(y) rnorm(N, y, 1))

pca_df <- prcomp(Y)$x[,1:2] %>% 
  as_data_frame() %>% 
  mutate(pst, x)

ggplot(pca_df, aes(x = PC1, y = PC2, color = pst)) +
  geom_point() +
  scale_colour_viridis()

ggplot(pca_df, aes(x = PC1, y = PC2, color = factor(x))) +
  geom_point() +
  scale_colour_brewer(palette = "Set1")

m <- as.matrix(model.matrix(~ 0 + factor(x)))
fit <- phenopath(Y, m)

qplot(pst, trajectory(fit), color = factor(x))

beta_fit <- t(fit$m_beta)

plot(beta[,1], beta_fit[,1])
plot(beta[,2], beta_fit[,2])
plot(beta[,3], beta_fit[,3])

plot(lambda, fit$m_lambda)

x_mat2 <- scale(model.matrix(~ factor(x))[,-1])

s_fit <- t(fit$s_beta)

sig <- significant_interactions(fit, 3)


