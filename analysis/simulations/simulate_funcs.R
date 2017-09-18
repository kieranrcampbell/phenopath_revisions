



sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

#' All cells change uniformly over pseudotime (given noise)
f_pseudotime <- function(pst, k, t0 = 0.5, mu0 = 4) {
  2 * mu0 * sigmoid(k * (pst - t0))
}

f_pseudotime_sample <- function(pst, k, t0 = 0.5, mu0 = 4) {
  2 * mu0 * sigmoid(sample(c(-1, 1), 1) * k * (pst - t0))
}


#' Here expression begins off on both branches and 
#' turns on on one branch
fx_1 <- function(pst, x, k, t0 = 0.5, mu0 = 4) {
  k <- abs(k) # Make sure k is positive
  
  y <- f_pseudotime(pst, k, t0, mu0)
  
  which_off <- sample(c(0, 1), 1)
  y[x == which_off] <- 0
  y
}

#' Here expression begins on on both branches and turns off 
#' on one branch
fx_2 <- function(pst, x, k, t0 = 0.5, mu0 = 4) {
  k <- -abs(k) # Make sure k is negative
  
  y <- f_pseudotime(pst, k, t0, mu0)
  
  which_off <- sample(c(0, 1), 1)
  y[x == which_off] <- 2 * mu0
  y
}

#' Here expression starts on/off for each branch and converges
fx_3 <- function(pst, x, k, t0 = 0.5, mu0 = 4) {
  abs_k <- abs(k) 
  
  y_off <- f_pseudotime(pst, -abs_k, t0, mu0 / 2) + mu0
  y_on <- f_pseudotime(pst, abs_k, t0, mu0 / 2)
  
  which_off <- sample(c(0, 1), 1)
  
  y <- rep(0, length(pst))
  y[x == which_off] <- y_off[x == which_off]
  y[x != which_off] <- y_on[x != which_off]
  
  y
}

#' Here expression starts together and divides
fx_4 <- function(pst, x, k, t0 = 0.5, mu0 = 4) {
  abs_k <- abs(k) 
  
  y_off <- f_pseudotime(pst, -abs_k, t0, mu0 / 2)
  y_on <- f_pseudotime(pst, abs_k, t0, mu0 / 2) + mu0
  
  which_off <- sample(c(0, 1), 1)
  
  y <- rep(0, length(pst))
  y[x == which_off] <- y_off[x == which_off]
  y[x != which_off] <- y_on[x != which_off]
  
  y
}
