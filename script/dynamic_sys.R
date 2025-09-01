# Dynamic systems

# Author: Junsong Lu
# Version: 2025-05-06

# Libraries
library(tidyverse)

# Sources

# Parameters

# ============================================================================

set.seed(123)
n_sims <- 5000
max_t  <- 50
w      <- 0.35
sigma_e<- 0.9
c_const<- 0.1

a_mat <- matrix(NA, nrow=n_sims, ncol=max_t+1)
b_mat <- matrix(NA, nrow=n_sims, ncol=max_t+1)
a_mat[,1] <- rnorm(n_sims,0,1)
b_mat[,1] <- 0

for(t in 2:(max_t+1)){
  a_prev <- a_mat[,t-1]; b_prev <- b_mat[,t-1]
  b_mat[,t] <- c_const + w*a_prev + rnorm(n_sims,0,sigma_e)
  a_mat[,t] <- c_const + w*b_prev + rnorm(n_sims,0,sigma_e)
}

cor(a_mat[,51], b_mat[,50])

rho_t <- sapply(2:(max_t+1), function(tt) {
  cor(a_mat[, tt], b_mat[, tt])
})


