# Derived correlation

# Author: Junsong Lu
# Version: 2025-05-25

# Libraries
library(tidyverse)

# Sources

# Parameters

# =========================== basic functions ==================================


power_law_moment <- function(m, gamma, k_min, k_max = Inf) {
  if (gamma == m + 1) {
    stop("Moment diverges: gamma = m + 1 leads to a logarithmic divergence.")
  }
  
  if (is.infinite(k_max)) {
    if (gamma > m + 1) {
      # Approximate result for large k_max
      moment <- (gamma - 1) / (gamma - 1 - m) * k_min^m
    } else {
      warning("Moment diverges for gamma <= m + 1 with infinite k_max.")
      moment <- Inf
    }
  } else {
    # Normalization constant
    C <- (gamma - 1) * k_min^(gamma - 1) / (1 - (k_min / k_max)^(gamma - 1))
    
    # Compute moment using definite integral
    moment <- C * (k_max^(m + 1 - gamma) - k_min^(m + 1 - gamma)) / (m + 1 - gamma)
  }
  
  return(moment)
}

power_law_moment(2, 3.1, 4, Inf)
power_law_moment(1, 3.1, 4, Inf)

# ======================== estimate expected distance ==========================

# first network
N <- 500
m_link <- 3  

set.seed(111)
ba_model <- sample_pa(
  n = N,            
  m = m_link,             
  directed = FALSE
)

plot(ba_model, vertex.size = 0.5, vertex.label = NA)

degree(ba_model)
diameter(ba_model)

deg <- degree(ba_model)

hist(degree(ba_model))
max(deg)
min(deg)

m <- sum(deg) / 2

c1 <- mean(deg)
c1

c2 <- mean(deg^2) - c1
c2


distance_probability <- function(s, t, c1, c2, m) {
  exponent <- s + t - 2
  term <- (c2 / c1)^exponent * c1^2
  prob <- (1 - (c2/c1)^2 / (2 * m))^term
  return(prob)
}

distance_probability(2, 1, c1, c2, m)

# expected link between two specific nodes

scaling <- c2 / c1
scaling

basic_prob <- function(k1, k2, m){
  return(k1 * k2 / (2 * m))
}

basic_prob(20, 20, m)

# Create a grid of k1 and k2 values from 1 to 20
k1_values <- seq(1, 25, by = 1)
k2_values <- seq(1, 25, by = 1)

# Generate a data frame with the grid of k1 and k2, and the calculated value for basic_prob * scaling
results <- expand.grid(k1 = k1_values, k2 = k2_values)
results$basic_prob_scaled <- basic_prob(results$k1, results$k2, m) * scaling  # Assume m = 1 for simplicity

# Plot using ggplot2
exp_link <- ggplot(results, aes(x = k1, y = k2, z = basic_prob_scaled)) +
  geom_tile(aes(fill = basic_prob_scaled), color = "white") +  # Use geom_tile for a heatmap-style plot
  scale_fill_viridis_c() +  # Use a color scale
  labs(x = "k1", y = "k2", fill = "Value") +
  theme_minimal() +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 10, vjust = 0.85),
        legend.position = "top",
        legend.key.width = unit(0.85, "cm"))

N2 <- 200
m_link2 <- 2  

set.seed(111)
ba_model2 <- sample_pa(
  n = N2,            
  m = m_link2,             
  directed = FALSE
)

deg2 <- degree(ba_model2)
m2 <- sum(deg2) / 2

c1_2 <- mean(deg2)
c1_2

c2_2 <- mean(deg2^2) - c1_2
c2_2

scaling2 <- c2_2 / c1_2

k1_values_2 <- seq(1, 25, by = 1)
k2_values_2 <- seq(1, 25, by = 1)

# Generate a data frame with the grid of k1 and k2, and the calculated value for basic_prob * scaling
results_2 <- expand.grid(k1 = k1_values_2, k2 = k2_values_2)
results_2$basic_prob_scaled <- basic_prob(results_2$k1, results_2$k2, m) * scaling2  # Assume m = 1 for simplicity

exp_link2 <- ggplot(results_2, aes(x = k1, y = k2, z = basic_prob_scaled)) +
  geom_tile(aes(fill = basic_prob_scaled), color = "white") +  # Use geom_tile for a heatmap-style plot
  scale_fill_viridis_c(option = "plasma") +  
  labs(x = "k1", y = "k2", fill = "Value") +
  theme_minimal() +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 10, vjust = 0.85),
        legend.position = "top",
        legend.key.width = unit(0.85, "cm"))


plot_exp_link <- cowplot::plot_grid(exp_link, exp_link2, ncol = 2)
ggsave("exp_link.png", plot_exp_link, width = 7, height = 6, 
       dpi = 700, bg = "white")
