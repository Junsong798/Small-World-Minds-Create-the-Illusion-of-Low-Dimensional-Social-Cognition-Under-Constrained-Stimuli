# Dynamic interactive network

# Author: Junsong Lu
# Version: 2025-07-04

# Libraries
library(tidyverse)
library(igraph)

# Sources

# Parameters

# =============================== functions ====================================

extract_neighbors <- function(id){
  neighbor_id <- which(adj_mat[id, ] != 0)
  return(neighbor_id)
}

net_input <- function(node, env_input, act, noise){
  # act_net is the matrix of current activation
  # access neighbor nodes' id and weights
  neighbor_id <- which(adj_mat[node, ] != 0)
  neighbors_w <- adj_mat[node, neighbor_id]
  #neighbors_act <- act[neighbor_id]
  neighbors_act <- pmax(act[neighbor_id], 0)
  
  # access environmental input
  ext <- env_input
  
  # calculate net input
  net_in <- sum(neighbors_act %*% neighbors_w) + ext
  net_in <- net_in + rnorm(1, 0, noise)
  
  #print(paste("Net input:", net_in))
  
  return(net_in)
}

update_act <- function(net_input, current_act){
  # differential updating
  if(net_input > 0){
    delta_a <- I * (M - current_act) * net_input - D * (current_act - r) 
  }else{
    delta_a <- I * (current_act - m) * net_input - D * (current_act - r)
  }
  
  return(delta_a)
}


# ========================= Barabasi-Albert Model ==============================

N <- 100
m_link <- 3   

n_steps <- 500
noise <- 0.05
I <- 0.4
D <- 0.1
r <- 0
M <- 1
m <- 0

set.seed(111)
ba_model <- sample_pa(
  n = N,            
  m = m_link,             
  directed = FALSE
)

plot(ba_model, vertex.size = 8)

diam <- diameter(ba_model, directed = FALSE)
diam

deg <- degree(ba_model)

hist(degree(ba_model))
max(deg)
min(deg)

sum(deg)

# first moment
mu_1 <- mean(deg)
mu_1

# second moment
mu_2 <- mean(deg^2)
mu_2

(mu_2 - mu_1) / mu_1


# ============================ simulate BA model ===============================

adj_mat <- as_adjacency_matrix(ba_model)
adj_mat[adj_mat != 0] <- runif(sum(adj_mat != 0), min = -0.2, max = 0.4)  
adj_mat <- as.matrix(adj_mat)

# initial activation

act_net <- rep(r, times = N)

# environmental inputs
env_input <- matrix(0.0, nrow = n_steps + 1, ncol = N)

env_input[1, ] <- c(0.5, rep(0, times = N - 1))

# random environmental input
env_input <- matrix(rnorm((n_steps + 1) * N, 0, 1), nrow = n_steps + 1, ncol = N)

# save activation patterns

act_array <- matrix(NA, nrow = n_steps + 1, ncol = N)
act_array[1, ] <- act_net

# simulations

pb <- txtProgressBar(min = 0, max = n_steps, style = 3, char = "=")


for (t in 1:n_steps) {
  
  setTxtProgressBar(pb, t)
  
  # extract current activation matrix
  current_act_net <- act_array[t, ]
  
  for (n in 1:N) {
    # for each node, calculate activation
    net_in <- net_input(n, env_input[t, n], act_array[t, ], noise)
    delta_a <- update_act(net_in, current_act_net[n])
    # save activation to the matrix
    act_array[t + 1, n] <- delta_a + current_act_net[n]
  }
  
}

colnames(act_array) <- paste0("node", 1:N)

act_dat <- as_tibble(act_array) %>%
  mutate(time = row_number())

node1_neighbors <- extract_neighbors(1)
sub_nodes <- node1_neighbors[1:10]
sub_nodes

# nodes with distance = 2

dist_matrix <- distances(ba_model, v = V(ba_model)[1], 
                         to = V(ba_model), mode = "all")

which(dist_matrix == 3)


# print neighor links
adj_mat[1, sub_nodes]

act_plot <- ggplot(act_dat %>%
         filter(time <= 100) %>%
         pivot_longer(cols = all_of(c("node1", paste0("node", sub_nodes))),
                      names_to = "node",
                      values_to = "activation")) +
  geom_line(aes(time, activation, group = node, color = node),
            linewidth = 1) +
  geom_point(aes(time, activation, group = node, color = node),
             alpha = 0.7) +
  scale_color_brewer(palette = "Set3") +
  theme_bw() +
  labs(x = "Time", y = "Activation") +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.position = "none")

ggsave("act_trend.png", plot = act_plot, width = 7, height = 6, dpi = 700)

cor.test(act_dat$node1[2:100], act_dat$node10[3:101])
cor.test(act_dat$node1[1:101], act_dat$node2[1:101])
cor.test(act_dat$node1[2:101], act_dat$node11[2:101])

cor.test(act_dat$node1[2:100], act_dat$node52[2:100])

# ====================== simulate multiple stimuli =============================

# generate multiple stimuli

n_stimuli <- 100
stimuli_list <- runif(n_stimuli)

act_array_list <- array(NA, dim = c(n_stimuli, n_steps + 1, N))
dim(act_array_list)

for (i in 1:n_stimuli) {
 
  # initial activation
  act_net <- rep(r, times = N)
  
  # environmental inputs
  env_input <- matrix(0.0, nrow = n_steps + 1, ncol = N)
  
  env_input[1, ] <- c(stimuli_list[i], rep(0, times = N - 1))
  
  # save activation patterns
  act_array <- matrix(NA, nrow = n_steps + 1, ncol = N)
  act_array[1, ] <- act_net
  
  # simulations
  
  for (t in 1:n_steps) {
    # extract current activation matrix
    current_act_net <- act_array[t, ]
    
    for (n in 1:N) {
      # for each node, calculate activation
      net_in <- net_input(n, env_input[t, n], act_array[t, ], noise)
      delta_a <- update_act(net_in, current_act_net[n])
      # save activation to the matrix
      act_array[t + 1, n] <- delta_a + current_act_net[n]
    }
    
  }
  
  colnames(act_array) <- paste0("node", 1:N)
  
  act_array_list[i, , ] <- act_array
}

#save(act_array_list, file = "./data/output/simu_multi_stimuli.RData")
load("./data/output/simu_multi_stimuli.RData")

# calculate mean activation for each node in each iteration

mean_act_node <- matrix(NA, nrow = n_stimuli, ncol = N)

for (i in 1:n_stimuli) {
  
  tmp_act_array <- act_array_list[i, , ]
  tmp_act_array <- apply(tmp_act_array[c(2:51),], 2, mean)
  mean_act_node[i, ] <- tmp_act_array
}

mean_act_node
cor.test(mean_act_node[,1], mean_act_node[,21])

plot(mean_act_node[,1], mean_act_node[,5])

bi_cor_plot <- ggplot(as_tibble(mean_act_node)) +
  geom_jitter(aes(V1, V21), shape = 21, stroke = 1, alpha = 0.85,
             size = 3.5, color = "white", fill = "#CF91A3",
             width = 0.01,
             height = 0.01) +
  geom_smooth(aes(V1, V21), method = "lm", color = "#995A88",
              linewidth = 1.25, alpha = 0.25) +
  theme_bw() +
  labs(x = "Attractiveness", y = "Intelligence") +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12))

ggsave("bi_cor.png", plot = bi_cor_plot, width = 7, height = 6, dpi = 700)

act_dat_tmp <- as_tibble(act_array_list[3,,]) %>%
  mutate(time = row_number()) %>%
  rename_with(~ paste0("node", 1:N), .cols = -time)

ggplot(act_dat_tmp %>%
         pivot_longer(cols = all_of(c("node1", paste0("node", sub_nodes))),
                      names_to = "node",
                      values_to = "activation")) +
  geom_line(aes(time, activation, group = node, color = node)) +
  theme_bw()

# factor analysis

fa_dat <- mean_act_node[, sub_nodes]
cor_mat <- cor(fa_dat)

mean(cor_mat[lower.tri(cor_mat)])
sd(cor_mat[lower.tri(cor_mat)])

png("cor_mat_transparent.png", 
    width = 7, height = 6, 
    units = "in", res = 700,
    bg = "transparent")
corrplot::corrplot(
  cor_mat, 
  method = "color",
  order = "hclust",
  type = "upper",
  tl.pos = "n",
  col = colorRampPalette(c("steelblue", "white", "darkred"))(200)
)
dev.off() 

png("cor_mat_constrained_env.png", 
    width = 7, height = 6, 
    units = "in", res = 700,
    bg = "transparent")
corrplot::corrplot(
  cor_mat, 
  method = "color",
  order = "hclust",
  type = "upper",
  tl.pos = "n",
  col = colorRampPalette(c("#52863C", "white", "#D33866"))(200)
)
dev.off() 

library(psych)

fa.parallel(fa_dat, fa="fa")

# KMO and MSA
KMO(fa_dat)

Nfacs <- 1
 
fit <- factanal(fa_dat, Nfacs, rotation = "promax")

print(fit, digits=2, cutoff=0.3, sort=TRUE)

# ============== simulate multiple stimuli (random input) ======================

# generate multiple stimuli

n_stimuli <- 100
stimuli_list <- runif(n_stimuli)

act_array_list <- array(NA, dim = c(n_stimuli, n_steps + 1, N))
dim(act_array_list)

for (i in 1:n_stimuli) {
  
  # initial activation
  act_net <- rep(r, times = N)
  
  # environmental inputs
  env_input <- matrix(rnorm((n_steps + 1) * N, 0, 1), 
                      nrow = n_steps + 1, ncol = N)
  
  # save activation patterns
  act_array <- matrix(NA, nrow = n_steps + 1, ncol = N)
  act_array[1, ] <- act_net
  
  # simulations
  
  for (t in 1:n_steps) {
    # extract current activation matrix
    current_act_net <- act_array[t, ]
    
    for (n in 1:N) {
      # for each node, calculate activation
      net_in <- net_input(n, env_input[t, n], act_array[t, ], noise)
      delta_a <- update_act(net_in, current_act_net[n])
      # save activation to the matrix
      act_array[t + 1, n] <- delta_a + current_act_net[n]
    }
    
  }
  
  colnames(act_array) <- paste0("node", 1:N)
  
  act_array_list[i, , ] <- act_array
}

#save(act_array_list, file = "./data/output/simu_multi_stimuli_random.RData")
load("./data/output/simu_multi_stimuli_random.RData")

# calculate mean activation for each node in each iteration

mean_act_node <- matrix(NA, nrow = n_stimuli, ncol = N)

for (i in 1:n_stimuli) {
  
  tmp_act_array <- act_array_list[i, , ]
  tmp_act_array <- apply(tmp_act_array[c(2:51),], 2, mean)
  mean_act_node[i, ] <- tmp_act_array
}

mean_act_node
cor.test(mean_act_node[,1], mean_act_node[,21])

plot(mean_act_node[,1], mean_act_node[,5])

# bi_cor_plot <- ggplot(as_tibble(mean_act_node)) +
#   geom_point(aes(V1, V21), shape = 21, stroke = 1.34, alpha = 0.85,
#              size = 3.5, color = "#8B0000", fill = "white") +
#   geom_smooth(aes(V1, V21), method = "lm", color = "#8B0000",
#               linewidth = 1.25,
#               alpha = 0.2) +
#   theme_bw() +
#   labs(x = "Attractiveness", y = "Intelligence") +
#   theme(axis.title.x = element_text(face = "bold", size = 12),
#         axis.title.y = element_text(face = "bold", size = 12))
# 
# ggsave("bi_cor.png", plot = bi_cor_plot, width = 7, height = 6, dpi = 700)

act_dat_tmp <- as_tibble(act_array_list[3,,]) %>%
  mutate(time = row_number()) %>%
  rename_with(~ paste0("node", 1:N), .cols = -time)

ggplot(act_dat_tmp %>%
         pivot_longer(cols = all_of(c("node1", paste0("node", sub_nodes))),
                      names_to = "node",
                      values_to = "activation")) +
  geom_line(aes(time, activation, group = node, color = node)) +
  theme_bw()

# factor analysis

fa_dat <- mean_act_node[, sub_nodes]
cor_mat <- cor(fa_dat)

mean(cor_mat[lower.tri(cor_mat)])
sd(cor_mat[lower.tri(cor_mat)])

png("cor_mat_random_env.png", 
    width = 7, height = 6, 
    units = "in", res = 700,
    bg = "transparent")
corrplot::corrplot(
  cor_mat, 
  method = "color",
  order = "hclust",
  type = "upper",
  tl.pos = "n",
  col = colorRampPalette(c("#52863C", "white", "#D33866"))(200)
)
dev.off() 


library(psych)

fa.parallel(fa_dat, fa="fa")

# KMO and MSA
KMO(fa_dat)

Nfacs <- 1

fit <- factanal(fa_dat, Nfacs, rotation = "promax")

print(fit, digits=2, cutoff=0.3, sort=TRUE)
