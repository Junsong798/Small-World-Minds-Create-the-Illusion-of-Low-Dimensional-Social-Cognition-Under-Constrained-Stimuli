# Dynamic interactive network

# Author: Junsong Lu
# Version: 2025-04-25

# Libraries
library(tidyverse)
library(igraph)

# Sources

# Parameters

# ======================= parameters & functions ===============================

n_steps <- 200
n_node <- 20
noise <- 0.2
I <- 0.4
D <- 0.1
r <- 0
M <- 1
m <- -0.2

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

# ================================ network =====================================

# adjacent matrix

adj_mat <- matrix(runif(n_node^2, 0, 0.5), nrow = n_node, ncol = n_node)

dy_net <- graph_from_adjacency_matrix(adj_mat, mode = "directed",
                                      weighted = TRUE)

plot(dy_net)

E(dy_net)$weight

# ================================= simulation =================================

# initial activation

act_net <- rep(r, times = n_node)

# environmental inputs
env_input <- matrix(0.0, nrow = n_steps + 1, ncol = n_node)

env_input[1, ] <- c(0.45, rep(0, times = n_node - 1))

# save activation patterns

act_array <- matrix(NA, nrow = n_steps + 1, ncol = n_node)
act_array[1, ] <- act_net

# simulations

pb <- txtProgressBar(min = 0, max = n_steps, style = 3, char = "=")


for (t in 1:n_steps) {
  
  setTxtProgressBar(pb, t)
  
  # extract current activation matrix
  current_act_net <- act_array[t, ]
  
  for (n in 1:n_node) {
    # for each node, calculate activation
    net_in <- net_input(n, env_input[t, n], act_array[t, ], noise)
    delta_a <- update_act(net_in, current_act_net[n])
    # save activation to the matrix
    act_array[t + 1, n] <- delta_a + current_act_net[n]
  }
  
}

colnames(act_array) <- paste0("node", 1:n_node)

act_dat <- as_tibble(act_array) %>%
  mutate(time = row_number())

ggplot(act_dat %>%
         pivot_longer(cols = node1:node4,
                      names_to = "node",
                      values_to = "activation")) +
  geom_line(aes(time, activation, group = node, color = node)) +
  theme_bw()

cor.test(act_dat$node1[2:100], act_dat$node2[3:101])
cor.test(act_dat$node1, act_dat$node3)

# ========================= Barabasi-Albert Model ==============================

N <- 100
m_link <- 3   

n_steps <- 500
noise <- 0.1
I <- 0.4
D <- 0.1
r <- 0
M <- 1
m <- -0.2

set.seed(111)
ba_model <- sample_pa(
  n = N,            
  m = m_link,             
  directed = FALSE
)

plot(ba_model, vertex.size = 5, vertex.label = NA)

diam <- diameter(ba_model, directed = FALSE)
diam

deg <- degree(ba_model)

hist(degree(ba_model))
max(deg)
min(deg)

sum(deg)

# 一阶矩（度的均值）
mu_1 <- mean(deg)
mu_1

# 二阶矩（度平方的均值）
mu_2 <- mean(deg^2)
mu_2

(mu_2 - mu_1) / mu_1

# 进一步推导，就是看reach 到k，某个度的节点，相关是多少，不需要最小相关

# ============================ simulate BA model ===============================

adj_mat <- as_adjacency_matrix(ba_model)
adj_mat[adj_mat != 0] <- runif(sum(adj_mat != 0), min = -0.2, max = 0.4)  
adj_mat <- as.matrix(adj_mat)
print(adj_mat[1:10, 1:10])

# initial activation

act_net <- rep(r, times = N)

# environmental inputs
env_input <- matrix(0.0, nrow = n_steps + 1, ncol = N)

env_input[1, ] <- c(0.45, rep(0, times = N - 1))

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

save(act_array, file = "./data/output/activation.RData")
load("./data/output/activation.RData")

colnames(act_array) <- paste0("node", 1:N)

act_dat <- as_tibble(act_array) %>%
  mutate(time = row_number())

node1_neighbors <- extract_neighbors(1)
sub_nodes <- node1_neighbors[1:5]
sub_nodes

# nodes with distance = 2

dist_matrix <- distances(ba_model, v = V(ba_model)[1], 
                         to = V(ba_model), mode = "all")

which(dist_matrix == 1)


# print neighor links
adj_mat[1, sub_nodes]

sd(act_dat$node1)

ggplot(act_dat %>%
         pivot_longer(cols = all_of(c("node1", paste0("node", sub_nodes))),
                      names_to = "node",
                      values_to = "activation")) +
  geom_line(aes(time, activation, group = node, color = node)) +
  theme_bw()

cor.test(act_dat$node1[2:100], act_dat$node10[3:101])
cor.test(act_dat$node1[1:101], act_dat$node2[1:101])
cor.test(act_dat$node1[2:101], act_dat$node11[2:101])

cor.test(act_dat$node1[2:100], act_dat$node45[2:100])

# 如果正反馈很多，即使noise很大，最后相关还是很高
# 这里算错了，应该是改变ev1, 然后计算不同trial之间的相关

lag_corr <- act_dat %>%
  mutate(node1_t = node1,
         node2_t_plus_1 = lead(node5, 1)) %>%
  filter(!is.na(node2_t_plus_1)) %>%   
  summarise(correlation = cor(node1_t, node2_t_plus_1, method = "spearman"))

print(lag_corr)

act_dat %>%
  mutate(node2_t_plus_1 = lead(node2, 1)) %>%
  filter(!is.na(node2_t_plus_1)) %>%
  ggplot(aes(x = node1, y = node2_t_plus_1)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Lag-1 Correlation: node1(t) vs node2(t+1)",
       x = "node1 at time t",
       y = "node2 at time t+1")


# ego centric network

ego_net <- make_ego_graph(ba_model, order = 1, nodes = 1)[[1]]

plot(ego_net,
     edge.color = ifelse(E(ego_net)$weight > 0, "green", "red"), 
     edge.width = abs(E(ego_net)$weight) * 3,    
     vertex.color = ifelse(V(ego_net)$name == "node1", "gold", "lightblue"), 
     vertex.size = 15,
     vertex.label.cex = 0.8,
     main = "Egocentric Network of Node1 (Green: +, Red: -)")
