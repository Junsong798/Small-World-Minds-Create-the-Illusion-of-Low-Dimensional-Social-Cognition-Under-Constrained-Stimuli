# Network simulations with multiple stimuli

# Author: Junsong Lu  
# Version: 2025-05-07

# Libraries
library(tidyverse)
library(igraph)

# Sources

# Parameters

# ================================ functions ===================================

extract_neighbors <- function(id){
  neighbor_id <- which(adj_mat[id, ] != 0)
  return(neighbor_id)
}

net_input <- function(node, env_input, act, noise){

  neighbor_id <- which(adj_mat[node, ] != 0)
  neighbors_w <- adj_mat[node, neighbor_id]
  neighbors_act <- pmax(act[neighbor_id], 0)
  
  # access environmental input
  ext <- env_input
  
  # calculate net input
  net_in <- sum(neighbors_act %*% neighbors_w) + ext
  net_in <- net_in + rnorm(1, 0, noise)
  
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

# ============================ Barabasi-Albert Model ===========================

N <- 50  
m_link <- 4 

set.seed(111)
ba_model <- sample_pa(
  n = N,            
  m = m_link,             
  directed = FALSE
)

plot(ba_model, vertex.size = 5, vertex.label = NA)


diam <- diameter(ba_model, directed = FALSE)
diam

# =============================== simple version ===============================

N <- 10 
m_link <- 2 

set.seed(111)
ba_model <- sample_pa(
  n = N,            
  m = m_link,             
  directed = FALSE
)

plot(ba_model, vertex.size = 5, vertex.label = NA)


diam <- diameter(ba_model, directed = FALSE)
diam

E(ba_model)$weight <- runif(ecount(ba_model), min = -0.4, max = 4)

E(ba_model)$color <- ifelse(E(ba_model)$weight > 0, "blue", "red")
E(ba_model)$width <- abs(E(ba_model)$weight) * 2 

plot(ba_model, 
     vertex.size = 5, 
     vertex.label = NA,
     edge.color = E(ba_model)$color,
     edge.width = E(ba_model)$width,
     layout = layout_with_fr(ba_model, weights = abs(E(ba_model)$weight)))

# ================================= simulation =================================

# parameters
n_steps <- 100
n_stimuli <- 50
noise <- 0.5

I <- 0.4
D <- 0.1
r <- 0
M <- 1
m <- -0.2

adj_mat <- as_adjacency_matrix(ba_model)
adj_mat[adj_mat != 0] <- runif(sum(adj_mat != 0), min = -0.2, max = 0.4)  
adj_mat <- as.matrix(adj_mat)
print(adj_mat[1:10, 1:10])

# initial activation
act_net <- rep(r, times = N)

# environmental inputs

sim_env <- function(i, var_stimuli) {
  env_input <- matrix(0.0, nrow = n_steps + 1, ncol = N)
  a0 <- rnorm(1, mean = 0, sd = var_stimuli)
  env_input[1, ] <- c(a0, rep(0, times = N - 1))
  return(env_input)
}

env_input_list <- map2(1:n_stimuli, rep(0.3, n_stimuli), sim_env)

env_input_list[[1]]

# simulations

single_simu <- function(n_steps, env_input, noise){
  # n_step: a scalar
  # env_input: a (n_step + 1) by N (# of node) matrix
  # noise: a scalar indicating noise variance, 0.5
  
  # save activation 
  act_array <- matrix(NA, nrow = n_steps + 1, ncol = N)
  act_array[1, ] <- act_net
  
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
  return(act_array)
}

simu_dat <- tibble(n_steps = n_steps,
                   env_input_list = env_input_list,
                   noise = noise) %>%
  mutate(act = pmap(list(n_steps, env_input_list, noise),
                    single_simu,
                    .progress = TRUE))

# find neighbors
extract_neighbors(1)
sub_nodes <- extract_neighbors(1)[1:5]
sub_nodes

# find k-distance neighbors
dist_matrix <- distances(ba_model, v = V(ba_model)[1], 
                         to = V(ba_model), mode = "all")

which(dist_matrix == 2)

adj_mat[1, sub_nodes]

simu_val <- simu_dat %>%
  mutate(node1_99 = map_dbl(act, ~ .[98, 1]),
         node2_100 = map_dbl(act, ~ .[99, 2]))

simu_val

ggplot(simu_val) +
  geom_point(aes(node1_99, node2_100))

cor.test(simu_val$node1_99, simu_val$node2_100)

ego_net <- make_ego_graph(ba_model, order = 1, nodes = 1)[[1]]

plot(ba_model)
