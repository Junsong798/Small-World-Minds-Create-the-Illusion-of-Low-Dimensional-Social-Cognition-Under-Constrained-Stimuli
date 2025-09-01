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

N <- 40  
m_link <- 3 

set.seed(111)
ba_model <- sample_pa(
  n = N,            
  m = m_link,             
  directed = FALSE
)

plot(ba_model, vertex.size = 5)


diam <- diameter(ba_model, directed = FALSE)
diam

E(ba_model)$weight <- rnorm(ecount(ba_model), 0, 0.25)
E(ba_model)$weight <- runif(ecount(ba_model), 0, 0.5)

# =========================== spread activation ================================
library(spreadr)

V(ba_model)

initial_df <- data.frame(node = '4', activation = 50, stringsAsFactors = FALSE)
initial_df

initial_df <- data.frame(node = c('4', '5', '44', '81', '12'), 
                         activation = c(50, 10, 25, 25, 35), stringsAsFactors = FALSE)
initial_df

result <- spreadr::spreadr(start_run = initial_df, decay = 0.05,
                           retention = 0.6, suppress = 0,
                           network = ba_model, time = 30) 

head(result)

a1 <- data.frame(node = '4', activation = 20, time = 0) # add back initial activation at t = 0
result_t0 <- rbind(result,a1)
ggplot(data = result_t0, aes(x = time, y = activation, color = node, group = node)) +
  geom_point() + geom_line() + ggtitle('unweighted, undirected network')

dat_act <- as_tibble(result) %>%
  pivot_wider(names_from = node,
              values_from = activation)

dat_act %>%
  rstatix::cor_test(`4`, `15`)

# extract neighbors

adj_mat <- as_adj(ba_model)

extract_neighbors <- function(id){
  neighbor_id <- which(adj_mat[id, ] != 0)
  return(neighbor_id)
}

dist_matrix <- distances(ba_model, v = V(ba_model), 
                         to = V(ba_model), mode = "all")

# batch correlations

batch_corr <- function(target_node){
  # target node is a number
  neighbor_id <- which(adj_mat[target_node, ] != 0)
  neighbor_id <- as.character(neighbor_id)
  
  node_list <- 1:N
  node_list <- node_list[-target_node]
  target_node <- as.character(target_node)
  node_list <- as.character(node_list)
  
  corr_list <- map_dbl(node_list, ~ {
    cor_result <- dat_act %>% rstatix::cor_test(all_of(target_node), all_of(.x))
    cor_result$cor
    })
  return(corr_list)
}

batch_corr(4)

dat_dis <- tibble(corr = batch_corr(4),
                  dist = dist_matrix[4, ])