# Figures

# Author: Junsong Lu
# Version: 2025-07-11

# Libraries
library(tidyverse)
library(igraph)

# Sources

# Parameters

# =========================== Barabasi-Albert ==================================

N <- 1000
m_link <- 2  

set.seed(111)
ba_model <- sample_pa(
  n = N,            
  m = m_link,             
  directed = FALSE
)

plot(ba_model, vertex.size = 0.5, vertex.label = NA)

write_graph(ba_model, file = "ba_1000.gml", format = "gml")


# dynamic figure

max_nodes <- 100             # Total number of nodes to grow
steps <- c(2, 3, 4, 5, 10, 20, 50, 100)  # Growth steps to visualize

# Generate full graph first
set.seed(1)
full_graph <- sample_pa(n = max_nodes, 
                        power = 1,
                        m = m_link,
                        directed = FALSE)

# Compute layout for the full graph (fixed layout)
layout_fixed <- layout_with_fr(full_graph)
size_list <- c(12, 12, 12, 12, 9, 9, 5, 5)

# Set up plotting area
png(filename = "growing_pattern3.png", width = 8, height = 6, 
    units = "in", res = 700)

par(mfrow = c(2, 4), mar = c(0.25, 0.25, 1, 0.25))  # Grid layout

# Plot subgraphs at different steps using the same layout
for (t in seq_along(steps)) {
  i <- which(steps == steps[t])
  sub_g <- induced_subgraph(full_graph, 1:steps[t])
  layout_sub <- layout_fixed[1:steps[t], , drop = FALSE]  # Always returns a matrix
  par(cex.main = 1.35)
  
  degrees <- degree(sub_g)
  
  plot(sub_g,
       layout = layout_sub,
       vertex.size = degrees * (1.8 + 4.5 * (1 - floor((t - 1) / 4))),#size_list[i],
       vertex.label = NA,
       vertex.color = "#A75263", #"#61678B",
       vertex.frame.color = "#BCB9D8",  # Change vertex frame color
       vertex.frame.width = 1,     # Ensure frame width is applied
       #edge.curved = 0.2,
       edge.width = 1,
       edge.color = "grey80",
       main = paste("t = ", steps[t]))
  #box(lwd = 1.2)
}
dev.off()
par(mfrow = c(1, 1))

# ======================= preferential attachment ==============================

# Function to simulate Chinese Restaurant Process
crp_simulation <- function(n_customers, alpha = 1) {
  tables <- integer(n_customers)  # Table assignments for each customer
  n_tables <- 0
  
  for (i in 1:n_customers) {
    # Count the number of customers at each table
    table_counts <- table(tables[tables != 0])  # Exclude empty tables
    table_probs <- table_counts / (i - 1 + alpha)  # Probability of choosing existing tables
    new_table_prob <- alpha / (i - 1 + alpha)  # Probability of choosing a new table
    
    # Combine the probabilities (existing tables + new table)
    table_probs <- c(table_probs, new_table_prob)
    
    # Sample a table based on the probabilities
    sampled_table <- sample(1:(n_tables + 1), size = 1, prob = table_probs)
    
    if (sampled_table == (n_tables + 1)) {
      n_tables <- n_tables + 1  # A new table is chosen
      tables[i] <- n_tables
    } else {
      tables[i] <- names(table_counts)[sampled_table]
    }
  }
  return(tables)
}

# Simulate the CRP process for 50 customers
n_customers <- 500
alpha <- 1
customer_table_assignments <- crp_simulation(n_customers, alpha)

# Create a data frame for ggplot
crp_data <- data.frame(
  Customer = 1:n_customers,
  Table = customer_table_assignments
)

# Plot using ggplot2
ggplot(crp_data, aes(x = Table, y = Customer)) +
  geom_point() +
  theme_bw() +
  labs(title = "A draw from a Chinese restaurant process",
       x = "Tables",
       y = "Customers") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# Function to simulate Preferential Attachment
m <- 1
pa_simulation <- function(n_links, m) {
  #n_links <- 50
  select_node <- rep(0, times = n_links)
  # count the # links per node
  n_node <- 0
  
  for (i in 1:n_links) {
    i = 3
    link_counts <- table(select_node[select_node != 0])
    link_prob <- link_counts / (m * (i - 1) + 1)
    link_self <- 1 / (m * (i - 1) + 1)
    
    link_prob <- c(link_prob, link_self)
    
    sample_link <- sample(1:i, size = 1, prob = link_prob)
    
    if (sample_link == i) {
      select_node[i] <- sample_link
      n_node <- n_node + 1
    } else {
      select_node[i] <- sample_link
    }
  }
  return(select_node)
}

pa_simulation(50)



pa_simulation <- function(n_links) {
  # 边界情况处理
  if (n_links == 0) return(numeric(0))
  if (n_links == 1) return(1)  # 仅1条边：节点2连接节点1
  
  # 初始化
  select_node <- rep(NA, n_links)  # 存储每条边的连接目标
  deg <- numeric(2 * n_links + 1)  # 度数向量（索引=节点ID）
  
  # 初始网络：2个节点，1条边（节点1↔节点2）
  select_node[1] <- 1  # 第一条边：节点2连接到节点1
  deg[1] <- 1          # 节点1度数=1
  deg[2] <- 1          # 节点2度数=1
  n_node <- 2           # 当前节点数=2
  
  # 模拟新增边（从第2条边开始）
  for (t in 2:n_links) {
    # 1. 计算概率 = 节点度数 / 总度数
    total_deg <- sum(deg[1:n_node])
    probs <- deg[1:n_node] / total_deg
    
    # 2. 按概率选择目标节点
    target <- sample(1:n_node, size = 1, prob = probs)
    select_node[t] <- target  # 记录连接目标
    
    # 3. 更新度数
    deg[target] <- deg[target] + 1  # 目标节点度数+1
    n_node <- n_node + 1            # 创建新节点
    deg[n_node] <- 1                # 新节点初始度数=1
  }
  
  return(select_node)
}

# 模拟50条边的网络
set.seed(123)  # 确保可重复性
network_links <- pa_simulation(50)

# =========================== communities ======================================

set.seed(123)

# Number of nodes
num_nodes <- 60

# Generate a random graph with 3 communities
block_sizes <- c(13, 13, 14, 10, 10)  # Adjusted to sum to 40

# Generate a random graph with 3 communities
com_mat <- matrix(rep(0.1, times = 5 * 5), nrow = 5)
diag(com_mat) <- 0.8
g <- sample_sbm(num_nodes, pref.matrix = com_mat,
                block.sizes = block_sizes)


# Add weights to the edges
E(g)$weight <- runif(ecount(g), min = 0.1, max = 1)

# Get the adjacency matrix
adj_matrix <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)

# Display the adjacency matrix
svg("weighted_network.svg", width = 10, height = 10, bg = NA)
plot(g, vertex.size = 10, vertex.label = NA, 
     vertex.color = "#FF9300",
     edge.width = E(g)$weight * 1,
     vertex.frame.color = "white", vertex.frame.width = 2.5)
dev.off()

community_colors <- rep(c("#DE6B69", "#61708B", "#87D1D1", "#A2549B",
                          "#E7C55A"), times = block_sizes)
plot(g, vertex.size = 10, vertex.label = NA, 
     vertex.color = community_colors,
     edge.width = E(g)$weight * 1,
     vertex.frame.color = "white", vertex.frame.width = 2.5)


layout <- layout_with_fr(g)  # This provides a layout for the graph

# Set up the PNG output
png("combined_network2.png", width = 12, height = 6, units = "in", res = 700)

# Set up a 1x2 plot layout (two plots side by side)
par(mfrow = c(1, 2), mar = c(2, 2, 1, 1))  # Adjust margins for better spacing

# First plot: Default color
plot(g, vertex.size = 10, vertex.label = NA, 
     vertex.color = "#FF9300",
     edge.width = E(g)$weight * 1,
     vertex.frame.color = "white", vertex.frame.width = 2.5,
     layout = layout)  # Use the same layout for both plots

# Second plot: Community-based color
#community_colors <- rep(c("#DE6B69", "#61708B", "#87D1D1"), times = block_sizes)
plot(g, vertex.size = 10, vertex.label = NA, 
     vertex.color = community_colors,
     edge.width = E(g)$weight * 1,
     vertex.frame.color = "white", vertex.frame.width = 2.5,
     layout = layout)  # Use the same layout for both plots

# Close the PNG device
dev.off()

# ============================ 2-hop network ===================================

set.seed(1)
n <- 100
sample_net <- sample_pa(n = n,
                        m = 2,
                        directed = FALSE)

plot(sample_net, vertex.size = 7, vertex.label = NA)

first_node <- 2
one_hop <- neighbors(sample_net, 2)
two_hop <- unlist(map(one_hop, ~ neighbors(sample_net, .)))
two_hop <- two_hop[!(two_hop %in% one_hop)]
two_hop <- two_hop[two_hop != 2]
two_hop <- unique(two_hop)

palf <- colorRampPalette(c("gray80", "dark red")) 
color_list <- palf(10)

vertex_colors <- rep("gray90", n)  
vertex_colors[one_hop] <- "#A75A5A"  
vertex_colors[two_hop] <- "#BD9E9E"  
vertex_colors[first_node] <- "#8B0000" 

neighbor_ind <- vector("numeric", length = n)
for (i in 1:n) {
  if (i %in% one_hop){
    neighbor_ind[i] <- 1
  } else if(i %in% two_hop) {
    neighbor_ind[i] <- 2
  } else{
    neighbor_ind[i] <- 3
  }
}

neighbor_ind[2] <- 0

neighbor_ind
V(sample_net)$neighbor_ind <- neighbor_ind

l <- layout_with_fr(sample_net)

write_graph(sample_net, file = "ba_activ.gml", format = "gml")

plot(sample_net,
     layout = l,
     vertex.label = NA,
     vertex.frame.color = "white",
     vertex.frame.width = 1,
     vertex.color = vertex_colors,
     vertex.size = 8)

degree(sample_net)

w <- c(6,6,5,2,2,1)
b <- c(3,3,2,4,4,3)
cor.test(w,b)

# sample calculation 1

log(100)/log(log(100))

# 参数设置
ki <- 10
kj <- 15
N <- 100
k_mean <- 5.88
k2_mean <- 54.38

a <- 0  # 下限
b <- 0.5   # 上限

n <- 3  # 最大路径长度

# 计算常量
M <- N * k_mean
r <- (k2_mean - k_mean) / k_mean

# 计算 E[w^d] for d = 1 to n
E_wd <- function(d, a, b) {
  (b^(d + 1) - a^(d + 1)) / ((d + 1) * (b - a))
}

# 向量化计算
d_vals <- 1:n
E_wd_vals <- sapply(d_vals, E_wd, a = a, b = b)
r_vals <- r^(d_vals - 1)

# 计算总和
effect <- (ki * kj / M) * sum(r_vals * E_wd_vals)

# 计算 marginal effect: 第3项减去第2项的贡献
effect_d2 <- (ki * kj / M) * sum(r_vals[1:2] * E_wd_vals[1:2])
marginal_d3 <- effect - effect_d2

# 打印结果
cat("Expected total effect:", effect, "\n")
cat("Marginal effect of increasing path length to 3:", marginal_d3, "\n")

# sample calculation 2

corr_AXd <- function(w, d, sigma_A, sigma_eps) {
  if (abs(w^2 - 1) < 1e-8) {
    # 特殊情况：w^2 = 1
    corr <- sigma_A / sqrt(sigma_A^2 + d * sigma_eps^2)
  } else {
    numerator <- w^d * sigma_A
    denominator <- sqrt(w^(2 * d) * sigma_A^2 +
                          sigma_eps^2 * (1 - w^(2 * d)) / (1 - w^2))
    corr <- numerator / denominator
  }
  return(corr)
}

w <- 0.5
d <- 3
sigma_A <- 1
sigma_eps <- 0.25

corr <- corr_AXd(w, d, sigma_A, sigma_eps)
cat("Corr(A, X_d) =", corr, "\n")

w <- 0.5
d <- 2
sigma_A <- 1
sigma_eps <- 0.25

corr <- corr_AXd(w, d, sigma_A, sigma_eps)
cat("Corr(A, X_d) =", corr, "\n")

