# Assortative mating I

# Author: Junsong Lu
# Version: 2025-04-24

# Libraries
library(tidyverse)

# Sources

# Parameters

# ============================== parameters ====================================

popsize <- 200
mr <- 10 * (.006/100) # mutation rate variability added to traits and preferences each generation
generations <- 1000
selstrength <- 0.10 # selection strength (bioenergy range of fitness)
ratio_male <- 0.5
n_trait <- 10

# =============================== functions ====================================

# Gene Generation

gene_generate <- function(popsize){
  if (popsize > 0) {
    # generate traits
    traits_gene <- matrix(runif(popsize * n_trait * 100, 0, 0.08), nrow = popsize)
    colnames(traits_gene) <- paste('trait_gene', 1:ncol(traits_gene), sep = "")
    
    return(as.data.frame(traits_gene))
  }
}

# Agent generation

# traits_gene <- gene_generate(200)
# sex <- rep(c(0, 1), times = popsize / 2)
# sex <- 1

agent_generate <- function(traits_gene, sex){

  # Determine each agent's traits by summing their relevant genes and then name them appropriately.
  traits <- sapply(1:n_trait, function(x) rowSums(traits_gene[, (1+(x-1)*100):(x*100)]))
  colnames(traits) <- paste('trait',1:ncol(traits), sep = "")
  
  # differential energy calculation
  mate_value <- vector("numeric", length = popsize)
  
  if (sex == 0) {
    # male's mate value does not depend on attractiveness
    mate_value <- rowMeans(n_trait - (8 - traits[, 2:n_trait]))
  } else {
    # female's mate value is solely determined by attractiveness
    mate_value <- n_trait - (8 - traits[, 1])
  }
  
  return(as.data.frame(cbind(traits, sex, mate_value)))
}

#traits_gene <- as.data.frame(cbind(traits, sex, mate_value))

# Attraction

#male_mv <- 5.97
#female_mv <- 9.26

attraction <- function(prefer, trait){
  # we assume that attraction is based on mate value, the closer the
  # mate value between two individuals, the more likely they will become couples
  
  # here, the perceive trait of a single person is the input, then, traits from
  # all available targets are used as the 2nd input
  max_dist <- sqrt(8^2)
  #raw_dist <- abs(prefer - trait)
  
  #Calculates the Euclidean distance between the agent's preferences and each mate's traits
  #att <- 10 * (-raw_dist + max_dist) / max_dist
  att <- apply(trait, 1, function(x) 10 * (-abs(prefer - x) + max_dist) / max_dist)
  
  return(att)
}

#attraction(male_mv, female_mv)

# Reproduction

reproduce_genes <- function(mother, father, mr){
  
  #Offspring inherit trait and preference genes randomly from either of their parents
  traits_gene <- rbinom(1000, 1, .5)
  traits_gene <- ifelse(traits_gene == 1, as.numeric(father), as.numeric(mother))
  
  #Add a set amount of variability to genes
  traits_gene <- traits_gene + rnorm(length(traits_gene), 0, mr)
  
  return(traits_gene)
}


# Offspring Reproduction

#repro <- 200
#offspring_gene <- traits_gene[1:repro,]

reproduce_offspring <- function(offspring_gene){
  # empty matrix for genes
  offspring <- matrix(NA, nrow = nrow(repro), 2 + n_trait)
  traits <- sapply(1:n_trait, function(x) rowSums(offspring_gene[, (1+(x-1)*100):(x*100)]))
  
  # calculate traits
  offspring[, 1:10] <- traits
  colnames(offspring[, 1:10]) <- paste('trait',1:n_trait, sep = "")
  
  # sex
  offspring[, 11] <- sample(rep(c(0, 1), each = (nrow(offspring) / 2)))
  
  # mate value
  offspring[, 12] <- n_trait - (8 - offspring[, 1]) # female method
  sex_ind <- offspring[, 11] == 0
  offspring[sex_ind, 12] <- rowMeans(n_trait - (8 - offspring[sex_ind, 2:10]))
  
  return(as.data.frame(offspring))
}

# ============================ simulation (single) =============================

# generate agent genes
male_genes <- gene_generate(popsize/2)
female_genes <- gene_generate(popsize/2)

# generate agents
males <- agent_generate(male_genes, 1)
females <- agent_generate(female_genes, 0)

# give each agent a unique PIN
males$PIN <- sample(1:nrow(males), nrow(males), replace = F)
females$PIN <- sample((nrow(males) + 1):(nrow(males) + nrow(females)),
                      nrow(females),
                      replace = F)

male_genes$PIN <- males$PIN
female_genes$PIN <- females$PIN

## ============================ save statistics ================================

# matrix for factor analysis results
facvar <- matrix(NA, 2, (generations + 1))

datam <- males
dataf <- females

# extracts a general factor from traits
mfa_result <- psych::fa(datam[,1:n_trait], nfactors = 1, rotate = "varimax", fm = "ml")
ffa_result <- psych::fa(dataf[,1:n_trait], nfactors = 1, rotate = "varimax", fm = "ml")

# calculate the proportion of variance accounted for by the general factor
var_explained_m <- sum(mfa_result$loadings[, 1]^2) / n_trait
var_explained_f <- sum(ffa_result$loadings[, 1]^2) / n_trait

facvar[1,1] <- var_explained_m
facvar[2,1] <- var_explained_f

## =========================== start simulations ===============================

pb <- txtProgressBar(
  min = 0, 
  max = generations,  
  style = 3,
  char = "="  
)

for (g in 1:generations) {
  
  setTxtProgressBar(pb, g)
  ##### computing attraction #####
  
  ## males
  ## Calculate how attracted each male is to each female (trait 2-9)
  male_att_matrix <- t(apply(males,1,function(q) attraction(q[12], females[12])))
  
  
  ## females
  ## Calculate how attracted each female is to each male (only trait 1 attractiveness)
  female_att_matrix <- t(apply(females,1,function(q) attraction(q[12], males[12])))
  
  
  ##### mate selection #####
  
  #Mutual attraction matrix
  #Computes the mutual attraction matrix--how mutually attracted each possible couple would be.
  mut_att_matrix <- as.data.frame(male_att_matrix * t(female_att_matrix))
  
  #Renames the columns to be the female PIN values
  colnames(mut_att_matrix) <- females$PIN
  
  #Renames the rows to be male PIN values
  rownames(mut_att_matrix) <- males$PIN
  
  #A blank dataframe for storing the paired couples
  # 100 * 46, 46 = 23 +23, 23 = 10 trait + 10 prefer + gender + id + energy
  pairs <- matrix(NA, nrow = min(nrow(males), nrow(females)),
                  ncol = (ncol(males) + ncol(females)))
  pairs <- as.data.frame(pairs)
  
  #Iterates through the mutual attraction matrix
  for (p in 1:(popsize/2)) {
    #Determines the location of the maximum mutual attraction value in the matrix
    # 这里打印了3个值，第一个是行名，然后是坐标；但实际上是2个值，只有坐标，所以下一步除以2就是看几个最大值
    max_mut <- which(mut_att_matrix == max(mut_att_matrix, na.rm = T),
                     arr.ind = T)
    
    #If there is more than one location equal to the maximum, picks a random location
    if (length(max_mut)/2 > 1) {
      max_mut <-max_mut[sample(nrow(maxmut), 1, replace = F), ]
    }
    
    #Determines which male and which female were actually paired
    # 抽出单独的男女，一行，23列
    paired_male <- males[males$PIN == rownames(mut_att_matrix)[max_mut[1]], ]
    paired_female <- females[females$PIN == colnames(mut_att_matrix)[max_mut[2]], ]
    
    #Places the paired male and female in the pairs dataframe
    pairs[p, ] <- cbind(paired_male, paired_female)
    
    #Removes the paired male and female from the mutual attraction matrix
    mut_att_matrix[max_mut[1], ] <- NA
    mut_att_matrix[, max_mut[2]] <- NA
    
  }
  
  #Just names the pairs matrix appropriately so it's easier to read. Probably unnecessary.
  colnames(pairs) <- c(colnames(males), colnames(females))
  
  ##### reproduce #####
  
  #Determine each couple's summed mate values
  pairs_mv <- (pairs[, 12] + pairs[, 25])
  
  #Rescale these values just to make imposing selection strength easier
  pairs_mv <- pairs_mv - min(pairs_mv)
  pairs_mv <- 100 * pairs_mv/max(pairs_mv)
  
  #Add the selstrength constant, which controls the strength of selection on energy values
  pairs_mv <- pairs_mv + (100 / selstrength)
  
  #Determines how often each couple gets to reproduce by sampling randomly from couples in proportion to their energy values
  # 这里抽样popsize个确保下一代的人数也是popsize
  repro <- pairs[sample(nrow(pairs), popsize, replace = T, prob=pairs_mv), ]
  
  #Generate offspring genes based on the repro dataframe
  # 这个是200 * 1000的大矩阵
  off_genes <- t(apply(repro, 1, function(x) reproduce_genes(female_genes[female_genes$PIN == x[26], ],
                                                             male_genes[male_genes$PIN == x[13], ], mr)))
  
  #Convert genes to a dataframe and label for ease of use later
  off_genes <- as.data.frame(off_genes)
  colnames(off_genes) <- colnames(female_genes)[1:1000]
  
  #Produce offspring agents based on these genes
  offspring <- reproduce_offspring(off_genes)
  
  #Rename everything for ease of use
  colnames(offspring) <- colnames(repro)[1:12]
  
  ##### death #####
  #Give each offspring a PIN
  offspring$PIN <- sample(1:nrow(offspring), nrow(offspring), replace = F)
  
  #Save the offspring over the parent dataframes. This kills the parents
  males <- subset(offspring, offspring$sex == 0)
  females <- subset(offspring, offspring$sex == 1)
  
  #Save the genes of the offspring
  male_genes <- off_genes[offspring$sex == 0, ]
  female_genes <- off_genes[offspring$sex == 1, ]
  
  #Store the agent PINs in the gene dataframes
  male_genes$PIN <- males$PIN
  female_genes$PIN <- females$PIN
  
  ##### data saving #####
  
  #Split the last pairs dataframe into males and females
  datam <- pairs[,1:10]
  dataf <- pairs[,14:23]
  
  # extracts a general factor from traits
  mfa_result <- psych::fa(datam, nfactors = 1, rotate = "varimax", fm = "ml")
  ffa_result <- psych::fa(dataf, nfactors = 1, rotate = "varimax", fm = "ml")
  
  # calculate the proportion of variance accounted for by the general factor
  var_explained_m <- sum(mfa_result$loadings[, 1]^2) / n_trait
  var_explained_f <- sum(ffa_result$loadings[, 1]^2) / n_trait
  
  facvar[1, g+1] <- var_explained_m
  facvar[2, g+1] <- var_explained_f
  
}

save(facvar, file = "./data/output/d_factor.RData")


rownames(facvar) <- c("male", "female")
colnames(facvar) <- c(1, 2:(generations+1))

ggplot(facvar %>%
         t() %>%                        
         as_tibble(rownames = "generation") %>%     
         mutate(generation = as.numeric(generation)) %>%
         pivot_longer(cols = c(male, female),
                      values_to = "var_explained",
                      names_to = "gender")) +
  geom_line(aes(generation, var_explained, group = gender, color = gender)) +
  theme_bw()

mean(facvar[2,])
