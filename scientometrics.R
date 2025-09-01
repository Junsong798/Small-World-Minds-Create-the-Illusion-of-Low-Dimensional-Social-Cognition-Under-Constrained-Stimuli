# Scientometrics of bayesian cognition

# Author: Junsong Lu
# Version: 2024-07-19

# Libraries
library(tidyverse)
library(bibliometrix)


# Sources

# Parameters

# ========================= input data =========================================

file <- paste0("data/wos/savedrecs", " (", 1:8, ").txt")

file <- c("data/wos/savedrecs.txt", file)
file

dat <- convert2df(file, dbsource = "wos", format = "plaintext")

names(dat)

# ========================== data analysis =====================================

#options(width=160)
results <- biblioAnalysis(dat)
summary(results, k = 10, pause = F, width = 130)

plot(x=results, k=10, pause=F)

results

## most cited papers

CR <- citations(dat, field = "article", sep = ";")
cbind(CR$Cited[1:20])

NetMatrix <- biblioNetwork(dat, analysis = "co-citation", 
                           network = "references", sep = ";")

net <- networkPlot(NetMatrix, n = 100, 
                   Title = "Co-Citation Network", 
                   type = "fruchterman", 
                   size.cex=TRUE, size=10, 
                   remove.multiple=FALSE, 
                   labelsize=1,edgesize = 7, edges.min=5)

# =========================== biblioshiny ======================================

biblioshiny()


