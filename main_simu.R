rm(list = ls())
library(igraph)
library(ks)
library(ggplot2)

setwd('/home/wcy/scMRF/R_github/R/GitHub-main/simu') # your path
source("./code/activePPI.R")
load("./data/simu_data_0.Rdata")

# load data and calculate Kernel Density
G.exp.n.sum <- simu_data$exp
gra <- simu_data$gra
Energe.single <- matrix(0,nrow = dim(G.exp.n.sum)[1],ncol = dim(G.exp.n.sum)[2])
for(i in 1:dim(G.exp.n.sum)[1])
{
  fit <- density(G.exp.n.sum[i,], bw = 'SJ',kernel = 'gaussian', n = 1024)
  ii<- findInterval(as.numeric(G.exp.n.sum[i,]), fit$x) # 
  Energe.single[i,] <- (fit$y[ii])
}
rownames(Energe.single) <- rownames(G.exp.n.sum)

# Pre-calculate distances between proteins
if(TRUE)
{
  dist.list <- c("canberra", "cos","euclidean", "manhattan", "minkowski")
  for(method1 in dist.list)
  {
    if(method1 == "cos")
    { 
      GL <- rownames(Energe.single)
      prob.pro <- matrix(0,nrow=dim(Energe.single)[1],ncol=dim(Energe.single)[1])
      for(Gi1 in 1:(length(GL)-1))
      {
        for(Gii1 in (Gi1+1):length(GL))
        {
          Gi <- GL[Gi1]
          Gii <- GL[Gii1]
          prob.pro[Gi1,Gii1] <- 1-sum(t(Energe.single[Gi,])*Energe.single[Gii,])/sqrt((sum(Energe.single[Gi,]^2))*sum(Energe.single[Gii,]^2))
          prob.pro[Gii1,Gi1] <- prob.pro[Gi1,Gii1]
        }
      }
    } else{if(method1 %in% dist.list){prob.pro <- as.matrix(dist(Energe.single, method = method1))}
    }
    rownames(prob.pro) <- rownames(Energe.single)
    colnames(prob.pro) <- rownames(Energe.single)
    save(prob.pro,file = paste0('./dist/',method1,'_corr.Rdata'))
  }
}

# run activePPI (calculate network activity)
for (dis in dist.list)
{
  load(paste0('./dist/',dis,'_corr.Rdata'))
  prob.pro <- prob.pro-min(prob.pro)
  prob.pro <- prob.pro/max(prob.pro)
  result_all <- c()
  for(j in 1:30)
  {
    result <- activePPI_network(gra, G.exp.n.sum,Energe.single, prob.pro) 
    result_all <- c(result_all,result$pvalue)
  }
  print(sprintf('distance = %s, pvalue: mean = %.3f, std = %.3f',dis,floor(mean(result_all)*1e3)/1e3,floor(sd(result_all)*1e3)/1e3))
}




