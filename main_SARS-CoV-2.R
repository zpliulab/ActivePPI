rm(list = ls())
.libPaths("/home/dell/miniconda3/envs/scMRF/lib/R/library")
library(igraph)
library(ks)
library(doParallel)
library(parallel)
require(doSNOW)

setwd('/home/wcy/scMRF/R_github/R/GitHub-main') # your path

# load data
FINAL_pathway_crop <- readRDS("./SARS-CoV-2/pathway-related/FINAL_pathway_crop.rds")
load("./SARS-CoV-2/MS-related/SARS.data.Rdata") 
source("./code/activePPI.R")


# load data and completion Matrix
data1[data1 == ""] = NA
protein <- impute::impute.knn(as.matrix(data1),k = 5, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
protein <-protein[["data"]]
protein <- apply(protein, 2, as.numeric)
rownames(protein) <- rownames(data1)
protein <- scale(protein)


# calculate all protein kde
Energe.single <- matrix(0,nrow = dim(protein)[1],ncol = dim(protein)[2])
for(i in 1:dim(protein)[1])
{
  fit <- density(x = protein[i,], bw = 'SJ',kernel = 'gaussian', n = 1024)
  ii<- findInterval(as.numeric(protein[i,]), fit$x) # 
  Energe.single[i,] <- fit$y[ii]}
rownames(protein) <- toupper(rownames(protein))
rownames(Energe.single) <- toupper(rownames(protein))


# Pre-calculate distances between proteins
cal_dist = FALSE
if(cal_dist)
{
  dist.list <- c("euclidean", "manhattan", "canberra", "minkowski",'cos')
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
        }}
    } else{
      if(method1 %in% dist.list)
      {prob.pro <- as.matrix(dist(Energe.single, method = method1))
      } else{ prob.pro <- 1 - abs(cor(t(Energe.single), method = method1))}
    }
    rownames(prob.pro) <- rownames(Energe.single)
    colnames(prob.pro) <- rownames(Energe.single)
    save(prob.pro,file = paste0('./SARS-CoV-2/dist/',method1,'_SARS_KS_corr.Rdata'))
  }
}


#Load variables, create output files, and enable parallelism
cls <- makeSOCKcluster(50)
registerDoSNOW(cls)
dist.list <- c("euclidean", "manhattan", "canberra", "minkowski",'cos')

dis <- dist.list[1]
load(paste0('./SARS-CoV-2/dist/',dis,'_SARS_KS_corr.Rdata'))
prob.pro <- prob.pro-min(prob.pro)
prob.pro <- prob.pro/max(prob.pro)
gene <- function(Protein.Pathway){return(toupper(unique(c(Protein.Pathway[[1]][["from"]],Protein.Pathway[[1]][["to"]]))))}
x = data.frame('pathway','pvalue','dist') #create output files
write.table(x,file=paste0('./SARS-CoV-2/result/SARS_result_',dis,'.csv'),sep=",",
            row.names = F,col.names = F,append=F)

pb <- txtProgressBar(max=length(Protein.Pathway), style=3, char = "*",) # 设置进度条
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)


# run activePPI for all pathway
clique.dentisy <- foreach(x = c(1:length(Protein.Pathway)), .options.snow=opts) %dopar% 
  activePPI_pathway(Protein.Pathway[x], 
                    protein[gene(Protein.Pathway[x]),], 
                    Energe.single[gene(Protein.Pathway[x]),], 
                    corr = prob.pro[gene(Protein.Pathway[x]),gene(Protein.Pathway[x])],
                    parall = FALSE, 
                    rtxt=TRUE, 
                    rtxtfile=paste0('./SARS-CoV-2/result/SARS_result_',dis,'.csv'))
names(clique.dentisy) <- names(Protein.Pathway) 

# free memory
close(pb)
stopCluster(cls)
env <- foreach:::.foreachGlobals
rm(list=ls(name=env), pos=env)

# save result
save.image(paste0("./SARS-CoV-2/result/KDE+clique+ks+",dis,"+sars.Rdata"))


#----------------- summary results -----------------
setwd('/home/wcy/scMRF/R_github/R/GitHub-main/SARS-CoV-2')
source("./code/util.R")
load("./SARS-CoV-2/pathway-related/FINAL_pathway_crop.RData")
dist.list <- c("euclidean", "manhattan", "canberra", "minkowski",'cos')
dis <- dist.list[1]


# pathway related-SARS-CoV-2
gold_pathway <- openxlsx::read.xlsx('../SARS-CoV-2//pathway-related/SARS.msigdbr.info.xlsx')
all_pathway.ID <- names(Protein.Pathway)
gold_pathway <- gold_pathway[gold_pathway$ID %in% all_pathway.ID,]

# load activePPI results
clique.dentisy <- read.table(paste0('../SARS-CoV-2/result/SARS_result_',dis,'.csv'),header = TRUE, sep = ',')

# BH correction and extraction components
result_processed <- process_result(clique.dentisy)
clique.dentisy <- result_processed$clique.dentisy

# reorder, extract the rank of the pathway related to SARS-CoV-2
clique.dentisy <- order_result(result_processed)

gold_pathway$value <- clique.dentisy[gold_pathway$ID,'adjpvalue']
gold_pathway$obj <- clique.dentisy[gold_pathway$ID,'obj']
gold_pathway$order <- clique.dentisy[gold_pathway$ID,'order']

save(clique.dentisy,gold_pathway,file = paste0("./SARS-CoV-2/result/FINAL_SARS_activePPI_ks",dis,".Rdata"))

print(gold_pathway)




