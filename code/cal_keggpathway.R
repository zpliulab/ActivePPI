rm(list = ls())
setwd("/home/wcy/scMRF/R/SARS-Cov-2")
library(KEGGREST)
library(base)

hsa.pathway.database = readRDS('hsa.pathway.database.rds')

# 下面代码是将每个小网络写入表格中
lpathway <- matrix(0, ncol = length(hsa.pathway.database), nrow = 1)
for (i in 1:length(hsa.pathway.database))
{
  cat(i,' ')
  A <- hsa.pathway.database[[i]][[1]][["GENE"]]
  B <- matrix(0, ncol = length(A)/2, nrow = 1)
  for (j in 1:(length(A)/2))
  {
    B[j] <- as.numeric(A[(j*2)-1])
  }
  colnames(B)<-NULL
  rownames(B)<-NULL
  
  if (i == 1)
  {lpathway[i] <- length(A)/2
  write.table(B," hsa.pathway.csv", col.names = FALSE,row.names = FALSE, quote=F, sep = ",")}
  else if (!is.null(A))
  {
    lpathway[i] <- length(A)/2
    write.table(B," hsa.pathway.csv", col.names = FALSE,row.names = FALSE, append=T, quote=F, sep = ",")}
}

# 下面代码是查看数据中存在的TF、Gene、或者miRNA名字

load('SARS.data.Rdata')
Gene_name <- toupper(rownames(data1))

# 下面代码是将每个小网络与ReguNetwork对应起来
human_network <- read.table(file='./Regnetwork/2022.human.source',sep='\t')
human_network_TF_ID <- as.matrix(human_network[,2]) #  因为有RNA，所以这是一个字符串
human_network_Gene_ID <- as.matrix(human_network[,4]) 

# 正式循环每个小网络
d <- 1
path_network <- c()
path_network.name <- c()
for (i1 in 1:345)
{
  A <- hsa.pathway.database[[i1]][[1]][["GENE"]]   
  if (!is.null(A) && ((length(A)/2)>3)  && ((length(A)/2)<200) )
  {
    B_ID <- matrix(0, ncol = length(A)/2, nrow = 1)  # B_ID是KEGG第i1个子网
    B_name <- matrix(0, ncol = length(A)/2, nrow = 1) # B_ID是KEGG第i1个子网对应的名字
    for (j in 1:(length(A)/2))
    {
      B_ID[j] <- as.numeric(A[(j*2)-1])
      B_name[j] <- toupper(substring(A[(j*2)],1,gregexpr(';', A[(j*2)])[[1]][1]-1))
    }
    colnames(B_ID)<-NULL
    rownames(B_ID)<-NULL
    # 先删除小网络中多余的基因（根据IAV数据）
    B_GENE <- matrix(0, ncol = length(B_name), nrow = 1)
    for(i2 in 1:length(B_GENE))   #匹配IAV中的基因名和小网络的基因名
    {
      if (length(grep(paste("^",B_name[i2],"$",sep =''),Gene_name))>0)
      {
        B_GENE [i2] <- grep(paste("^",B_name[i2],"$",sep =''),Gene_name)
      }
    }
    B_GENE_fitter_ID <- B_ID[B_GENE>0]      #网络B中 基因ID号
    B_GENE_fitter_name <- B_name[B_GENE>0]  #网络B中 基因名字
    B_GENE_fitter_S1 <- B_GENE[B_GENE>0]    #网络B中 数据的位置
    if((length(B_GENE_fitter_ID)/length(B_ID))<0.5)
    {print('出现了太少的匹配结果!') }
    networknl <- c()

    for(i in 1:length(B_GENE_fitter_ID)) #遍历每个基因  看看是不是TF
    {
      number <- grep(paste("^",B_GENE_fitter_ID[i],"$",sep =''),human_network_TF_ID)
      if (length(number)>0)  #找到了一个TF
      {
        for (z in 1:length(number))  #对于每个小网络中的基因 查看其是否为靶基因
        {
          networkn <- c()
          number2 <- grep(paste("^",human_network_Gene_ID[number[z]],"$",sep =''),B_GENE_fitter_ID)
          if (length(number2)>0)
          {
            networkn[1] <- B_GENE_fitter_name[i]           #调控基因
            networkn[3] <- B_GENE_fitter_S1[i]             #调控基因数据位置
            networkn[2] <- B_GENE_fitter_name[number2[1]]  #靶基因
            networkn[4] <- B_GENE_fitter_S1[number2[1]]    #靶基因数据位置
            networkn[5] <- 1
          #  networkn[8] <- paste('第',i1,'子网络', seq='')
            networknl <- rbind(networknl,networkn)
            print(paste('第',d,'个关联'))
            d <- d+1
          }
        }
      }
    }
    if(!is.null(networknl))
    {
      if(dim(networknl)[1]>2)
      {
        colnames(networknl)<-c('from','to','data.of.from', 'data.of.to','score')
        path_network <- c(path_network,list(data.frame(networknl)))
        path_network.name <-c(path_network.name,hsa.pathway.database[[i1]][[1]][["ENTRY"]][["Pathway"]])
        
      }
    }
  }
}
names(path_network) <- path_network.name
# 将path_network写入表格每个sheet中(openxlsx)
save(path_network,file = "Kegg_SARS_network.RData")
save.image("cal_pathway.RData")




