#################################################
# The main content of activePPI algorithm
#
# Function 1 (estimating pathway activity):
#   result <- activePPI_pathway(Protein.Pathway, protein, Energe.single, corr = prob.pro, parallel = FALSE, verbose = FALSE, rtxt = FALSE, rtxtfile='SARS_result.csv')
# Input:
#   Protein.Pathway : a list containing gene set networks.It must contain "from" and "to" columns, eg.
#       form to
#       A B
#       A C
#       B D
#   protein:   a (n*m) protein mass spectrometry information matrix
#   Energe.single:   Kernel density of each protein (n*m)
#   corr :   distance matrix between proteins (n*n)
#   parallel (default = FALSE):   whether to run in parallel (unstable, it is recommended to run one by one, you can add a foreach parallel outside)
#   verbose (default = FALSE):   whether to load the progress bar
#   rtxt (default=FALSE):    whether to output documentation
#   rtxtfile: output file name
#
# Output:
#  result: a list, including pvalue (result$pvalue), 
#         joint probability density of the pathway (result$dis), 
#         energy of each clique in each network (result$re), 
#         optimal network architecture (result$gra)
#
# Function 2 (estimating network activity, this is a test on a simulated example):
#   result <- activePPI_network(gra, protein, Energe. single, prob. pro)
#
# Input:
#   gra: a network of igraph structures
#   The rest of the settings are consistent with the function one
#
# Output:
#   result: a list, including pvalue (result$pvalue), 
#          joint probability density of each network (result$prob), 
#          energy of each clique in each network (result$prob.list),
#          optimal network architecture (result$ gra)
##################################################
activePPI_pathway <- function(Protein.Pathway, protein,Energe.single, corr = prob.pro, 
                    parall = FALSE,verbose = FALSE, rtxt = FALSE, rtxtfile='SARS_result.csv')
{
methdod = 'kde'
Energe.singleN <- Energe.single
df.network<-data.frame(from = toupper(Protein.Pathway[[1]][,'from']),to = toupper(Protein.Pathway[[1]][,'to']))
df.network <- df.network[df.network[,1] %in% rownames(protein),]
df.network <- df.network[df.network[,2 ] %in% rownames(protein),]
nodes <- unique(c(df.network[,1],df.network[,2]))
df.network[,3] <- paste0(df.network[,1],'-',df.network[,2])
D1 <- duplicated(df.network[,3])
df.network <- df.network[D1 == FALSE,] 
df.network <- df.network[,-3] 
gra<-igraph::graph_from_data_frame(d=df.network,vertices=nodes,directed=FALSE)
gra1 <- as.matrix(igraph::get.adjacency(gra))
gra1[gra1>1] <- 1
gra <- igraph::graph_from_adjacency_matrix(gra1,diag = FALSE, mode = 'undirected')
gra1 <- as.matrix(igraph::get.adjacency(gra))
gra.list <- list()

# rewire多次网络结构
number <- 1000
k <- 1
d <- 1
iter <- igraph::vcount(gra)
if(iter > 100)
{
  number <- 100
}
if((iter > 100)&&(names(Protein.Pathway)=='hsa05171'))
{
  number <- 1000
}

while(k <= number)
{
  if(k != 1)
  {
    bug <<- 0
    tryCatch(
       {gra<-igraph::rewire(gra,with=igraph::keeping_degseq(loop=TRUE,niter=iter))}, 
      error = function(e) {bug <<- 1}
    )
    if(bug == 1)
    { print(paste0(names(Protein.Pathway),'不能重构网络，算法终止！'))
      return(list(pvalue=1,dis=Inf))
      }
  }
  
  gra2 <- as.matrix(igraph::get.adjacency(gra))
  if(k == 1 | length(which((gra2-gra1)!=0))>0 )
  {
    d <- 1
    gra.list <- c(gra.list,list(gra))
    k <- k +1
  }else{d <- d+1}
  
  if(d>1e4) 
  { print(paste0(names(Protein.Pathway),' 无法重构出更多的网络结构，算法终止！'))
    return(list(pvalue=1,dis=Inf))
  }
}

# 判断网络结构的rewire数量是否合格
if(length(gra.list) < number)
{
  return(list(pvalue=1,dis=Inf))
}

# 根据parall参数和节点数量决定是否使用并行计算
result <- c()
mean.dist <- c()
if(!parall || length(V(gra))<100)
{
  if(verbose)
  {## 新建一个进度条弹窗
   pbs <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                        max = number, # Maximum value of the progress bar
                        style = 3,    # Progress bar style (also available style = 1 and style = 2)
                        width = 50,   # Progress bar width. Defaults to getOption("width")
                        char = "=")   # Character used to create the bar 
  }
  for(i in 1:length(gra.list))
  {
    gra <- gra.list[[i]]
    glist <- names(igraph::V(gra))
    result1 <- caldentisy(gra,protein[glist,],Energe.singleN[glist,],methdod, prob.pro=prob.pro)
    mean.dist <- c(mean.dist,mean(result1$prob.list))
    result <- c(result,result1$clique.dentisy)
    if(verbose)
    {setTxtProgressBar(pbs,i)}
  }
  if(verbose)
  {close(pbs)}
  
}else{
  # 并行计算时  需要开启30个内核
  
  library(doParallel)
  library(parallel)
  require(doSNOW)
  cls <- makeSOCKcluster(30)
  registerDoSNOW(cls)
  source("/home/wcy/scMRF/R/caldentisy.R")
  clusterExport(cls, "caldentisy")
  clusterExport(cls, "prob.pro")
  pbs <- txtProgressBar(max=length(Protein.Pathway), style=3, char = "*",) # 设置进度条
  progress <- function(n) setTxtProgressBar(pbs, n)
  opts <- list(progress=progress)
  result1 <- foreach(x = c(1:length(gra.list)), .options.snow=opts) %dopar% 
            caldentisy(gra.list[[x]],
                       protein[names(igraph::V(gra.list[[x]])),],
                       Energe.singleN[names(igraph::V(gra.list[[x]])),],
                       methdod, prob.pro=prob.pro)
  close(pbs) 
  stopCluster(cls)
  env <- foreach:::.foreachGlobals 
  rm(list=ls(name=env), pos=env)
  for(i in 1:length(result1))
  {
    mean.dist <- c(mean.dist,mean(result1[[i]]$prob.list))
    result <- c(result,result1[[i]]$clique.dentisy)
  }
}

# 计算完成 总结pvalue 
result <- as.matrix(result)
pvalue <- (length(which(result< result[1]))+1)/number
if((pvalue == 0)||(pvalue == 0.01 && number == 100))
{
  pvalue <- 0.001
}
# while dist=Inf, then pvalue=1
if(is.infinite(mean.dist[1]))
{
  pvalue <- 1
}
# write text
if(rtxt)
{
  x = data.frame(names(Protein.Pathway),pvalue,mean.dist[1])
  write.table(x,file = rtxtfile,sep=",",
              row.names = F,col.names = F,append=T)
}

# end
return(list(pvalue = pvalue, dis = mean.dist[1],re = result, gra = gra.list))
}



activePPI_network <- function(gra,protein,Energe.single,prob.pro)
{
  method = 'kde'
  result2 <- list()
  result3 <- list()
  Energe.singleN <- Energe.single
  gra1 <- igraph::get.adjacency(gra)
  number <- 1000
  k = 1
  result <- c()
  while(k <= number)
  {
    if(k != 1)
    {
      bug <<- 0
      tryCatch(
        {gra<-igraph::rewire(gra,with=igraph::keeping_degseq(loop=TRUE,niter=igraph::vcount(gra)))}, 
        error = function(e) {bug <<- 1}
      )
      if(bug == 1)
      {return(pvalue <- 1)}
    }
    gra2 <- igraph::get.adjacency(gra)
    if(k == 1 | !any(mapply(identical,list(gra1),list(gra2)))) 
    {
      k <- k +1
      gra <- graph_from_adjacency_matrix(
        as.matrix(gra2),
        mode = "undirected",
        weighted = NULL,
        diag = TRUE,
        add.colnames = NULL,
        add.rownames = NA
      )
      result1 <- caldentisy(gra,protein,Energe.singleN,method,prob.pro)
      result <- c(result,result1$clique.dentisy) 
      result2 <- c(result2,list(result1$prob.list)) 
      result3 <- c(result3,list(gra))
    }
  }
  result <- as.matrix(result)
  pvalue <- (length(which((result<result[1]) == TRUE))+1)/number
  return(list(pvalue = pvalue,prob = result,prob.list = result2,gra = result3))
}

caldentisy <- function(network,protein,Energe.single,methdod = 'kde',prob.pro)
{
  dist.list <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  prob.list <- c()
  network <- igraph::as.undirected(network, mode ="mutual",edge.attr.comb="c")
  VI <- igraph::vcount(network)
  if(VI < 50)
  {
    Cliq <- igraph::max_cliques(network)     
  }else{
    Cliq <- igraph::max_cliques(network,max = 10)     
  }
  
  if (length(Cliq) > 0)
  {
    # order Cliq
    GN <- as.double(lapply(Cliq, length))
    
#    GN <- order(-1*GN)
#    Cliq <- Cliq[GN]
    Cliq <- Cliq[1:min(length(Cliq),1000)]
    C1 <- length(Cliq)
    
    # from Cliq get GENE list
    GENE.list <- list()
    for (i in 1:length(Cliq))
    {GENE.list <- c(GENE.list,list(rownames(as.matrix(Cliq[[i]]))))}
    
    
    if(methdod == 'kde')
    {
      clique.dentisy <- 0
      # 计算任意两个基因的密度
      for (i in 1:length(GENE.list))
      {
        GL <- GENE.list[[i]]
        if (length(GL)==1)
        {prob <- 1
        }else{
          prob = c()
          for(Gi1 in 1:(length(GL)-1))
          {
            for(Gii1 in (Gi1+1):length(GL))
            {
              Gi <- GL[Gi1]
              Gii <- GL[Gii1]
              prob <- c(prob,prob.pro[Gi1,Gii1])
            }
          }
          prob <- log(1+prob/C1)
          prob <- mean(prob)
        } 
        clique.dentisy <- clique.dentisy+prob
        prob.list <- c(prob.list,prob)
        
      }
    }
    
  }else{clique.dentisy <- Inf
  prob.list = Inf}
  return(list(clique.dentisy = exp(-1/clique.dentisy),prob.list = prob.list))
}
