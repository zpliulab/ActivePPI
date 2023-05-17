make_data <- function(nl = 0.01, show = FALSE,saved = FALSE)
{
# nl: noise radio (0-1) 
  
library(igraph)
library(ks)
library(ggplot2)

#--------------------------   load network
adj <- read.csv('./simu/data/simexp_new.csv',header = F)
rownames(adj) <- paste0('P',1:10)
colnames(adj) <- paste0('P',1:10)
adj.prob <- adj/10

for(i in 1:10)
{
  adj.prob[i,i] <- 1 - sum(adj.prob[i,])
}

adj[adj>0] <- 1

adj.list <- data.frame(matrix(ncol = 3, nrow = 1))
d = 1
for (i1 in 1:9)
{
  for (j1 in (i1+1):10)
  {
    if (adj[rownames(adj)[i1],rownames(adj)[j1]])
    {
      adj.list[d,] <- cbind(rownames(adj)[i1],rownames(adj)[j1],adj.prob[rownames(adj)[i1],rownames(adj)[j1]]*10)
      d <- d+1
    }
  }
}
colnames(adj.list) <- c('from','to','value')

#openxlsx::write.xlsx(adj.list,file = './result/simu_data_list.xlsx')

col = circlize::colorRamp2(breaks = c(0,0.02, 0.05, 0.09, 1), colors = c('red','#FF9D9D', '#FFF3F3','#A1CEFF','#2937B3'))

if (show)
{
  pdf('./result/simu_adj.pdf',width=5.5, height=5)
  corrplot::corrplot(as.matrix(adj.prob), method = "number", type = "lower",
                     is.corr = F,diag=T,
                     col="black",tl.pos = "lt",
                     tl.col = "black", tl.cex = 1, tl.srt = 45)
  corrplot::corrplot(as.matrix(adj.prob), method = "circle", type = "upper",col.lim=c(0,1),
                     is.corr = F,tl.pos = "lt",diag=T,add = T,
                     col=colorRampPalette(c('red','#FF9D9D', '#FFF3F3','#A1CEFF','#2937B3'))(100),
                     tl.col = "black", tl.cex = 1, tl.srt = 45)
  dev.off()
}

gra <- graph_from_adjacency_matrix( as.matrix(adj),
                                    mode = "undirected",
                                    weighted = NULL,
                                    diag = TRUE,
                                    add.colnames = NULL,
                                    add.rownames = NA)

#--------------------------   update data
G.exp.n.sum <- c()
Ntime <- 50 # iterations
for (i in 1:30)
{
  G.exp.n.sum2 <- c()
  G.exp.n <- as.matrix(rnorm(10,mean=0,sd=1))
  rownames(G.exp.n) <- paste0('G',1:10)
  G.exp.n <- t(as.matrix(G.exp.n))
  for(i in 1:Ntime)
  {
    G.exp <- G.exp.n  %*%  as.matrix(adj.prob)
    G.exp <- G.exp + rnorm(10,mean=0,sd=nl)
    G.exp.n.sum2 <- cbind(G.exp.n.sum2,t(G.exp.n))
    G.exp.n <- G.exp
  }
  G.exp.n.sum2 <-G.exp.n.sum2[,Ntime]
  G.exp.n.sum <- cbind(G.exp.n.sum,G.exp.n.sum2)
}
colnames(G.exp.n.sum) <- paste0('Sample',1:dim(G.exp.n.sum)[2])

#------------------------------- save data
simu_data<-c()
simu_data$exp <- G.exp.n.sum
simu_data$gra <- gra
if(saved)
{save(simu_data,file= paste0('./data/simu_data_',nl,'.Rdata'))}

#------------------------------- show data
if(show)
{ G.exp.n.sum2 <- c()
  G.exp.n <- as.matrix(rnorm(10,mean=0,sd=1))
  rownames(G.exp.n) <- paste0('G',1:10)
  G.exp.n <- t(as.matrix(G.exp.n))
  for(i in 1:Ntime)
  {
    G.exp <- G.exp.n  %*%  as.matrix(adj.prob)
    G.exp <- G.exp + rnorm(10,mean=0,sd=nl)
    G.exp.n.sum2 <- cbind(G.exp.n.sum2,t(G.exp.n))
    G.exp.n <- G.exp
  }
  colnames(G.exp.n.sum2) <- paste0('T',1:Ntime)
  rname <- as.matrix(rownames(G.exp.n.sum2))
  G.exp.n.sum2 <- cbind(as.data.frame(G.exp.n.sum2),rname)
  G.exp.n.sum2 <- reshape2::melt(G.exp.n.sum2, id = 'rname')
  colnames(G.exp.n.sum2)  <-  c("Protein","Time",'Value')
  pdf('./result/convergence.pdf',width=8, height=7)
  ggplot(G.exp.n.sum2,aes(x = Time,y = Value,group=Protein,color=Protein))+
    geom_point(size=2)+
    geom_line(position = position_dodge(0.1),cex=1.2)+
    ylab('Expression value')+
    ggsci::scale_color_npg(breaks = paste0('G',1:10))+
    theme_classic()+
    theme(axis.text.x = element_text(size = 12, angle = 90,vjust = 0.5,colour = "black"),
          axis.text.y = element_text(size = 12,colour = "black"),
          axis.title.y = element_text(size = 12,colour = "black"), 
          legend.text = element_text(size = 12,colour = "black"), 
          plot.background = element_rect( fill = "white"),
          panel.background = element_rect( fill = "white"),
          axis.line = element_line(colour = "black"))
  dev.off()
  pdf('./result/simu_exp.pdf',width=10, height=4)
  pheatmap::pheatmap(simu_data$exp,
                     cluster_cols = F, cluster_rows = F, scale = "none",
                     treeheight_col = 0, treeheight_row = 0,
                     display_numbers = F,
                     border_color = "black",color = colorRampPalette(c('red', 'white', 'blue'))(100))
  dev.off()
}
return(simu_data)
}