
process_result <- function(clique.dentisy)
{
  clique.dentisy <- as.data.frame(clique.dentisy)
  colnames(clique.dentisy) <- c('pathway','pvalue','obj')
  clique.dentisy[which(clique.dentisy$obj == Inf),'pvalue'] <- 1
  rownames(clique.dentisy) <- clique.dentisy$pathway
  clique.dentisy <- clique.dentisy[,-1]
  clique.dentisy.obj <- clique.dentisy$obj
  clique.dentisy.pvalue <- clique.dentisy$pvalue
  clique.dentisy.adjpvalue <- as.data.frame(p.adjust(clique.dentisy.pvalue,method = "BH"))
  rownames(clique.dentisy.adjpvalue) <- rownames(clique.dentisy)
  clique.dentisy.pvalue <- as.data.frame(clique.dentisy.pvalue)
  rownames(clique.dentisy.pvalue) <- rownames(clique.dentisy)
  clique.dentisy$adjpvalue <- clique.dentisy.adjpvalue[,1]
  return(list(clique.dentisy = clique.dentisy,
              adjpvalue = clique.dentisy.adjpvalue, 
              pvalue = clique.dentisy.pvalue, 
              obj = clique.dentisy.obj))
}

order_result <- function(clique.dentisy)
{
  clique.dentisy.obj <- clique.dentisy$obj
  clique.dentisy.pvalue <- clique.dentisy$pvalue 
  clique.dentisy.adjpvalue <- clique.dentisy$adjpvalue 
  clique.dentisy <- clique.dentisy$clique.dentisy
  
  clique.dentisy.uni <- sort(unique(clique.dentisy[,1]))
  clique.dentisy$order <- rep(0,dim(clique.dentisy)[1])
  ord <- 0
  for(i in clique.dentisy.uni)
  {
    la <- which(clique.dentisy[,1] == i)    
    obj.la <- clique.dentisy.obj[la]        
    la.uni <- sort(unique(obj.la))          
    ord1 <- 0
    for (la.uni1 in la.uni)
    {
      la.uni2 <- intersect(which(clique.dentisy.obj == la.uni1), y = la)        
      clique.dentisy$order[la.uni2] <- ord + ord1 + 1                          
      ord1 <- ord1 + length(la.uni2)
    }
    ord <- ord + length(which(clique.dentisy[,1]  == i))
  }
  return(clique.dentisy) 
}