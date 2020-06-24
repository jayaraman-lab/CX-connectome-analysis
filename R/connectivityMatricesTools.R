
cos_dist <- function(mat){
  sim <- mat / sqrt(rowSums(mat * mat))
  sim <- sim %*% t(sim)
  D_sim <- as.dist(1 - sim)
}

sqrt_cos_dist <- function(mat){
  sim <- sqrt(mat) / sqrt(rowSums(mat))
  sim <- sim %*% t(sim)
  D_sim <- as.dist(1 - sim)
}

cor_dist <- function(mat){
  connCor <- cor(t(mat),method="spearman")
  d <- as.dist((1-connCor)/2)
}

bin_dist <- function(mat,threshold=0.01){
  connCor <- dist(mat>threshold,method="binary")
}

plot_dist <- function(dd,order=TRUE){
  ddM <- as.matrix(dd)
  if (order){
    hcl <- hclust(dd)
    ddM <- ddM[hcl$order,hcl$order]
  }
  ggplot(reshape2::melt(ddM)) + geom_tile(aes(x=Var1,y=Var2,fill=value)) + theme_minimal_grid()+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5)) + xlab("") + ylab("")+coord_fixed() +xlab("") + ylab("")
}

library(ggiraph)

plot_dist_interactive <- function(dd,order=TRUE){
  ddM <- as.matrix(dd)
  if (order){
    hcl <- hclust(dd)
    ddM <- ddM[hcl$order,hcl$order]
  }
  ggplot(reshape2::melt(ddM)) + geom_tile_interactive(aes(x=Var1,y=Var2,fill=value,
                                                          data_id=paste0(Var1,Var2),
                                                          tooltip=paste(Var1,"to",Var2,". \nDist=",value))) + theme_minimal_grid()+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5)) + xlab("") + ylab("")+coord_fixed() +xlab("") + ylab("")
}



hClust_connMat <- function(connMat,rowSelect=rownames(connMat),dist){
  connMatS <- connMat[rownames(connMat) %in% rowSelect,]
  connCor <- cor(connMatS,method="spearman")
  connCor[is.na(connCor)] <- 1
  d <- as.dist(1-connCor)
  hc <- hclust(d)
  hc
}