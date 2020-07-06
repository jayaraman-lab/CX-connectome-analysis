
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
  as.dist((1-connCor)/2)
}

bin_dist <- function(mat,threshold=0.01){
  dist(mat>threshold,method="binary")
}

plot_dist <- function(dd,order=TRUE){
  ddM <- as.matrix(dd)
  if (is.logical(order)){
    if (order){
      hcl <- hclust(dd)
      ddM <- ddM[hcl$order,hcl$order]
    }
    }else{
    ddM <- ddM[order,order]
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

df_conn_mat <- function(connMat,rowOrd=1:nrow(connMat),colOrd=1:ncol(connMat),thresh=0){
  connMat <- connMat[rowOrd,colOrd]
  connMat[connMat<=thresh] <- NA
  reshape2::melt(connMat,na.rm=TRUE)
}

path2weight <- function(nBag,to,from,ROI,stat="weightRelative"){
  nBag$outputs <- filter(nBag$outputs,roi==ROI)
  nBag$inputs <- filter(nBag$inputs,roi==ROI)
  dplyr::bind_rows(lapply(to,function(dest){
    step2conn <- filter(nBag$outputs,type.to == dest)
    dplyr::bind_rows(lapply(from,function(orig){
      step1conn <- filter(nBag$inputs,(type.from == orig) & (type.to %in% step2conn$type.from))
      step2connL <- filter(step2conn,type.from %in% step1conn$type.to)
      out <- data.frame(from=step1conn$type.from,
                        inter_type=step1conn$type.to,
                        to=step2conn$type.to[match(step1conn$type.to,step2conn$type.from)],
                        step1=step1conn[[stat]],
                        step2=step2conn[[stat]][match(step1conn$type.to,step2conn$type.from)])
    out[[stat]] <- out$step1 * out$step2
    out}))
  }))
}

fullBagPathWeight <- function(nBag,ROI,stat="weightRelative",from=NULL,to=NULL,select=TRUE){
  nBag$outputs <- filter(nBag$outputs,roi==ROI)
  nBag$inputs <- filter(nBag$inputs,roi==ROI)
  if (!is.null(from)){nBag$inputs <- filter(nBag$inputs,type.from %in% from)}
  if (!is.null(to)){nBag$outputs <- filter(nBag$outputs,type.to %in% to)}
  
  res <- inner_join(nBag$inputs,nBag$outputs,by=c("type.to"="type.from",
                                          "databaseType.to"="databaseType.from",
                                          "previous.type.to"="previous.type.from",
                                          "supertype.to1"="supertype.from1",
                                          "supertype.to2"="supertype.from2",
                                          "supertype.to3"="supertype.from3",
                                          "roi"="roi"),suffix=c(".inp",".out")) %>% rename(
                                            "type.inter"="type.to",
                                            "type.to"="type.to.out",
                                            "databaseType.inter"="databaseType.to",
                                            "databaseType.to"="databaseType.to.out",
                                            "previous.type.inter"="previous.type.to",
                                            "previous.type.to"="previous.type.to.out",
                                            "supertype.inter1"="supertype.to1",
                                            "supertype.to1"="supertype.to1.out",
                                            "supertype.inter2"="supertype.to2",
                                            "supertype.to2"="supertype.to2.out",
                                            "supertype.inter3"="supertype.to3",
                                            "supertype.to3"="supertype.to3.out"
                                          )
  res[[paste0(stat,"_pathway")]] <- res[[paste0(stat,".inp")]] *  res[[paste0(stat,".out")]] 
  if (select){res <- select(res,contains("pathway") | contains(".from") | contains(".inter") | contains(".to") | contains(stat))}
  res
}


hClust_connMat <- function(connMat,rowSelect=rownames(connMat),dist){
  connMatS <- connMat[rownames(connMat) %in% rowSelect,]
  connCor <- cor(connMatS,method="spearman")
  connCor[is.na(connCor)] <- 1
  d <- as.dist(1-connCor)
  hc <- hclust(d)
  hc
}