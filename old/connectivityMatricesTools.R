

df_conn_mat <- function(connMat,rowOrd=1:nrow(connMat),colOrd=1:ncol(connMat),thresh=0){
  connMat <- connMat[rowOrd,colOrd]
  connMat[connMat<=thresh] <- NA
  reshape2::melt(connMat,na.rm=TRUE)
}


hClust_connMat <- function(connMat,rowSelect=rownames(connMat),distance=cos_dist){
  connMatS <- connMat[rownames(connMat) %in% rowSelect,]
  d <- distance(connMatS)
  hc <- hclust(d)
  hc
}

connectivity_clusters <- function(connTable,slctROIs,unit=c("type","neuron","supertype1","supertype2","supertype3","databaseType"),clusterSide=c("inputs","outputs"),distance=cos_dist,knownStats=FALSE,...){
  unit <-match.arg(unit)
  clusterSide <- match.arg(clusterSide)
  if (unit=="neuron") {unit <- ""}else{
    unit=paste0(unit,".")}
  from <- paste0(unit,"from")
  to <- paste0(unit,"to")
  if (clusterSide=="inputs"){stat <- ifelse(knownStats,"knownOutputContribution","OutputContribution")}else{
    stat <- ifelse(knownStats,"knownWeightRelative","weightRelative")
  }
  connMat <- connectivityMatrix(connTable,slctROIs = slctROIs,allToAll = FALSE,from = from,to=to,value = stat,ref=clusterSide)
  connDist <- distance(connMat)
  connHC <- hclust(connDist)
  connCl <- cutree(connHC,...)
  list(distance=connDist,hc=connHC,clusters=connCl)
}
