########################################################################################
# Functions to analyze tangential/contextual inputs and their effects on CX networks   #
########################################################################################

# From DanTE: Get edges and nodes and plot graph from a given connection matrix
graphConTab_old <- function(conTab,xyLookup,textRepel,guideOnOff){
  # Get the table of nodes (types)
  nodes = data.frame(name = unique(c(conTab$type.from,conTab$type.to)))
  nodes$superType <- nodes$name %>% as.character %>% supertype()
  
  # Position the nodes according to the lookup table
  nodes$x <- sapply(nodes$name, function(x) xyLookup$x[match(x,xyLookup$type)])
  nodes$y <- sapply(nodes$name, function(x) xyLookup$y[match(x,xyLookup$type)])
  
  # Assign colors to the supertypes
  sTs <- xyLookup$type %>% as.character() %>% supertype() %>% unique() %>% sort() %>% as.factor()
  if (guideOnOff){
    sTScale <- scale_colour_discrete(drop=TRUE,limits = levels(sTs))
  } else {
    sTScale <- scale_colour_discrete(drop=TRUE,limits = levels(sTs),guide = FALSE)
  }
  sTScale_edge <- scale_edge_colour_discrete(drop=TRUE,limits = levels(sTs),guide = FALSE)
  
  # Get the edges from the connection table
  edges <- conTab[which((conTab$type.from %in% nodes$name) & (conTab$type.to %in% nodes$name)),] %>%
    mutate(to = sapply(type.to, function(f) which(f == nodes$name)),
           from = sapply(type.from, function(f) which(f == nodes$name)))
  edges$superType <- edges$type.from %>% as.character %>% supertype()
  
  # Get the mean weights between types in the connection table
  edges_Mean <- unique(edges[,c('from','to','superType')])
  meanWs <- c()
  for (i in 1:nrow(edges_Mean)){
    meanWs <- append(meanWs,
                     mean(edges[which((edges$from == edges_Mean$from[i]) & (edges$to == edges_Mean$to[i])),]$weightRelative))
  }
  edges_Mean$weightRelative <- meanWs
  
  # Plot the network
  graph <- tbl_graph(nodes,edges_Mean)
  gg <-
    ggraph(graph,layout="manual",x=nodes$x,y=nodes$y) + 
    geom_edge_diagonal(aes(width=weightRelative,color=superType),alpha=0.25,
                       strength=1,
                       arrow = arrow(length = unit(1, "cm")),
                       end_cap = circle(1, 'cm')) + 
    geom_edge_loop(aes(direction=45,span=90,width=weightRelative,color=superType,strength=0.1),alpha=0.25) +
    geom_node_point(aes(color=superType),size=8) + 
    sTScale_edge +
    sTScale +
    geom_node_text(aes(label=name),angle=40,size=6,repel = textRepel) +
    theme_classic() + theme(legend.text=element_text(size=12),legend.title=element_text(size=12),
                            axis.line=element_blank(),axis.text.x=element_blank(),
                            axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),axis.title.y=element_blank(),) + 
    coord_fixed(ratio = 1,clip="off")
  
  return(gg)
}


# From DanTE: Get edges and nodes and plot graph from a given connection matrix
graphConTabPolyChrome <- function(conTab,xyLookup,textRepel,guideOnOff){
  # Get the table of nodes (types)
  nodes = data.frame(name = xyLookup$type)
  nodes$superType <- nodes$name %>% as.character %>% supertype()
  
  # Position the nodes according to the lookup table
  nodes$x <- sapply(nodes$name, function(x) xyLookup$x[match(x,xyLookup$type)])
  nodes$y <- sapply(nodes$name, function(x) xyLookup$y[match(x,xyLookup$type)])
  
  # Assign colors to the supertypes
  pal <- supertype2Palette()
  if (guideOnOff){
    sTScale <- scale_colour_manual(values = pal$pal, drop=TRUE,limits = pal$breaks)
  } else {
    sTScale <- scale_colour_manual(values = pal$pal, drop=TRUE,limits = pal$breaks,guide = FALSE)
  }
  sTScale_edge <- scale_edge_colour_manual(values = pal$pal, drop=TRUE,limits = pal$breaks,guide = FALSE)
  
  # Get the edges from the connection table
  edges <- conTab[which((conTab$type.from %in% nodes$name) & (conTab$type.to %in% nodes$name)),] %>%
    mutate(to = sapply(type.to, function(f) which(f == nodes$name)),
           from = sapply(type.from, function(f) which(f == nodes$name)))
  edges$superType <- edges$type.from %>% as.character %>% supertype()
  
  # Get the mean weights between types in the connection table
  edges_Mean <- unique(edges[,c('from','to','superType')])
  meanWs <- c()
  for (i in 1:nrow(edges_Mean)){
    meanWs <- append(meanWs,
                     mean(edges[which((edges$from == edges_Mean$from[i]) & (edges$to == edges_Mean$to[i])),]$weightRelative))
  }
  edges_Mean$weightRelative <- meanWs
  
  # Plot the network
  graph <- tbl_graph(nodes,edges_Mean)
  gg <-
    ggraph(graph,layout="manual",x=nodes$x,y=nodes$y) + 
    #geom_edge_diagonal(aes(width=weightRelative,color=superType),alpha=0.5,
    geom_edge_diagonal(aes(width=weightRelative),alpha=0.5,                   
                       strength=0.2,
                       arrow = arrow(length = unit(0.5, "cm")),
                       end_cap = circle(0.5, 'cm')) + 
    #geom_edge_loop(aes(direction=45,span=90,width=weightRelative,color=superType,strength=0.1),alpha=0.5) +
    geom_node_point(aes(color=superType),size=4) + 
    sTScale_edge +
    sTScale +
    #geom_node_text(aes(x = xyLookup$xTxt, y = xyLookup$yTxt, label=name),size=3,repel = textRepel) +
    geom_node_text(aes(label=name),size=3,repel = textRepel) +
    theme_classic() + theme(legend.text=element_text(size=6),legend.title=element_text(size=8),
                            axis.line=element_blank(),axis.text.x=element_blank(),
                            axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),axis.title.y=element_blank(),) + 
    coord_fixed(ratio = 1,clip="off")
  
  return(gg)
}

# Modified example code from DanTE
# Cluster connection matrix by the correlation of the weightRelative connectivity to targets and the inverse (clustering targets based on their inputs)
plotCorrMatCluster <- function(PlotDir,Type2TypeConnTab,Type2TypeConnTabName){
  library(reshape2)
  # Select only the relevant columns
  Data4Clust <- Type2TypeConnTab %>% select(type.from,type.to,weightRelative)
  
  # Recast it into a matrix
  Data4Clust <- dcast(Data4Clust,type.from~type.to) # rows are type.from, columns are type.to except the first column has the type.from names, values are weightRelative
  Data4Clust[is.na(Data4Clust)] <- 0
  rownames(Data4Clust) <- Data4Clust$type.from
  Data4Clust <- Data4Clust %>% select(tail(colnames(Data4Clust),ncol(Data4Clust)-1)) # chop off the first column
  
  # Plot the correlation matrix among type.to; Cluster and rearrange the Type2TypeConnTab
  Type2TypeConnTab_HCbyTo <- plotCorrClusterByCol(PlotDir,Type2TypeConnTab,Type2TypeConnTabName,Data4Clust,'type.to')

  # Plot the correlation matrix among type.from; Cluster and rearrange the Type2TypeConnTab
  Data4Clust <- t(Data4Clust) # transpose Data4Clust to plot correlations between type.from
  Type2TypeConnTab_HCbyFrom <- plotCorrClusterByCol(PlotDir,Type2TypeConnTab,Type2TypeConnTabName,Data4Clust,'type.from')
  
  return(list(Type2TypeConnTab_HCbyFrom,Type2TypeConnTab_HCbyTo))
}

plotCorrClusterByCol <- function(PlotDir,Type2TypeConnTab,Type2TypeConnTabName,Data4Clust,clusterBy){
  # Plot the correlation matrix among the clusterBy types
  library(corrplot)
  postscript(file = paste(PlotDir,"/",Type2TypeConnTabName,"_corrPlotBy",clusterBy,".eps",sep=""),width=24,height=24)
  corrplot(cor(Data4Clust), order="hclust", hclust.method="ward.D2", tl.cex=0.25, tl.col="black") # cor computes the correlation matrix between the columns of the input matrix
  dev.off()
  
  # Calculate the dissimilarity matrix based on the clusterBy types
  Data4ClustSc <- scale(Data4Clust) # centers/scales the columns of Data4Clust (ROMAIN: MIGHT INTRODUCE ARTIFACTS)
  d <- dist(t(Data4ClustSc), method = "euclidean") # compute the distance matrix between rows (clusterBy types)
  # Perform hierarchical clustering and plot the dendrogram
  hc <- hclust(d, method = "ward.D2" )
  dend1 <- as.dendrogram(hc)
  postscript(file = paste(PlotDir,"/",Type2TypeConnTabName,"_clusterBy",clusterBy,".eps",sep=""),width=36,height=24)
  par(cex=0.25)
  plot(dend1)
  dev.off()
  
  # Order the Type2TypeConnTab according to the clustering of the clusterBy types
  Type2TypeConnTab_hc <- Type2TypeConnTab
  Type2TypeConnTab_hc[,clusterBy] <- factor(Type2TypeConnTab_hc[[clusterBy]], levels = hc$labels[hc$order], ordered=TRUE)
  Type2TypeConnTab_hc <- Type2TypeConnTab_hc %>% arrange(.data[[clusterBy]])
  # Type2TypeConnTab_hc[,clusterBy] <- as.vector(Type2TypeConnTab_hc[,clusterBy])
  
  # Plot and save the Type2TypeConnTab_hc based on the connectionMeasure of "weightRelative"
  plotType2TypeConnTab_hc <- plotConnectivityMatrix(Type2TypeConnTab_hc,byGroup="type",connectionMeasure="weightRelative")
  print(plotType2TypeConnTab_hc)
  ggsave(paste(Type2TypeConnTabName,"_connMat_clusterBy",clusterBy,".eps",sep=""), plot=plotType2TypeConnTab_hc, device="eps", path=PlotDir, scale=1, 
         width=24, height=16, units="in", dpi=300, limitsize=FALSE)
  
  return(Type2TypeConnTab_hc)
}

# From Hannah: Function to plot a distance matrix
HHplot_dist <- function(dd,order=TRUE){
  ddM <- as.matrix(dd)
  if (order){
    hcl <- hclust(dd)
    ddM <- ddM[hcl$order,hcl$order]
  }
  ggplot(reshape2::melt(ddM)) + geom_tile(aes(x=Var1,y=Var2,fill=value)) +
    theme_classic() +
    scale_fill_gradient2(low="black", mid="grey", high="white", 
                         midpoint =0.5, limits=c(0,1)) +
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5)) + 
    xlab("") + ylab("")+coord_fixed() +xlab("") + ylab("")
}

# Modified from Hannah: Cluster and plot the cosine distance matrix from a connectivity table data frame
cosDistClusterPlot <- function(PlotDir,Type2TypeConnTab,Type2TypeConnTabName){
  Type2TypeConnTab_cosDistClusterByInp <- cosDistClusterPlotBySide(PlotDir,Type2TypeConnTab,Type2TypeConnTabName,"inputs")
  Type2TypeConnTab_cosDistClusterByOut <- cosDistClusterPlotBySide(PlotDir,Type2TypeConnTab,Type2TypeConnTabName,"outputs")
  return(list(Type2TypeConnTab_cosDistClusterByInp,Type2TypeConnTab_cosDistClusterByOut))
}

cosDistClusterPlotBySide <- function(PlotDir,Type2TypeConnTab,Type2TypeConnTabName,clusterBy){
  Type2TypeConnMatBySide <- connectivityMatrix(Type2TypeConnTab,unique(Type2TypeConnTab$roi),allToAll=FALSE,from="type.from",to="type.to",value="weightRelative",ref=clusterBy)
  Type2TypeConnMatBySide_CosDist <- cos_dist(Type2TypeConnMatBySide)
  Type2TypeConnMatBySide_CosDistPlot <- HHplot_dist(Type2TypeConnMatBySide_CosDist, order=TRUE)
  print(Type2TypeConnMatBySide_CosDistPlot)
  ggsave(paste0(Type2TypeConnTabName,"_cosDistClusterBy",clusterBy,".eps"), plot=Type2TypeConnMatBySide_CosDistPlot, device="eps", path=PlotDir, scale=1, 
         width=24, height=24, units="in", dpi=300, limitsize=FALSE)
  return(Type2TypeConnMatBySide_CosDist)
}

# From DanTE: Function to plot a dendrogram
dendPlot <- function(hc,rotate){
  library(ggdendro)
  # Pull of the data
  dend <- as.dendrogram(hc)
  dend_data <- dendro_data(dend, type = "rectangle")
  # Add a supertype category
  HClabels <- dend_data$labels
  HClabels$supertype <- HClabels$label %>% as.character() %>% supertype() %>% as.factor()
  p <- ggplot(dend_data$segments)
  if (rotate){
    p <- p + geom_segment(aes(x = y, y = x, xend = yend, yend = xend)) +
      geom_text(data = HClabels, aes(y, x, label = label,color=supertype),
                hjust = 1, size = 1)
  } else {
    p <- p + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_text(data = HClabels, aes(x, y, label = label,color=supertype),
                hjust = 1, angle = 90, size = 1)
  }
  p <- p + theme_classic() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line=element_blank())
  return(p)
}


