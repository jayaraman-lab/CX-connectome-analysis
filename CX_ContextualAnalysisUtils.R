########################################################################################
# Functions to analyze tangential/contextual inputs and their effects on CX networks   #
########################################################################################

# Modified from supertype2Palette() function in Romain's neuprintrExtra package
CXplusPalette <- function(){
  s2 <- c("vDelta","v\u0394","hDelta","h\u0394","Delta7","\u03947","EL","EPG","EPGt","ExR","FBt","FC","FR","FS","LNO","SPS-PB","LPsP","P","PEG","PEN","PFGs","PFL","PFN","PFR","ER","SA","MBON","CRE","SIP","SMP","SLP","LAL","IB","LH","LHP","LHPV","LHPD","ATL")
  pal <- paletteer::paletteer_d("Polychrome::palette36")[c(35,35,32,32,28,28,8,12,33,6,10,9,3,25,18,21,30,31,34,16,27,7,26,1,15,36,4,5,13,17,29,14,22,11,11,11,11,23)]
  names(pal) <- s2
  list(pal=pal,breaks=s2)
}

# Read tables of PN to MBON with input modality percentages and return a data frame
getPN2MBONinputTable <- function(){
  PN2MBON_InputPc <- read.csv('/Users/danc/Dropbox (HHMI)/WorkDrop/Data/FIBSEM_Analysis/MB/PN2MBON.csv', header = FALSE)
  PN2ATMBON_InputPc <- read.csv('/Users/danc/Dropbox (HHMI)/WorkDrop/Data/FIBSEM_Analysis/MB/PN2ATMBON.csv', header = FALSE)
  PNtypes <- c('UniGlom', 'MultiGlom', 'ThermoHygro', 'Visual')
  colnames(PN2MBON_InputPc) <- PNtypes
  colnames(PN2ATMBON_InputPc) <- PNtypes
  MBONtypeNames <- c('MBON01', 'MBON02', 'MBON03', 'MBON04', 'MBON05', 'MBON06', 'MBON07', 'MBON09', 'MBON11', 'MBON12', 'MBON13', 'MBON14', 'MBON15', 'MBON16', 'MBON17', 'MBON18', 'MBON19a', 'MBON19b', 'MBON20', 'MBON21', 'MBON22', 'MBON23')
  atypicalMBONtypeNames <- c('MBON10', 'MBON24', 'MBON25', 'MBON26', 'MBON27', 'MBON29', 'MBON30', 'MBON31', 'MBON32', 'MBON33', 'MBON34', 'MBON35', 'MBON28')
  allMBONtypeNames <- c(MBONtypeNames,atypicalMBONtypeNames)
  rownames(PN2MBON_InputPc) <- MBONtypeNames
  rownames(PN2ATMBON_InputPc) <- atypicalMBONtypeNames
  PN2allMBON_InputPct <- rbind(PN2MBON_InputPc,PN2ATMBON_InputPc)
  PN2allMBON_InputPct <- PN2allMBON_InputPct %>% mutate(Olfactory = UniGlom + MultiGlom) # add a column for the sum of UniGlom and MultiGlom as Olfactory
  PN2allMBON_InputPct <- PN2allMBON_InputPct %>% mutate(MBON = allMBONtypeNames) # add a column for matching MBON type names
  
  PN2allMBON_InputPct[[17,"MBON"]] <- "MBON19" # combine the 2 MBON19 types
  PN2allMBON_InputPct[[18,"MBON"]] <- "MBON19"
  
  return(PN2allMBON_InputPct)
}

# Filter and analyze the MBON to CX direct and indirect tables based on input modality
filterConnTabsByInputMod <- function(PlotDir,inputModTab,inputMod,inputThresh,filtCol,directConnTab,indirectConnTab1,indirectConnTab2){
  
  inputModTabFiltered <- inputModTab %>% filter(!!sym(inputMod) > inputThresh) # filter inputModTab based on inputMod and inputThresh
  directConnTabFiltered <- directConnTab %>% filter(databaseType.from %in% inputModTabFiltered[[filtCol]]) # filter directConnTab based on filtered inputModTab
  indirectConnTab1Filtered <- indirectConnTab1 %>% filter(databaseType.from %in% inputModTabFiltered[[filtCol]]) # filter indirectConnTab1 based on filtered inputModTab
  indirectConnTab2Filtered <- indirectConnTab2 %>% filter(type.from %in% indirectConnTab1Filtered$type.to) # filter indirectConnTab2 based on filtered indirectConnTab1
  
  # Cluster by cosine distance and plot
  directConnTabFiltered_CosDist <- list()
  if (nrow(directConnTabFiltered)>1){
    directConnTabFiltered_CosDist <- cosDistClusterPlot(PlotDir,directConnTabFiltered,paste0(inputMod,filtCol,"directConnTab_Type2Type"))
    directConnTabFiltered <- directConnTabFiltered_CosDist[[3]]
  }
  
  indirectConnTab1Filtered_CosDist <- cosDistClusterPlot(PlotDir,indirectConnTab1Filtered,paste0(inputMod,filtCol,"indirectConnTab1_Type2Type"))
  indirectConnTab1Filtered <- indirectConnTab1Filtered_CosDist[[3]]
  indirectConnTab2Filtered_CosDist <- cosDistClusterPlot(PlotDir,indirectConnTab2Filtered,paste0(inputMod,filtCol,"indirectConnTab2_Type2Type"))
  indirectConnTab2Filtered <- indirectConnTab2Filtered_CosDist[[3]]
  
  # Combine the filtered direct and indirect tables
  directConns <- directConnTabFiltered %>% select(type.from,type.to,weightRelative)
  indirectConns1 <- indirectConnTab1Filtered %>% select(type.from,type.to,weightRelative)
  indirectConns2 <- indirectConnTab2Filtered %>% select(type.from,type.to,weightRelative)
  drctIndrctComboTable <- bind_rows(directConns,indirectConns1,indirectConns2, .id = NULL)
  
  # Re-plot pathways
  fromTypes <- unique(c(unique(as.vector(directConns$type.from)),unique(as.vector(indirectConns1$type.from))))
  numFrom <- length(fromTypes)
  midNodes <- unique(as.vector(indirectConns2$type.from))
  numMidNodes <- length(midNodes)
  allTargets <- unique(c(unique(as.vector(directConns$type.to)),unique(as.vector(indirectConns2$type.to))))
  numAllTargets <- length(allTargets)
  types <- unique(c(fromTypes,midNodes,allTargets))
  numTypes <- length(types)
  
  xyLookup = data.frame(type = types, x = c(rep(-1,times = numFrom), rep(0,times = numMidNodes), rep(1,times = numAllTargets)), 
                        y = c(seq(-1,1,length.out = numFrom), seq(0,2,length.out = numMidNodes), seq(-1,1.5,length.out = numAllTargets)))
  
  drctIndrctComboPath <- graphConTabPolyChrome(drctIndrctComboTable,xyLookup,FALSE,TRUE) # graph the TypeToType ConnTable using the lookupTable
  drctIndrctComboPath <- drctIndrctComboPath + scale_y_reverse()
  print(drctIndrctComboPath)
  ggsave(paste0(inputMod,filtCol,"drctIndrctComboPath.svg"), plot=drctIndrctComboPath, device="svg", path=PlotDir, scale=1, width=30, height=90, units="in", dpi=300, limitsize=FALSE)
  
  # Cluster by cosine distance and plot
  directConnTabFiltered_CosDist <- list()
  if (nrow(directConnTabFiltered)>1){
    directConnTabFiltered_CosDist <- cosDistClusterPlot(PlotDir,directConnTabFiltered,paste0(inputMod,filtCol,"directConnTab_Type2Type"))
  }
  
  indirectConnTab1Filtered_CosDist <- cosDistClusterPlot(PlotDir,indirectConnTab1Filtered,paste0(inputMod,filtCol,"indirectConnTab1_Type2Type"))
  indirectConnTab2Filtered_CosDist <- cosDistClusterPlot(PlotDir,indirectConnTab2Filtered,paste0(inputMod,filtCol,"indirectConnTab2_Type2Type"))
  
  return(list(inputModTabFiltered,directConnTabFiltered,indirectConnTab1Filtered,indirectConnTab2Filtered,drctIndrctComboTable,directConnTabFiltered_CosDist,indirectConnTab1Filtered_CosDist,indirectConnTab2Filtered_CosDist))
}

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
                       strength=0.5,
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
  
  # Assign colors to the supertypes using the CXplusPalette() function for an extended palette
  pal <- CXplusPalette()
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
    geom_edge_diagonal(aes(width=weightRelative,color=superType),alpha=0.5,
    #geom_edge_diagonal(aes(width=weightRelative),alpha=0.5,                   
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

# Modified from Hannah: Cluster and plot the cosine distance matrix from a connectivity table data frame
cosDistClusterPlot <- function(PlotDir,Type2TypeConnTab,Type2TypeConnTabName,plotFacet=TRUE){
  # Apply cosine distance clustering separately to the from and to side
  Type2TypeConnTab_cosDistClusterByInp <- cosDistClusterPlotBySide(PlotDir,Type2TypeConnTab,Type2TypeConnTabName,"inputs")
  Type2TypeConnTab_cosDistClusterByOut <- cosDistClusterPlotBySide(PlotDir,Type2TypeConnTab_cosDistClusterByInp[[1]],Type2TypeConnTabName,"outputs")

  # Plot the factorized Type2TypeConnTab_hc based on the connectionMeasure of "weightRelative" and grouped by the factors
  Type2TypeConnTab_hc <- Type2TypeConnTab_cosDistClusterByOut[[1]]
  Type2TypeConnTab_hc <- Type2TypeConnTab_hc %>% arrange(cluster.from,type.from,cluster.to,type.to) # re-arrange the rows of Type2TypeConnTab_hc
  plotType2TypeConnTab_hc <- plotConnectivityMatrix(Type2TypeConnTab_hc,byGroup="type",connectionMeasure="weightRelative")

  # plotType2TypeConnTab_hc <- plotType2TypeConnTab_hc + scale_x_discrete(breaks=levels(Type2TypeConnTab_hc$type.to)) + scale_y_discrete(breaks=levels(Type2TypeConnTab_hc$type.from))
  
  # Facet the matrix
  if (plotFacet){
    # plotType2TypeConnTab_hc <- plotType2TypeConnTab_hc + facet_grid(as.formula("cluster.from ~ cluster.to"),scale="free",space="free")
    plotType2TypeConnTab_hc <- plotType2TypeConnTab_hc + facet_grid(rows=vars(cluster.from),cols=vars(cluster.to),scale="free",space="free") + 
      theme(axis.text.x = element_blank(),axis.text.y = element_blank(), 
            strip.placement = "outside", #strip.background = element_rect(fill=NA, colour="grey50"),
            #strip.text.y.left = element_text(angle = 0),
            #strip.text.x.bottom = element_text(angle = 90),
            strip.background = element_blank(), #remove background for facet labels
            panel.border = element_rect(colour = "grey", fill = NA, size=0.3), #add grey border
            panel.spacing = unit(0.05, "lines")) #space between facets
  }
  
  # Color the row and column labels
  yConnColors <- Type2TypeConnTab_cosDistClusterByInp[[6]]
  xConnColors <- Type2TypeConnTab_cosDistClusterByOut[[6]]
  plotType2TypeConnTab_hc <- plotType2TypeConnTab_hc + theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5,colour = xConnColors),axis.text.y = element_text(colour = yConnColors))

  print(plotType2TypeConnTab_hc)
  ggsave(paste0(Type2TypeConnTabName,"_cosDistClustConnMat",".eps"), plot=plotType2TypeConnTab_hc, device="eps", path=PlotDir, scale=1, 
         width=24, height=16, units="in", dpi=300, limitsize=FALSE)
  
  return(list(Type2TypeConnTab_cosDistClusterByInp,Type2TypeConnTab_cosDistClusterByOut,Type2TypeConnTab_hc))
}

cosDistClusterPlotBySide <- function(PlotDir,Type2TypeConnTab,Type2TypeConnTabName,InpOrOutp){
  Type2TypeConnMatBySide <- connectivityMatrix(Type2TypeConnTab,unique(Type2TypeConnTab$roi),allToAll=FALSE,from="type.from",to="type.to",value="weightRelative",ref=InpOrOutp)
  Type2TypeConnMatBySide_CosDist <- cos_dist(Type2TypeConnMatBySide)
  
  # Segment the clusters at cut = 0.8
  hcl <- hclust(Type2TypeConnMatBySide_CosDist)
  cut <- 0.8
  # tplt <- plot(hcl,hang= -.5,cex = 0.8) # + rect.hclust(hcl, h=cut, border = "red") # dendrogram
  clu.h <- cutree(hcl,h=cut) # cut tree/dendrogram from height 80
  typeClusters <- stack(clu.h) %>% rename(type=ind, clust=values) %>% arrange(clust) # convert the cluster vector assignment to a data frame
  
  # Reorganize Type2TypeConnTab based on the clustering
  Type2TypeConnTab_hc <- Type2TypeConnTab
  clusterBy <- ifelse(InpOrOutp=="inputs","type.from","type.to")
  clusterCol <- ifelse(InpOrOutp=="inputs","cluster.from","cluster.to")
  Type2TypeConnTab_hc[,clusterBy] <- factor(Type2TypeConnTab_hc[[clusterBy]], levels = hcl$labels[hcl$order], ordered=TRUE)
  Type2TypeConnTab_hc <- mutate(Type2TypeConnTab_hc, clusterVals = clu.h[as.vector(Type2TypeConnTab_hc[[clusterBy]])]) # add a column of cluster numbers matched to the clusterBy types as indexed by de-factorized type names
  colnames(Type2TypeConnTab_hc)[ncol(Type2TypeConnTab_hc)] <- clusterCol
  
  # Plot the cosine distance square matrix
  Type2TypeConnMatBySide_CosDistClustPlot <- HHplot_dist(Type2TypeConnMatBySide_CosDist, order=TRUE, colorAxis=TRUE, InpOrOutp)
  Type2TypeConnMatBySide_CosDistPlot <- Type2TypeConnMatBySide_CosDistClustPlot[[1]]
  Type2TypeConnMatBySide_CosDistClustColor <- Type2TypeConnMatBySide_CosDistClustPlot[[2]]
  print(Type2TypeConnMatBySide_CosDistPlot)
  ggsave(paste0(Type2TypeConnTabName,"_cosDistClusterBy",clusterBy,".eps"), plot=Type2TypeConnMatBySide_CosDistPlot, device="eps", path=PlotDir, scale=1, 
         width=24, height=24, units="in", dpi=300, limitsize=FALSE)
  
  Type2TypeConnTab_hc[,clusterBy] <- as.vector(Type2TypeConnTab_hc[[clusterBy]])
  return(list(Type2TypeConnTab_hc,Type2TypeConnMatBySide_CosDist,hcl,clu.h,typeClusters,Type2TypeConnMatBySide_CosDistClustColor))
}

# Modified from Hannah and neuprintrExtra: Function to plot a distance matrix
HHplot_dist <- function(dd,order=TRUE,colorAxis=TRUE,InpOrOutp,
                        axesPalette=paletteer::paletteer_d("Polychrome::palette36")[3:36],
                        theme=theme_classic()){
  ddM <- as.matrix(dd)
  if (order){
    hcl <- hclust(dd)
    ddM <- ddM[hcl$order,hcl$order]
    
    if (colorAxis){
      clusters <- cutree(hcl,h=0.8)
      if(InpOrOutp=="outputs"){axesPalette <- rev(axesPalette)} # assign color in the reversed order for outputs to avoid overlap with inputs
      connCols <- axesPalette[clusters[rownames(ddM)]]
      names(connCols) <- rownames(ddM)
    }
  }
  
  p <- ggplot(reshape2::melt(ddM)) + geom_tile(aes(x=Var1,y=Var2,fill=value)) + theme +
    scale_fill_gradient2(low="black", mid="grey", high="white", midpoint=0.5, limits=c(0,1)) +
    theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + xlab("") + ylab("") + coord_fixed()
  
  if (colorAxis){
    p <- p + theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,colour=connCols),axis.text.y = element_text(colour=connCols))
  }
  
  return(list(p,connCols))
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


