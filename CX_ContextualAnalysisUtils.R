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
plotOutputCorrMat <- function(){
  # Get the FBt downstream partners
  FBTypes <- neuprint_search("FB.*")$type %>% unique()
  FBTypes <- FBTypes[which(!is.na(FBTypes))]
  # Pull the connection table of postsynaptic partners for the FB tangential cells
  FBConnTab <- getConnectionTable(getTypesTable(FBTypes)$bodyid,"POST","FB")
  # Convert the connection table to a type to type table. type.from: FB tangential. type.to: their postsynaptic partners in the FB
  FBCT_T2T <- getTypeToTypeTable(FBConnTab)
  
  # Cluster the FB tangential cells by common inputs
  library(reshape2)
  # Select only the relevant columns
  Data4Clust <- FBCT_T2T %>% select(type.from,type.to,weightRelative)
  # Recast it into a matrix
  Data4Clust <- dcast(Data4Clust,type.from~type.to) # converts into a matrix. rows are type.from, columns are type.to (first column is the type.from names), values are weightRelative
  Data4Clust[is.na(Data4Clust)] <- 0 # populate non-connection with 0
  rownames(Data4Clust) <- Data4Clust$type.from # set type.from names as rownames
  Data4Clust <- Data4Clust %>% select(tail(colnames(Data4Clust),ncol(Data4Clust)-1)) # chop off the first column
  
  # Calculate the dissimilarity matrix
  Data4Clust <- scale(Data4Clust) # normalize and center the columns
  d <- dist(t(Data4Clust), method = "euclidean") # compute distance matrix across the rows of the input matrix, so need to transpose the original matrix to compute aross its columns
  # Perform hierarchical clustering
  hc <- hclust(d, method = "ward.D2" ) # can play around a little
  
  # Order the FB tangential cells' targets according to their clustering
  FBCT_T2T$type.to <- factor(FBCT_T2T$type.to, levels = hc$labels[hc$order]) #factorize type.to according to the hc$order; also rearrange the rows of FBCT_T2T?
  # Pull out the reordered matrix
  Data4Corr <- dcast(FBCT_T2T %>% select(type.from,type.to,weightRelative),type.from~type.to)
  Data4Corr[is.na(Data4Corr)] <- 0
  rownames(Data4Corr) <- Data4Corr$type.from
  Data4Corr <- Data4Corr %>% select(tail(colnames(Data4Corr),ncol(Data4Corr)-1)) # chop off the first column
  
  # Plot the correlation matrix
  library(corrplot)
  corrplot(cor(Data4Corr)) # plot the correlation matrix grouped by columns of the input matrix
}

