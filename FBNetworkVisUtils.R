############################################################
# Functions to generate FB network diagrams
############################################################

#Create x and y coordinates for a manual network layout of neurons - FB columnar types
xyLookupTableCol <- function(){
  superTypes <- c("PFN","Delta0","SAF","Delta6","PFR","PFG","FC","PFL","^FR","FS")
  lyer <- c(1,2,2,2,2,2,2,3,3,3)
  lout <- c("h","c","v","c","v","v","c","h","h","h")
  
  xyLookup <- data.frame(type = c(), x = c(), y = c())
  for (st in 1:length(superTypes)){
    types <- neuprint_search(paste0(superTypes[st],".*"))$type %>% unique() %>% sort()
    xs <- numeric(length(types)) +
      which(superTypes[which(lyer == lyer[st])] == superTypes[st]) - 0.5*length(which(lyer == lyer[st]))-0.5
    ys <- numeric(length(types)) + lyer[st]
    if (lout[st] == "h"){
      xs <- xs + seq(-0.5,0.5,length.out=length(types))
    } else if (lout[st] == "v") {
      ys <- ys + seq(-0.25,0.25,length.out=length(types))
    } else if (lout[st] == "c") {
      angs <- seq(-pi,pi,length.out=length(types)+1)
      angs <- angs[1:(length(angs)-1)]
      xs <- xs + 0.5*sin(angs)
      ys <- ys + 0.5*cos(angs)
    }
    if (lyer[st] == 1)
      xs <- seq(-2,2,length.out=length(types))
    if (lyer[st] == 3)
      xs <- seq(-0.5,0.5,length.out=length(types)) +
      2*which(superTypes[which(lyer == lyer[st])] == superTypes[st]) - length(which(lyer == lyer[st]))-1
    typeLookup <- data.frame(type = types, x = xs, y = ys)
    xyLookup <- rbind(xyLookup,typeLookup)
  }
  return(xyLookup)
  #ggplot(xyLookup,aes(x=x,y=y)) + geom_point() + scale_y_reverse() + geom_text(aes(label=type),hjust=1, vjust=1,size=3)
}


#Create x and y coordinates for a manual network layout of neurons - "inner" neurons
xyLookupTableInner <- function(){
  nronGps <- list(c("PFN"),
                  c("Delta0","Delta6","PFR_a","PFG","FC"),
                  c("PFL","^FR","FS","PFR_b"))
  lout <- c("h","c","h")
  
  xyLookup <- data.frame(type = c(), x = c(), y = c())
  for (st in 1:length(nronGps)){
    types <-  c()
    for (tp in 1:length(nronGps[[st]])){
      typesNow <- neuprint_search(paste0(nronGps[[st]][tp],".*"))$type %>% unique() %>% sort()
      types <- types %>% append(typesNow)
    }
    
    xs <- numeric(length(types))
    ys <- numeric(length(types)) + st
    
    if (lout[st] == "h"){
      xs <- xs + seq(-0.5,0.5,length.out=length(types))
    } else if (lout[st] == "c") {
      angs <- seq(-pi,pi,length.out=length(types)+1)
      angs <- angs[1:(length(angs)-1)]
      xs <- xs + 0.5*sin(angs)
      ys <- ys + 0.5*cos(angs)
    }
    
    typeLookup <- data.frame(type = types, x = xs, y = ys)
    xyLookup <- rbind(xyLookup,typeLookup)
  }
  return(xyLookup)
  #ggplot(xyLookup,aes(x=x,y=y)) + geom_point() + scale_y_reverse() + geom_text(aes(label=type),hjust=1, vjust=1,size=3)
}

# "Clean up" weight matrices
# For each presynaptic type, find the weight matrix with the corresponding postsynaptic type. If the number of connections is less than the number of neurons in the presynaptic type (i.e. one of those neurons doesn't have a partner) or vice versa, remove those entries from the weight matrix

cleanUpConTab <- function(conTab){
  preTypes <- conTab$type.from %>% unique()
  
  for (pre in 1:length(preTypes)){
    numPres = nrow(neuprint_search(paste0(preTypes[pre],".*")))
    postTypes <- conTab[which(conTab$type.from == preTypes[pre]),]$type.to %>% unique()
    for (post in 1:length(postTypes)){
      numPost = nrow(neuprint_search(paste0(postTypes[post],".*")))
      conSub <- conTab[which((grepl(preTypes[pre],conTab$type.from)) & 
                               (grepl(postTypes[post],conTab$type.to))),]
      if ((nrow(conSub) < 0.5*numPres) | (nrow(conSub) < 0.5*numPost)){
        conTab <- conTab[which(!(grepl(preTypes[pre],conTab$type.from)) | 
                                 !(grepl(postTypes[post],conTab$type.to))),]
      }
    }
  }
  return(conTab)
}

# Get edges and nodes and plot graph from a given connection matrix
graphConTab <- function(conTab,xyLookup,textRepel,guideOnOff){
  
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
    geom_edge_diagonal(aes(width=weightRelative,color=superType),alpha=0.5,
                       strength=0.5,
                       arrow = arrow(length = unit(1, "cm")),
                       end_cap = circle(1, 'cm')) + 
    geom_edge_loop(aes(direction=45,span=90,width=weightRelative,color=superType,strength=0.1),alpha=0.5) +
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

############################################################
# Functions for Delta 0/12 analyses
############################################################

# Cluster the synapse of an individual D12 or D0 neuron into up and down, right and left
DeltaSynClust <- function(bodyids, numClusts){
  
  xCent <- 26000
  zMax <- 24000
  
  # Get the synapse data
  synDat <- neuprint_get_synapses(bodyids)
  synDat <- synDat %>% mutate(type = neuprint_get_meta(bodyid)$type, partnerType = neuprint_get_meta(partner)$type) %>%
    mutate(x=as.numeric(x),y=as.numeric(y),z=as.numeric(z),prepost=as.logical(prepost))
  synDat$quad <- ""
  
  zRad <- zMax - mean(synDat$z)

  for (bid in 1:length(bodyids)){
    
    # Get the x and z positions of the synapses for a given bodyid
    kmeanMat = synDat %>% filter(bodyid == bodyids[bid]) %>% select(x,z)
    
    # Cluster them
    synClusts = kmeans(kmeanMat,numClusts, iter.max = 50, nstart = 10)
    
    # Find the center position and cluster into 
    cents = as.data.frame(synClusts$centers)
    cents$cluster <- rownames(cents)
    
    cents$RL <- "L"
    if (length(which(cents$x > xCent)) > 0)
      cents[which(cents$x > xCent),]$RL <- "R"
    cents$UD <- "U"
    if (length(which(sqrt(0.8*(cents$x-xCent)^2+(zMax-cents$z)^2) < zRad)) > 0)
      cents[which(sqrt(0.8*(cents$x-xCent)^2+(zMax-cents$z)^2) < zRad),]$UD <- "D"
    
    cents <- cents %>% mutate(quad = paste0(UD,RL))
    
    synDat[which(synDat$bodyid == bodyids[bid]),]$quad <- lapply(synClusts$cluster, function(x) cents$quad[match(x, cents$cluster)]) %>% unlist()
    
  }
  
  return(synDat)
}

# Convert a connection table to a full matrix, generate a correlation matrix across inputs or outputs, and plot it
corMatPlot <- function(connMat,preId,postId){
  
  
  corMat <- corMat(connMat,preId,postId)
  
  corMatPlot <- ggplot(corMat) + 
    theme_classic() + theme(axis.text.x = element_text(angle = 90))
  corMatPlot <- corMatPlot + geom_tile(aes(Var1,Var2,fill=value)) + 
    scale_fill_gradient(low="white", high="red") + coord_fixed(ratio=1)
  
  return(corMatPlot)
}

# Convert a connection table to a full matrix, generate a correlation matrix
corMat <- function(connMat,preId,postId){
  preIds <- connMat[preId] %>% unlist() %>% unique() %>% as.character() %>% sort()
  postIds <- connMat[postId] %>% unlist() %>% unique() %>% as.character() %>% sort()
  
  justWeights = matrix(0L,nrow=length(preIds),ncol=length(postIds))
  
  for (i in 1:length(preIds)){
    for (j in 1:length(postIds)){
      w = connMat[which(connMat[preId] == preIds[i] & connMat[postId] == postIds[j]),]$weightMean
      if (length(w) > 0)
        justWeights[i,j] = w
    }
  }
  rownames(justWeights) <- preIds
  colnames(justWeights) <- postIds
  
  library(reshape2)
  
  corMat <- cor(justWeights) %>% melt()
  #ggplot(melt(justWeights), aes(Var1,Var2, fill=value)) + geom_raster()
  
  return(corMat)
}