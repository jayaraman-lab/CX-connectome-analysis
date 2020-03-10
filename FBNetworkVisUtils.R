
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
                       strength=1,
                       arrow = arrow(length = unit(1, "cm")),
                       end_cap = circle(1, 'cm')) + 
    geom_edge_loop(aes(direction=45,span=90,width=weightRelative,color=superType,strength=0.1),alpha=0.5) +
    geom_node_point(aes(color=superType),size=15) + 
    sTScale_edge +
    sTScale +
    geom_node_text(aes(label=name),angle=40,size=12,repel = textRepel) +
    theme_classic() + theme(legend.text=element_text(size=36),legend.title=element_text(size=36),
                            axis.line=element_blank(),axis.text.x=element_blank(),
                            axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),axis.title.y=element_blank(),) + 
    coord_fixed(ratio = 1,clip="off")
  
  return(gg)
}

# Function to group Delta synapses into quadrants
D12DatQuads <- function(D12BodyIds){
  
  # Get the synapse info for the Delta 12
  D12Dat <- neuprint_get_synapses(D12BodyIds)
  
  # Get type info
  D12Dat <- D12Dat %>% mutate(type = neuprint_get_meta(bodyid)$type, partnerType = neuprint_get_meta(partner)$type) %>%
    mutate(x=as.numeric(x),y=as.numeric(y),z=as.numeric(z),prepost=as.logical(prepost))
  D12Type <- unique(D12Dat$type)
  
  # Divide into quadrants
  kmeanMat = D12Dat %>% select(x,z)
  if (!(D12Type %in% c("Delta0M","Delta0N"))){
    if (D12Type == "Delta0K"){
      synClusts = kmeans(kmeanMat,3, iter.max = 50, nstart = 10)
    } else {
      synClusts = kmeans(kmeanMat,4, iter.max = 50, nstart = 10)
    }
    
    cents = as.data.frame(synClusts$centers)
    cents$cluster <- rownames(cents)
    cents$RL <- "L"
    cents[which(cents$x >= sort(cents$x, decreasing=T)[2], arr.ind=TRUE),]$RL <- "R"
    cents$UD <- "D"
    cents[which(cents$z < sort(cents$z, decreasing=T)[2], arr.ind=TRUE),]$UD <- "U"
    cents <- cents %>% mutate(quad = paste0(UD,RL))
    
    D12Dat$quad <- lapply(synClusts$cluster, function(x) cents$quad[match(x, cents$cluster)]) %>% unlist()
    
    
  } else {
    synClustsRL = kmeans(kmeanMat,2, iter.max = 50, nstart = 10)
    centsRL = as.data.frame(synClustsRL$centers)
    centsRL$cluster <- rownames(centsRL)
    centsRL$RL <- "L"
    centsRL[which.max(centsRL$x),]$RL <- "R"
    
    D12Dat$RL <- lapply(synClustsRL$cluster, function(x) centsRL$RL[match(x, centsRL$cluster)])
    D12Dat$UD <- ""
    
    for (clust in 1:2){
      
      kmeanMat = D12Dat %>% filter(synClustsRL$cluster == clust) %>% select(x,y,z)
      synClustsUD = kmeans(kmeanMat,3, iter.max = 250, nstart = 200)
      
      clustIDs <- synClustsUD$cluster
      downClustInds <- which(clustIDs == as.numeric(which.max(synClustsUD$centers[,3])))
      
      UD <- paste0(vector(mode = "character",length=length(clustIDs)),"U")
      UD[downClustInds] <- "D" 
      
      D12Dat[which(synClustsRL$cluster == clust),]$UD <- UD
    }
    
    D12Dat <- D12Dat %>% mutate(quad = paste0(UD,RL))
    D12Dat <- D12Dat[ , !(names(D12Dat) %in% c("UD","RL"))]
  }
  
  return(D12Dat)
}

# Function to create bar plots of the D12 quadrants

D12QuadBarPlot <- function(D12Dat){
  quads <- c("UL","UR","DL","DR")
  synThresh <- 3
  
  sTs = as.factor(c("FBt","D0","D6","PFL","PFN","PFR","FC","FR","FS"))
  sTScale <- scale_fill_discrete(drop=TRUE,limits = levels(sTs))
  
  gList <- list()
  for (q in 1:length(quads)){
    D12QuadPost <- D12Dat %>% filter(quad == quads[q]) %>% filter(prepost == FALSE) %>% group_by(partnerType) %>% summarize(quadWeight = length(partner))
    D12QuadPost <- D12QuadPost %>% filter(!is.na(partnerType)) %>% filter(quadWeight > synThresh)
    D12QuadPost$supertype <- D12QuadPost$partnerType %>% supertype()
    gPost <- ggplot(D12QuadPost) + geom_bar(aes(x = as.factor(partnerType),weight = quadWeight,fill=as.factor(supertype))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("output") + sTScale
    gList[[2*q-1]] <- gPost
    
    D12QuadPre <- D12Dat %>% filter(quad == quads[q]) %>% filter(prepost == TRUE) %>% group_by(partnerType) %>% summarize(quadWeight = length(partner))
    D12QuadPre <- D12QuadPre %>% filter(!is.na(partnerType)) %>% filter(quadWeight > synThresh)
    D12QuadPre$supertype <- D12QuadPre$partnerType %>% supertype()
    gPre <- ggplot(D12QuadPre) + geom_bar(aes(x = as.factor(partnerType),weight = quadWeight,fill=as.factor(supertype))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("input") + sTScale
    gList[[2*q]] <- gPre
  }
  
  return(gList)
}

# Create a network diagram for the D12s
D12NetworkPlot <- function(D12Dat)
{
  quads <- c("UL","UR","DL","DR")
  quadXs <- c(0.75,-0.75,0.25,-0.25)
  quadYs <- c(-0.5,-0.5,0.5,0.5)
  synThresh <- 3
  
  sTs <- c("D0","D6","FBt","FS","FC","FR","SA") %>% as.factor()
  sTScale <- scale_colour_discrete(drop=TRUE,limits = levels(sTs))
  
  
  inputTypes <- D12Dat %>% filter(prepost == TRUE) %>%
    group_by(partnerType,quad) %>%
    summarize(quadWeight = length(partner)) %>%
    filter(quadWeight > synThresh) %>%
    select(partnerType) %>% unique() %>% unlist() %>% as.character()
  inputTypes[which(grepl("FB",inputTypes))] <- "FBt"
  inputTypes <-  unique(inputTypes)
  
  outputTypes <- D12Dat %>% filter(prepost == FALSE) %>%
    group_by(partnerType,quad) %>%
    summarize(quadWeight = length(partner)) %>%
    filter(quadWeight > synThresh) %>%
    select(partnerType) %>% unique() %>% unlist() %>% as.character()
  outputTypes[which(grepl("FB",outputTypes))] <- "FBt"
  outputTypes <- unique(outputTypes)
  
  nodes = data.frame(type = c(unique(D12Dat$quad),inputTypes,outputTypes),
                     x = c(sapply(unique(D12Dat$quad), function(f) quadXs[which(f == quads)]),
                           seq(-1,1,length.out=length(inputTypes)),
                           seq(-1,1,length.out=length(outputTypes))),
                     y = c(sapply(unique(D12Dat$quad), function(f) quadYs[which(f == quads)]),
                           (vector(mode="numeric",length = length(inputTypes))),
                           (vector(mode="numeric",length = length(outputTypes)) - 1)))
  nodes$superType <- nodes$type %>% as.character() %>% supertype() %>% as.factor()
  nodes$superType[1:4] <- NA
  
  edges <- data.frame(from=c(),to=c())
  for (q in 1:length(quads)){
    D12QuadInput <- D12Dat %>% filter(quad == quads[q]) %>%
      filter(prepost == TRUE) %>%
      group_by(partnerType) %>%
      summarize(quadWeight = length(partner))
    if (length(which(grepl("FB",D12QuadInput$partnerType)))>0)
      D12QuadInput[which(grepl("FB",D12QuadInput$partnerType)),]$partnerType <- "FBt"
    D12QuadInput <- D12QuadInput %>% group_by(partnerType) %>% summarize(quadWeight = sum(quadWeight)) %>%
      filter(!is.na(partnerType)) %>% filter(quadWeight > synThresh)
    D12QuadInput <- D12QuadInput %>%
      mutate(from = sapply(partnerType, function(f) length(quads)+which(f == inputTypes)),
             to = q)
    
    
    D12QuadOutput <- D12Dat %>% filter(quad == quads[q]) %>%
      filter(prepost == FALSE) %>% group_by(partnerType) %>% summarize(quadWeight = length(partner))
    if (length(which(grepl("FB",D12QuadOutput$partnerType)))>0)
      D12QuadOutput[which(grepl("FB",D12QuadOutput$partnerType)),]$partnerType <- "FBt"
    D12QuadOutput <- D12QuadOutput %>% group_by(partnerType) %>% summarize(quadWeight = sum(quadWeight)) %>%
      filter(!is.na(partnerType)) %>% filter(quadWeight > synThresh)
    D12QuadOutput <- D12QuadOutput %>% 
      mutate(to = sapply(partnerType, function(f) length(quads)+length(inputTypes) + which(f == outputTypes)),
             from = q)
    
    edges <- rbind(edges,
                   select(D12QuadInput,to,from,quadWeight),
                   select(D12QuadOutput,to,from,quadWeight))
  }
  
  
  # Plot the network
  graph <- tbl_graph(nodes,edges)
  
  # Generate an output pdf of the plots
  gg <-
    ggraph(graph,layout="manual",x=nodes$x,y=nodes$y) + 
    geom_edge_fan(aes(width=quadWeight),colour="grey",alpha=0.5,
                  strength=1,
                  arrow = arrow(length = unit(0.5, "cm")),
                  end_cap = circle(0.5, 'cm')) + 
    geom_edge_loop(colour="grey",aes(direction=10,span=10,width=quadWeight),alpha=0.5) +
    geom_node_point(aes(color=superType),size=5) + 
    sTScale +
    geom_node_text(aes(label=type),angle=40,size=6) +
    scale_y_reverse() + theme_classic() + theme(legend.text=element_text(size=36)) +
    theme(legend.text=element_text(size=12),legend.title=element_text(size=12),
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank()) + 
    coord_fixed(ratio = 1) 
  
  return(gg)
}
