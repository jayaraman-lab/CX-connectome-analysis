########################################################################################################################
# Visualization tools for connectivity tables
########################################################################################################################


### Connectivity matrix

plotConnectivityMatrix = function(myConTab, byGroup = "name", connectionMeasure="weightRelative", cmax=NULL) {
  #' Generate plot of connectivity matrix using ggplot
  #' @param myConTab Neuprint connection table
  #' @param synapseCutOff Minimum number of synapses between two partners to be considered a connection
  #' @param byType choose IDs or types
  #' @param connectionMeasure choose which weigh
  #' @return conmatPlot is a ggplot object
  # cool: low = "papayawhip", mid = "darkseagreen2", high = "steelblue4"
  # warm: 
  
  myConTab = myConTab %>% ungroup() %>% mutate(plotWeight = myConTab[[connectionMeasure]])
  
  if (byGroup == "id"){
    myConTab$nameid.from = paste(as.character(myConTab$name.from), as.character(myConTab$from), sep = "_")
    myConTab$nameid.to = paste(as.character(myConTab$name.to), as.character(myConTab$to), sep = "_")
  }
  if (is.null(cmax)){
    cmax = max(myConTab$plotWeight)
  }
  
  conmatPlot = ggplot(myConTab) + 
    theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="thistle", mid="blueviolet", high="black", 
                         midpoint =0.5*cmax, limits=c(0,cmax))  
  if (byGroup == "name"){
    conmatPlot =  conmatPlot + geom_tile(aes(name.to,name.from,fill=plotWeight))
  } else if(byGroup == "id"){
    conmatPlot =  conmatPlot + geom_tile(aes(nameid.to,nameid.from,fill=plotWeight))
  } else if(byGroup == "type"){
    conmatPlot =  conmatPlot + geom_tile(aes(type.to,type.from,fill=plotWeight))
  } else if(byGroup == "Glom_to_Col"){
    conmatPlot =  conmatPlot + geom_tile(aes(FBcol.to,PBglom.from,fill=plotWeight))
  } else if(byGroup == "Col_to_Glom"){
    conmatPlot =  conmatPlot + geom_tile(aes(PBglom.to,FBcol.from,fill=plotWeight))
  } else{
    conmatPlot =  conmatPlot + geom_tile(aes(nameid.to,nameid.from,fill=plotWeight))
  }
  
  return(conmatPlot)
}

addMatrixPlotLabs = function(conmatPlot, preName, postName, slctROI, connectionMeasure="weightRelative") {
  #' Add labels to a connection matrix plot
  #' @param preName String describing presynaptic neurons
  #' @param postName String describing postsynaptic neurons
  #' @param slctROI Bane of ROI used to construct connection table
  #' @return conmatPlot is a ggplot object
  conmatPlot = conmatPlot +
    xlab("Post-synaptic neuron")+ylab("Pre-synaptic neuron") + labs(fill=connectionMeasure) +
    labs(title=paste(preName,'to', postName, 'Connectivity within', slctROI, sep=' '))
  
  return(conmatPlot)
}

structureMatrixPlotByType = function(conmatPlot){
  conmatPlot = conmatPlot +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(), 
          strip.placement = "outside", strip.background = element_rect(fill=NA, colour="grey50")) +
    facet_grid(reorder(type.from, desc(type.from)) ~ type.to, space="free", scales="free",switch="both")
  return(conmatPlot)
}


### Graph
#Reroganize to make graph with types instead of bodyids
reorganizeGraphData = function(fromData, toData, fromTypeData, toTypeData, weightData,relWeightData, cutoff){
  graphData = data.frame(from = fromTypeData, to = toTypeData,
                         fromid = fromData, toid = toData, 
                         weight = weightData, relWeight = relWeightData)
  
  graphData = graphData %>% group_by(from) %>% mutate(relWeightSumFrom = sum(relWeight)) %>% 
    group_by(from, to) %>% 
    mutate(connectTypeCount = length(fromid), relWeightNormFrom = mean(relWeightSumFrom/connectTypeCount)) %>%
    ungroup() %>% filter(weight > cutoff)
  
  return(graphData)
}

getGraphNodes = function(graphData){
  graphData_noSelf = graphData %>% filter(as.character(from) != as.character(to))
  nodes  = union(unique(graphData_noSelf$from), unique(graphData_noSelf$to))
  return(nodes)
}

getNoSelfGraphData = function(graphData){
  graphData_noSelf = graphData %>% filter(as.character(from) != as.character(to))
  return(graphData_noSelf)
}

#Reroganize to make graph with types instead of bodyids
getSelfFBGraphData = function(graphData){
  graphData_toSelf = graphData %>% filter(as.character(from) == as.character(to))

  graphData_selfFB = full_join(data.frame("from" = graphData_toSelf$from,
                                          #"weight" = graphData_toSelf$weight,
                                          "relWeightSelf" = graphData_toSelf$relWeight),
                               data.frame("from" = getGraphNodes(graphData)))
  #graphData_selfFB$weight[is.na(graphData_selfFB$weight)] <- 0
  graphData_selfFB$relWeightSelf[is.na(graphData_selfFB$relWeightSelf)] <- 0

  
  return(graphData_selfFB)
}

# convenient graph plotting
constructConnectivityGraph = function(graphData, cutoff, vertexSize, selfFBscale, arrowSize, edgeNorm, colormap=NULL, useRandCol=FALSE){
  source("colorCodeLookup.R")
  
  connectGraph = graph_from_data_frame(graphData)
  
  nodeCols = colors()[30+seq(1, length(V(connectGraph)$name))]
  if(!useRandCol){
    for (i in seq(1, length(V(connectGraph)$name))) {
      ncol = colors()[colorValueLookup$col[colorValueLookup$type ==  getSimpleTypeNames(V(connectGraph)$name[i])]]
      if (length(ncol) > 0) {nodeCols[i] = ncol}
    }
    if (! is_null(colormap)){
      for (i in seq(1, length(V(connectGraph)$name))) {
        ncol = unique(colormap %>% filter(Type %in% V(connectGraph)$name[i]) %>% select(hex))
        if (length(ncol$hex) > 0) {nodeCols[i] = ncol$hex}
      }
    }
  }
  
  # The labels are currently node IDs. Setting them to NA will render no labels
  V(connectGraph)$label.color="black"
  V(connectGraph)$label.cex=0.8
  V(connectGraph)$label.dist=0

  V(connectGraph)$size = vertexSize
  V(connectGraph)$vertex.frame.color="gray"
  V(connectGraph)$color=nodeCols

  # Set edge width based on weight:

  E(connectGraph)$width = E(connectGraph)$relWeight/edgeNorm

  #change arrow size and edge color:
  E(connectGraph)$arrow.size = arrowSize
  
  connectGraph = delete_edges(connectGraph, E(connectGraph)[relWeight<cutoff])
  
  
  return(connectGraph)
}

customizeGraphEdges  = function(connectGraph){
  
  edge.start <- ends(connectGraph, es=E(connectGraph), names=F)[,1]
  edge.col = V(connectGraph)$color[edge.start]
  return(edge.col)
}