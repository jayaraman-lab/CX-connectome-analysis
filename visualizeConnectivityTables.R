########################################################################################################################
# Visualization tools for connectivity tables
########################################################################################################################


### Connectivity matrix

plotConnectivityMatrix = function(myConTab, byType, connectionMeasure="weightRelative") {
  #' Generate plot of connectivity matrix using ggplot
  #' @param myConTab Neuprint connection table
  #' @param synapseCutOff Minimum number of synapses between two partners to be considered a connection
  #' @param byType choose IDs or types
  #' @param connectionMeasure choose which weigh
  #' @return conmatPlot is a ggplot object
  # cool: low = "papayawhip", mid = "darkseagreen2", high = "steelblue4"
  # warm: 
  
  myConTab = myConTab %>% mutate(plotWeight = myConTab[[connectionMeasure]])
  
  if (!byType){
    myConTab$nameid = paste(as.character(myConTab$name), as.character(myConTab$bodyid), sep = "_")
    myConTab$partnerid = paste(as.character(myConTab$partnerName), as.character(myConTab$partner), sep = "_")
  }
  
  conmatPlot = ggplot(myConTab) + 
    theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="ivory", mid="peachpuff", high="black", limits=c(0,max(myConTab$plotWeight)))
  
  if (byType){
    conmatPlot =  conmatPlot + geom_tile(aes(partnerName,name,fill=plotWeight))
  }
  else{
    conmatPlot =  conmatPlot + geom_tile(aes(partnerid,nameid,fill=plotWeight))
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
    facet_grid(reorder(type, desc(type)) ~ partnerType, space="free", scales="free",switch="both")
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
                                          "weight" = graphData_toSelf$weight, 
                                          "relWeight" = graphData_toSelf$relWeightNormFrom),
                               data.frame("from" = nodes))
  graphData_selfFB$weight[is.na(graphData_selfFB$weight)] <- 0
  graphData_selfFB$relWeight[is.na(graphData_selfFB$relWeight)] <- 0
  
  return(graphData_selfFB)
}

# convenient graph plotting
constructConnectivityGraph = function(nodes, graphData_noSelf, graphData_selfFB, 
                                      cutoff, vertexSize, selfFBscale, arrowSize, edgeNorm, nodeCols){
  connectGraph = graph_from_data_frame(graphData_noSelf)
  connectGraph <- delete_edges(connectGraph, E(connectGraph)[relWeightNormFrom<cutoff])
  
  # The labels are currently node IDs. Setting them to NA will render no labels
  V(connectGraph)$label.color="black"
  V(connectGraph)$label.cex=0.8
  V(connectGraph)$label.dist=0

  V(connectGraph)$size = vertexSize + as.numeric(vertexSize*selfFBscale*graphData_selfFB$relWeight )
  V(connectGraph)$vertex.frame.color="gray"
  V(connectGraph)$color=nodeCols

  # Set edge width based on weight:
  E(connectGraph)$width <- E(connectGraph)$relWeightNormFrom/edgeNorm
  #change arrow size and edge color:
  E(connectGraph)$arrow.size <- arrowSize
  
  return(connectGraph)
}

customizeGraphEdges  = function(connectGraph){
  
  edge.start <- ends(connectGraph, es=E(connectGraph), names=F)[,1]
  edge.col = V(connectGraph)$color[edge.start]
  return(edge.col)
}