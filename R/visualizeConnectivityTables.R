########################################################################################################################
# Visualization tools for connectivity tables
########################################################################################################################


### Connectivity matrix

plotConnectivityMatrix = function(myConTab, byGroup = "name", connectionMeasure="weightRelative", cmax=NULL, postfixfrom = NULL, postfixto = NULL) {
  #' Generate plot of connectivity matrix using ggplot
  #' @param myConTab Neuprint connection table
  #' @param synapseCutOff Minimum number of synapses between two partners to be considered a connection
  #' @param byType choose IDs or types
  #' @param connectionMeasure choose which weigh
  #' @return conmatPlot is a ggplot object
  # cool: low = "papayawhip", mid = "darkseagreen2", high = "steelblue4"
  # warm: 
  
  myConTab = myConTab %>% ungroup() %>% mutate(plotWeight = myConTab[[connectionMeasure]])
  
  if (byGroup == "id" | byGroup == "colname" | byGroup == "rowname"){
    myConTab$nameid.from = paste(as.character(myConTab$name.from), postfixfrom, as.character(myConTab$from), sep = "_")
    myConTab$nameid.to = paste(as.character(myConTab$name.to), postfixto,as.character(myConTab$to), sep = "_")
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
  } else if(byGroup == "colname"){
    conmatPlot =  conmatPlot + geom_tile(aes(colname.to,nameid.from,fill=plotWeight))
  } else if(byGroup == "rowname"){
    conmatPlot =  conmatPlot + geom_tile(aes(nameid.to,rowname.from,fill=plotWeight))
  } else if(byGroup == "id_sort"){
    conmatPlot =  conmatPlot + geom_tile(aes(nameOrder.to,nameOrder.from,fill=plotWeight))
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
    xlab("postsynaptic neuron")+ylab("presynaptic neuron") + labs(fill=connectionMeasure) +
    labs(title=paste(preName,'to', postName, sep=' '))
  
  return(conmatPlot)
}

structureMatrixPlotByType = function(conmatPlot){
  conmatPlot = conmatPlot +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(), 
          strip.placement = "outside", strip.background = element_rect(fill=NA, colour="grey50")) +
    facet_grid(reorder(type.from, desc(type.from)) ~ type.to, space="free", scales="free",switch="both")
  return(conmatPlot)
}

structureMatrixPlotByType_lines = function(conmatPlot, yonly=FALSE){
  if (yonly){
    conmatPlot = conmatPlot +
      facet_grid(reorder(type.from, desc(type.from)) ~ ., space="free", scales="free",switch="y") +
      theme(axis.text.y = element_blank(), 
            strip.placement = "outside", #strip.background = element_rect(fill=NA, colour="grey50"),
            strip.text.y.left = element_text(angle = 0),
            strip.text.x.bottom = element_text(angle = 90),
            strip.background = element_blank()) #remove background for facet labels
    return(conmatPlot)
  }else{
    conmatPlot = conmatPlot +
      facet_grid(reorder(type.from, desc(type.from)) ~ type.to, space="free", scales="free",switch="both") +
      theme(axis.text.x = element_blank(),axis.text.y = element_blank(), 
            strip.placement = "outside", #strip.background = element_rect(fill=NA, colour="grey50"),
            strip.text.y.left = element_text(angle = 0),
            strip.text.x.bottom = element_text(angle = 90),
            strip.background = element_blank(), #remove background for facet labels
            panel.border = element_rect(colour = "grey", fill = NA, size=0.3), #add black border
            panel.spacing = unit(0.1, "lines")) #remove space between facets
    return(conmatPlot)
  }
}


### Bar plot
inOutContributionPlot <- function(data){
  source("R/paperTheme.R")
  ncol = length(unique(data$partner))
  data = data %>% ungroup() %>% group_by(partner, roi, dir, ref) %>% summarise(measure = sum(measure))
  bar = ggplot(data, aes(x=ref, y=measure, fill=partner)) + 
    geom_bar(position = "fill", stat="identity") + 
    theme_classic() + theme_paper() + theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2))
  
  if(ncol <= 31){
    bar = bar + scale_fill_manual(values=paletteer_d("Polychrome::palette36", n=5+length(unique(data$partner)))[c(-seq(5))])
  }else if(ncol <= 36){
    bar = bar + scale_fill_manual(values=paletteer_d("Polychrome::palette36", n=length(unique(data$partner))))
  }#else{ ...add a better colormap than the default}
  return(bar)
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
constructConnectivityGraph = function(graphData, cutoff, vertexSize, selfFBscale, arrowSize, edgeNorm, colormap=NULL){
  connectGraph = graph_from_data_frame(graphData)
  
  nodeCols = colors()[30+seq(1, length(V(connectGraph)$name))]
  if (! is_null(colormap)){
    for (i in seq(1, length(V(connectGraph)$name))) {
      ncol = unique(colormap %>% filter(Type %in% V(connectGraph)$name[i]) %>% select(hex))
      if (length(ncol$hex) > 0) {
        #print(as.character(ncol$hex))
        nodeCols[i] = as.character(ncol$hex)
      }
    }
  }

  # The labels are currently node IDs. Setting them to NA will render no labels
  V(connectGraph)$label.color="black"
  V(connectGraph)$label.cex=0.8
  V(connectGraph)$label.dist=0

  V(connectGraph)$size = vertexSize
  V(connectGraph)$vertex.label.family="Arial"
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