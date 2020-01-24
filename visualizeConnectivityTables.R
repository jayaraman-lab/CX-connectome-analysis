########################################################################################################################
# Visualization tools for connectivity tables
########################################################################################################################


### Connectivity matrix

plotConnectivityMatrix = function(myConTab, synapseCutOff = 3, byType) {
  #' Generate plot of connectivity matrix using ggplot
  #' @param myConTab Neuprint connection table
  #' @param synapseCutOff Minimum number of synapses between two partners to be considered a connection
  #' @param byType choose IDs or types
  #' @return conmatPlot is a ggplot object

  conmatPlot = ggplot(myConTab  %>% filter(weight > synapseCutOff)) + 
    theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low = "grey", mid = "tan", high = "maroon4", limits=c(0,max(myConTab$weight)))
  if (byType){
    conmatPlot =  conmatPlot + geom_tile(aes(partnerName,name,fill=weight))
  }
  else{
    conmatPlot =  conmatPlot + geom_tile(aes(partnerid,nameid,fill=weight))
  }
  
  return(conmatPlot)
}

addMatrixPlotLabs = function(conmatPlot, preName, postName, slctROI) {
  #' Add labels to a connection matrix plot
  #' @param preName String describing presynaptic neurons
  #' @param postName String describing postsynaptic neurons
  #' @param slctROI Bane of ROI used to construct connection table
  #' @return conmatPlot is a ggplot object
  conmatPlot = conmatPlot +
    xlab("Post-synaptic neuron")+ylab("Pre-synaptic neuron") + labs(fill="# synapses") +
    labs(title=paste(preName,'to', postName, 'Connectivity within', slctROI, sep=' '))
  
  return(conmatPlot)
}

structureMatrixPlotByType = function(conmatPlot){
  conmatPlot = conmatPlot +
    facet_grid(type ~ partnerType, space="free", scales="free", switch="both") +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),strip.placement = "outside", 
          strip.background = element_rect(fill=NA, colour="grey50"))
  return(conmatPlot)
}


### Graph
