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
  # cool: low = "papayawhip", mid = "darkseagreen2", high = "steelblue4"
  # warm: 

  conmatPlot = ggplot(myConTab  %>% filter(weight > synapseCutOff)) + 
    theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="ivory", mid="peachpuff", high="black", limits=c(0,max(myConTab$weight)))
  if (byType){
    conmatPlot =  conmatPlot + geom_tile(aes(name,partnerName,fill=weight))
  }
  else{
    conmatPlot =  conmatPlot + geom_tile(aes(nameid,partnerid,fill=weight))
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
    facet_grid(partnerType ~ type, space="free", scales="free", switch="both") +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),strip.placement = "outside", 
          strip.background = element_rect(fill=NA, colour="grey50"))
  return(conmatPlot)
}

structureMatrixPlotByPreType = function(conmatPlot){
  conmatPlot = conmatPlot +
    facet_grid(partnerid ~ type, space="free", scales="free", switch="both") +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),strip.placement = "outside", 
          strip.background = element_rect(fill=NA, colour="grey50"))
  return(conmatPlot)
}

structureMatrixPlotByPostType = function(conmatPlot){
  conmatPlot = conmatPlot +
    facet_grid(partnerType ~ nameid, space="free", scales="free", switch="both") +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),strip.placement = "outside", 
          strip.background = element_rect(fill=NA, colour="grey50"))
  return(conmatPlot)
}


### Graph
#Reroganize to make graph with types instead of bodyids
reorganizeGraphData = function(fromData, toData, weightData, cutoff){
  graphData = data.frame(from = fromData, to = toData, weight = weightData)
  graphData = graphData %>% group_by(from, to) %>% 
    summarise(weight = mean(weight, na.rm = TRUE)) %>% ungroup() %>% 
    filter(weight > cutoff)
  
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
  graphData_selfFB = full_join(data.frame("from" = graphData_toSelf$from, "weight" = graphData_toSelf$weight), data.frame("from" = nodes))
  graphData_selfFB$weight[is.na(graphData_selfFB$weight)] <- 0
  
  return(graphData_selfFB)
}

# convenient graph plotting
constructConnectivityGraph = function(nodes, graphData_noSelf, graphData_selfFB, cutoff,
                                      arrowSize, vertexSize, edgeNorm){
  connectGraph = graph_from_data_frame(graphData_noSelf)
  connectGraph <- delete_edges(connectGraph, E(connectGraph)[weight<cutoff])
  
  simpleTypesNodes = getSimpleTypeNames(nodes)
  
  nodeCols = seq(1, length(nodes))
  for (i in seq(1, length(simpleTypesNodes))) {
    nodeCols[i] = colors()[colorValueLookup$col[colorValueLookup$type == simpleTypesNodes[i]]]
  }
  
  # The labels are currently node IDs. Setting them to NA will render no labels
  V(connectGraph)$label.color="black"
  V(connectGraph)$label.cex=0.8
  V(connectGraph)$label.dist=0
  
  typeCounts = typeCounts %>% filter(type %in% V(connectGraph)$name)
  typeCounts = typeCounts[match(V(connectGraph)$name, typeCounts$type),]
  V(connectGraph)$size = vertexSize + as.numeric(vertexSize*graphData_selfFB$weight/max(c(1, max(graphData_selfFB$weight) ) ) )
  V(connectGraph)$vertex.frame.color="gray"
  V(connectGraph)$color=nodeCols
  
  # Set edge width based on weight:
  E(connectGraph)$width <- E(connectGraph)$weight/edgeNorm
  #change arrow size and edge color:
  E(connectGraph)$arrow.size <- arrowSize
  edge.start <- ends(connectGraph, es=E(connectGraph), names=F)[,1]
  edge.col = V(connectGraph)$color[edge.start]
  
  return(connectGraph)
}


# color code
getSimpleTypeNames <- function(mydata){
  simpleTypes = gsub("Delta0[[:alnum:]]*", "Delta0", mydata)
  simpleTypes = gsub("Delta12[[:alnum:]]*", "Delta12", simpleTypes)
  simpleTypes = gsub("Delta6[[:alnum:]]*", "Delta6", simpleTypes)
  simpleTypes = gsub("Delta0-Delta12[(]*[[:alnum:]]*[)]*", "Delta0-12", simpleTypes)
  simpleTypes = gsub("FB[[:alnum:]]*", "FB", simpleTypes)
  simpleTypes = gsub("FB_[[:alnum:]]*", "FB", simpleTypes)
  simpleTypes = gsub("FB-Q[[:alnum:]]*", "FB-Q", simpleTypes)
  simpleTypes = gsub("EB-Q[[:alnum:]]*", "EB-Q", simpleTypes)
  #simpleTypes = gsub("ExR[[:alnum:]]*", "ExR", simpleTypes)
  simpleTypes = gsub("PB[[:alnum:]]*", "PB", simpleTypes)
  simpleTypes = gsub("SLP-AB[[:alnum:]]*", "SLP-AB", simpleTypes)
  simpleTypes = gsub("R1[[:alnum:]]*_[[:alnum:]]*", "R1", simpleTypes)
  simpleTypes = gsub("R2[[:alnum:]]*_[[:alnum:]]*", "R2", simpleTypes)
  simpleTypes = gsub("R3a_[[:alnum:]]*", "R3a", simpleTypes)
  simpleTypes = gsub("R3d_[[:alnum:]]*", "R3d", simpleTypes)
  simpleTypes = gsub("R3m_[[:alnum:]]*", "R3m", simpleTypes)
  simpleTypes = gsub("R3p_[[:alnum:]]*", "R3p", simpleTypes)
  simpleTypes = gsub("R3w_[[:alnum:]]*", "R3w", simpleTypes)
  simpleTypes = gsub("R4d_[[:alnum:]]*", "R4d", simpleTypes)
  simpleTypes = gsub("R4m[[:alnum:]]*", "R4m", simpleTypes)
  simpleTypes = gsub("R5[[:alnum:]]*", "R5", simpleTypes)
  simpleTypes = gsub("R6[[:alnum:]]*", "R6", simpleTypes)
  #simpleTypes = gsub("TuBu[[:alnum:]]*", "TuBu", simpleTypes)
  simpleTypes = gsub("TuBu09_[[:alnum:]]*", "TuBu09", simpleTypes)
  simpleTypes = gsub("PDM21a_[[:alnum:]]*", "PDM21a", simpleTypes)
  simpleTypes = gsub("PEN_a\\(PEN1\\)", "PEN1", simpleTypes)
  simpleTypes = gsub("PEN_b\\(PEN2\\)", "PEN2", simpleTypes)
  simpleTypes = gsub("MC[[:alnum:]]*", "MC", simpleTypes)
  simpleTypes = gsub("TuTu[[:alnum:]]*", "TuTu", simpleTypes)
  simpleTypes = gsub("ADM06b_[[:alnum:]]*_[[:alnum:]]*", "ADM06b", simpleTypes)
  simpleTypes = gsub("MBON[[:alnum:]]*", "MBON", simpleTypes)
  simpleTypes = gsub("AVL[[:alnum:]]*_[[:alnum:]]*_pct", "AVL", simpleTypes)
  simpleTypes = gsub("AVL[[:alnum:]]*_pct", "AVL", simpleTypes)
  simpleTypes = gsub("AVM[[:alnum:]]*_pct", "AVM", simpleTypes)
  simpleTypes = gsub("ADL.*_pct", "ADL", simpleTypes)
  simpleTypes = gsub("PVL[[:alnum:]]*_pct", "PVL", simpleTypes)
  simpleTypes = gsub("PVM.*_pct", "PVM", simpleTypes)
  simpleTypes = gsub("PDL27e_pct", "PDL27e", simpleTypes)
  simpleTypes = gsub("PDL[12x,20s]{1}.*_pct", "PDLother", simpleTypes)
  simpleTypes = gsub("ADM11[[:alnum:]]*_pct", "ADM11", simpleTypes)
  simpleTypes = gsub("ADM06p_pct", "ADM06p", simpleTypes)
  simpleTypes = gsub("ADM06d_pct", "ADM06d", simpleTypes)
  simpleTypes = gsub("ADM03.*_pct", "ADM03", simpleTypes)
  simpleTypes = gsub("olfactory multi .*", "other", simpleTypes)
  simpleTypes = gsub("PDM09[[:alnum:]]*_pct", "PDM09", simpleTypes)
  simpleTypes = gsub("PDM26[[:alnum:]]*_pct", "PDMother", simpleTypes)
  simpleTypes = gsub("PDM14j.*_pct", "PDM14j", simpleTypes)
  simpleTypes = gsub("PDM14[m,r,d]{1}.*_pct", "PDM14other", simpleTypes)
  simpleTypes = gsub("PDM28[[:alnum:]]*_pct", "PDMother", simpleTypes)
  return(simpleTypes)
}

colorValueLookup = data.frame(
  type = c('R1', 'R2', 'R3a', 'R3d', 'R3m', 'R3p',  'R3w', 'R4d', 'R4m', 'R5', 'R6',
           'ExR1','ExR2','ExR3','ExR4','ExR5','ExR6','ExR7','ExR8',
           'TuBu01', 'TuBu02', 'TuBu03', 'TuBu04', 'TuBu05', 'TuBu06', 'TuBu07', 'TuBu08', 'TuBu09', 'TuBu10',
           'PDM21a', 'MC',  'TuTu',
           'EPG', 'EPGt', 'PEN1', 'PEN2', 'PEG', 'EQ5',
           'PDM14j', 'ADM06d','ADM06p', 'ADM06b','PDL27e','ADM06s', 'PDM09','PDMother','PDM14other','ADM03',
           'AVL', 'AVM','ADL','ADM11', 'MBON', 'PVL', 'PDM', 'PDLother','PVM',
           'FB', 'other' ),
  col = c( 367 ,   9,   34,    101,   32,    21,     58,    11,    12,    657,  517,
           468,   456,  467,   463,   464,   465,   466,    98,
           592,   591,   590,   589,   616,   617,   618,   619,   128,   130,
           259,   600,  26,
           499,    499,    143,   144,    573,  640, 
           639, 430, 431, 630, 452,632, 103, 104,105, 651,
           420, 420, 420, 420, 76, 535, 535,535,535,
           563, 651)
)