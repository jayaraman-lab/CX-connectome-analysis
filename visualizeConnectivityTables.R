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
  return(simpleTypes)
}


colorValueLookup = data.frame(
  type = c('R1', 'R2', 'R3a', 'R3d', 'R3m', 'R3p',  'R3w', 'R4d', 'R4m', 'R5', 'R6',
           'ExR1','ExR2','ExR3','ExR4','ExR5','ExR6','ExR7','ExR8',
           'TuBu01', 'TuBu02', 'TuBu03', 'TuBu04', 'TuBu05', 'TuBu06', 'TuBu07', 'TuBu08', 'TuBu09', 'TuBu10',
           'PDM21a', 'MC',
           'EPG', 'EPGt', 'PEN1', 'PEN2', 'PEG', 'EQ5'),
  col = c( 367 ,   9,   34,    101,   32,    21,     58,    11,    12,    657,  517,
           468,   456,  467,   463,   464,   465,   466,    98,
           592,   591,   590,   589,   616,   617,   618,   619,   128,   130,
           259,   600,
           499,    499,    143,   144,    573,    640)
)