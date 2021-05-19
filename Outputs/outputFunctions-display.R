library(graphlayouts)
library(tidygraph)
library(ggraph)
library(igraph)
library(neuprintrExtra)
library(neuprintr)
library(nat)
library(paletteer)
library(dplyr)
library(stringr)
library(patchwork)

## Palettes################################################################

#' A shortcut for the polychrome 36 palette
p36 <- paletteer_d("Polychrome::palette36")

#' A custom palette for supertypes level 3 to be used throughout the output figures
supertype3Palette <- c("Unassigned"="grey80",
                       "Terra incognita"=p36[2],
                       "EB Columnar"=supertype2Palette()$pal["EPG"],
                       "ER"=supertype2Palette()$pal["ER"],
                       "ExR"=supertype2Palette()$pal["ExR"],
                       "FB Tangential"=supertype2Palette()$pal["FBt"],
                       "FB Output"=supertype2Palette()$pal["FS"],
                       "FB Columnar"=supertype2Palette()$pal["PFL"],
                       "PB Input"=supertype2Palette()$pal["SPS-PB"],
                       "LNO"=supertype2Palette()$pal["LNO"],
                       "SA"=supertype2Palette()$pal["SA"],
                       "5HT"=p36[14],
                       "OA"=p36[16],
                       "Peptidergic"=p36[28],
                       "DAN"=p36[29],
                       "MBON"=p36[4],
                       "LH"=p36[31],
                       "DN"=p36[24],
                       "Fru"=p36[30],
                       "KC"=p36[13],
                       "Antennal lobe"=p36[20],
                       "Clock"=p36[6],
                       "Visual PNs"=p36[32],
                       "Other Sensory"=p36[3])

#' A palette for the supertypes customized (a mix of different granularities) used in the output section
customSupertypePalette <- c(
  "EB Columnar"=supertype2Palette()$pal["EPG"],
  "ER6"=supertype2Palette()$pal["ER"],
  "PFL3"=supertype2Palette()$pal[["PFL"]],
  "PFL2"="#cb8801",
  "PFL1"="#ffe5b3",
  "ExR2-6"="#66a3ff",
  "ExR8"="#0148b2",
  "ExR7"=supertype2Palette()$pal[["ExR"]],
  "ExR1"="#b3d1ff",
  "PFR_b"=supertype2Palette()$pal[["PFR"]],
  "PFR_a"="#b6afb3",
  "FC1"=supertype2Palette()$pal[["FC"]],
  "FC2"="#009978",
  "FS2"=supertype2Palette()$pal[["FS"]],
  "FS4"="#8c4073",
  "FS3"="#af508e",
  "FS1"="#e7cbdd",
  "FB5"=supertype2Palette()$pal[["FBt"]],
  "FB4"="#6e8415",
  "FB6-7"="#bedf3a",
  "FB8-9"="#e2f1a7",
  "FR2"=supertype2Palette()$pal[["FR"]],
  "FR1"="#c40812",
  "MBON"=supertype3Palette["MBON"]
)

## Custom connectivity plots#######################################################################################

#' Plot per columns/glomeruli connectivity matrices for columnar outputs
#' @param bag A \code{neuronBag} object containing the output connection of a columnar neuron type of interest
#' @param type The columnar type (or set of columnar types) one wants to analyze the outputs of
#' @param targetFilt A subset of the neuron's target one wants to restrict the figure to
#' @param grouping Whether to group per "glomerulus" or "column"
#' @param facetInputs A variable to facet the inputs on (usually not used)
#' @details The grouping is determined using the "_CX" or "_LX,RX" parts of the neuron names. 
#' This uses \code{plotConnectivity} internally.
#' 
plotGlomMat <- function(bag,
                        type,
                        targetFilt=NULL,
                        grouping=c("glomerulus","column"),
                        facetInputs=NULL){
  grouping <- match.arg(grouping)
  outRaw <- bag$outputs_raw
  if (!is.null(targetFilt)) outRaw <- outRaw %>% filter(type.to %in% targetFilt$type.to)
  outRaw <- outRaw %>% 
    filter(!type.to %in%  type.from & type.from %in% type) %>%
    arrange(type.to) %>% 
    mutate(glomerulus=factor(gsub("_","",str_extract(name.from,"_[L|R][1-9]_")),
                             levels=c(paste0("L",9:1),paste0("R",1:9))),
           column=factor(gsub("_","",str_extract(name.from,"_C[1-9]")),
                         levels=paste0("C",9:1))) %>% 
    group_by(glomerulus,type.from) %>% mutate(n_glom_from=length(unique(from))) %>%
    group_by(column,type.from) %>% mutate(n_col_from=length(unique(from))) %>% ungroup()
  
  
  groupingName <- ifelse(grouping=="glomerulus","glom","col")
  outPerGroup <- group_by(outRaw,!!(sym(grouping)),
                          type.from,type.to,supertype3.to,supertype2.to,to,name.to,roi) %>% 
    summarize(weightRelative=sum(weightRelative),
              name.from=paste0((!!(sym(grouping)))[1],
                               " (n=",(!!sym(paste0("n_",groupingName,"_from")))[1],")"),
              from=(!!(sym(grouping)))[1]) %>% ungroup() %>% 
    arrange(type.to,!!(sym(grouping)))
  
  if (grouping=="glomerulus"){
    orderIn <-  c(paste0("L",9:1),paste0("R",1:9))
  }else{
    orderIn <- paste0("C",9:1)
  }
  
  connMat <- plotConnectivity(outPerGroup,slctROI=outPerGroup$roi[1],
                              grouping="neuron",xaxis="outputs",
                              replacementLabels = "name",orderIn = orderIn,
                              legendName="relative weight",
                              theme=theme_paper_rects(),
                              facetOutputs="supertype2.to",facetInputs=facetInputs) + 
    xlab("post synaptic neuron") + ylab(paste(type,grouping)) +theme(strip.text.x=element_blank())
  
  
  groupCounts <- outPerGroup %>% 
    group_by(name.from,type.from,!!sym(grouping)) %>% 
    summarize(totalW=sum(weightRelative)) %>% ungroup() %>% 
    mutate(!!grouping := factor(!!(sym(grouping)),levels=orderIn)) %>% arrange(!!(sym(grouping))) %>%
    mutate(name.from = factor(name.from,levels=unique(name.from)))
  
  countPlot <- ggplot(groupCounts,aes(x=name.from,y=totalW,group=type.from)) + 
    geom_line() +theme_paper_hgrid() + coord_flip() + 
    theme_paper(axis.text.x=element_text(angle=90),
                axis.text.y=element_blank(),axis.line.y=element_blank(),
                axis.title.y=element_blank(),axis.ticks.y = element_blank()) + 
    ylim(c(0,NA)) + ylab("total relative weight")
  
  connMat + countPlot + plot_layout(guides="collect",widths=c(1,0.1))
}

## 3D synapses plots ##########################################################

#'3D Synapse plots with ROI renderings
#'@param synapses A data frame of synapses (as returned by \code{neuprint_get_synapses}, eventually enriched with metadata columns)
#'@param ROIs Regions of interest to display, as a vector of strings
#'@param saveFolder The name of the folder where to save the resulting image. If NULL, nothing is saved. Otherwise, the rgl window is closed after saving (to avoid accumulating windows in a script)
#'@param saveName The name of the file to save the rendering as (it will be saved as a .png)
#'@param synapsePalette A palette to apply to the synapses
#'@param synapseCluster A variable of the \code{synapses} data frame to map the palette onto
#'@param roiRef A ROI data frame as returned by \code{selectRoiSet} to be used as a reference for the ROI palette
#'@param roiPal A palette to color the ROIs with
#'@param alphaRois Alpha value for the ROI renderings
#'@param size The window size. Increasing the size is the easiest way to obtain higher resolutions images
#'@param ... To be passed to \code{nat::plot3d}
displaySynapses3D <- function(synapses=NULL,
                             ROIs,
                             saveFolder=NULL,
                             saveName="rendering",
                             synapsePalette=NULL,
                             synapseCluster="customContributor",
                             roiRef=selectRoiSet(getRoiTree(),
                                                 exceptions=list("LAL(R)"=4,"CRE(R)"=4)),
                             roiPal=customROIPalette(),
                             alphaRois=0.05,
                             windowSize=c(1500,1500),...){
  
  nopen3d()
  par3d(windowRect = c(30, 30, windowSize[1]+30, windowSize[2]+30))
  for (r in ROIs){
    locMesh <- neuprint_ROI_mesh(r)
    superROI <- roiRef$level0[match(r,roiRef$level2)]
    plot3d(locMesh,color=roiPal[superROI],
           alpha=alphaRois,
           specular="black",
           xlab="",ylab="",zlab="",box=FALSE,axes=FALSE,add=T)
  }
  par3d(scale=c(1,1,1))
  
  if (!is.null(synapses)){
    for (cl in unique(synapses[[synapseCluster]])){
      plot3d(filter(synapses,!!(sym(synapseCluster))==cl)[,c("x","y","z")],
             col=synapsePalette[cl],add=T,...)
    }
  }
  
  nview3d("ventral")
  if (!is.null(saveFolder)){
    rgl.snapshot(file.path(saveFolder,paste0(saveName,"-front.png")))
    nview3d("right",extramat=rotationMatrix(-pi/2, 1, 0, 0))
    rgl.snapshot(file.path(saveFolder,paste0(saveName,"-side.png")))
    rgl.close()
  }
}

#' The default palette used in those plots
customROIPalette <-  function(){
  roiH <- getRoiTree()
  roiOutLabels <- selectRoiSet(roiH,default_level = 0)
  rPal <- roisPalette(rois = roiOutLabels)
  rPal["CX"] <- "grey60"
  rPal
}

## GRAPHS and SUBGRAPHS

#' A simple network graph plot using \code{ggraph}
#' @param gr A \code{tidygraph} object
#' @param pal A palette to use for supertypes
#' @param colP The variable to map onto the supertype color code
#' @param loop Whether to show loop connections
#' @param statW The statistic to use for the edges width
#' @param ... To be passed to \code{ggraph} (useful to pass layouts for example)
standardGraph <- function(gr,pal,colP="customSupertype",loop=F,statW="weightRelative",widthRange=c(0.2,3),widthLimits=c(0.001,NA),...){
  sG <- ggraph(gr,...)+
    geom_edge_fan(aes(color=!!sym(paste0(colP,".from")),
                      width=!!sym(statW)),
                  end_cap = circle(4, "mm"),
                  linejoin = "mitre",linemitre=3,
                  arrow =arrow(length = unit(1, 'mm'),type = "closed")) +
    theme_paper_map() 
  if (loop){
    sG <- sG + geom_edge_loop(aes(color=!!sym(paste0(colP,".from")),
                                  width=weightRelative),end_cap = circle(3, "mm"),
                              linejoin = "mitre",linemitre=3,
                              arrow =arrow(length = unit(1, 'mm'),type = "closed"))
  }
  sG + 
    geom_node_point(aes(fill=!!sym(colP)),size=4,shape=21)+
    geom_node_text(aes(label=type))+
    scale_fill_manual(values=pal,name="supertype")+ 
    scale_edge_color_manual(values=pal) + 
    scale_edge_width(range=widthRange,limits=widthLimits,name="relative weight") + 
    guides(edge_color="none")
  
}


#' Get the graph reduced to the strong targets of a selection of output types
#' \code{plotSubgraph_subgraph} just returns the \code{tidygraph} object, 
#' \code{plotSubgraph} plots it using tidygraph
#' @param contributors The output neurons to consider (as a vector of strings)
#' @param influenceThreshold Only keep targets receiving at least \code{influenceThresholds} path weights from 
#' the output neurons of interest
#' @param colP which variable of the "targets" table should map to colors
#' @param pal A palette to use, mapping to the colP variable of the "targets" table, 
#' with the "contributors" being labeled as "Source"
#' @param conns A data frame of the main connections of the output neurons (\code{contributors} is a subset of the source types in that table)
#' @param targets A metadata data frame of the main targets in the same set as \code{conns}
#' @param graph A graph object to be reduced
#' @param ... To be passed to \code{standardGraph}
plotSubgraph_subgraph <- function(contributors,
                                  influenceThreshold=0.005,
                                  conns,
                                  targets,
                                  graph,
                                  colP="supertype"){
  mainT <- filter(conns, Path_weight>influenceThreshold & 
                    type.from %in% contributors & 
                    type.to %in% (graph %N>% as_tibble())$type) 
  
  sourceSG <- as_tbl_graph(induced_subgraph(graph,c(contributors,mainT$type.to))) %>% 
    activate(nodes) %>% 
    mutate(grouped = targets[[colP]][match(type,targets$type)]) %>% 
    mutate(grouped = ifelse(type %in% contributors,"Source",as.character(grouped))) %E>%
    mutate(grouped.from=.N()$grouped[match(type.from,.N()$type)])
  
  sourceSG
}

plotSubgraph <- function(contributors,
                         influenceThreshold=0.005,
                         pal=customSupertypePalette,
                         conns,
                         targets,
                         graph,
                         colP="supertype",
                         ...){
  sourceSG <- plotSubgraph_subgraph(contributors,influenceThreshold,conns,targets,graph,colP=colP)
  sourceSGPlot <- standardGraph(sourceSG,pal=pal,colP="grouped",...)
  sourceSGPlot
}

#' Return the intersection of the downstream neighborhood of one type with the 
#' upstream neighborhood of another type
#' @param type.from A neuron to explore the downstream neighborhood of
#' @param type.to A neuron to explore the upstream neighborhood of
#' @param graph A \code{tidygraph} object
#' @param order The neighborhood order (1 for direct neighbors, 2, for neighbors up to 2 steps away, etc)
#' @details Internally uses the \code{ego} function from \code{igraph}
getLocalSubgraph <- function(type.from,type.to,graph,order=2){
  to_searchOut <- unlist(sapply(ego(graph,nodes=type.from,mode="out",order=order),names))
  to_searchIn <- unlist(sapply(ego(graph,nodes=type.to,mode="in",order=order),names))
  to_search <- intersect(to_searchIn,to_searchOut)
  to_search <- c(type.from,type.to,to_search)
  graph <- igraph::induced_subgraph(graph,to_search)
  paths <- unlist(lapply(type.from,function(ty) all_simple_paths(graph,from=ty,to=type.to,mode="out")))
  nodes <- unique(names(paths))
  graph <- as_tbl_graph(induced_subgraph(graph,nodes))
}

## Motifs ###########################################################################

#' Extract and label motifs from the connectivity tables of two different regions
#' @param outDf Reference connectivity table (usually the "out of the CX table")
#' @param CXDf Connectivity table to look for motifs into
#' @param ty Reference types to consider 
#' @param wr Statistic to use in the reference table ("Path_weight" by default)
#' @return A table restricted to connections participating in motifs, labeled by motifs (column \code{linkCat})
getMotifsGraphDf <- function(outDf,
                             CXDf,
                             ty,
                             wr="Path_weight"){
  CXOuts <- unique(outDf$type.from)
  
  outDf <- filter(outDf,type.from %in% ty & !is.na(!!sym(wr))) %>% 
    rename(weightRelative=!!sym(wr)) %>% mutate(motif="out of CX pathways")
  
  if(nrow(outDf)>0){
    CXDf <- filter(CXDf, (type.from %in% unique(c(outDf$type.from,outDf$type.to))) & 
                     (type.to %in% unique(c(outDf$type.from,outDf$type.to)))) 
    
    if(nrow(CXDf)==0){CXDf <- mutate(CXDf,motif=character())}else{
      
      CXDf <- CXDf%>%
        mutate(motif=case_when(
          paste0(type.from,type.to) %in% paste0(outDf$type.from,outDf$type.to) ~ "CX parallel connection",
          paste0(type.from,type.to) %in% paste0(outDf$type.to,outDf$type.from) & 
            !paste0(type.from,type.to) %in% paste0(outDf$type.from,outDf$type.to) ~ "CX canonical feedback",
          (sapply(1:length(type.to),function(i){
            any(outDf[outDf$type.to == type.from[i],]$type.from %in% outDf[outDf$type.to == type.to[i],]$type.from)}) &
             !paste0(type.from,type.to) %in% paste0(outDf$type.from,outDf$type.to) &
             !paste0(type.from,type.to) %in% paste0(outDf$type.from,outDf$type.to)) ~"linked targets in CX"
        )
        )}
    
    fullDf <- rbind(select(outDf,any_of(names(CXDf))),select(CXDf,any_of(names(outDf)))) %>% 
      mutate(targetCat=ifelse(type.to %in% CXOuts,"CX output links","Other CX links"),
             sourceCat=ifelse(type.from %in% CXOuts,"CX output links","Other CX links"),
             category.to=ifelse(type.to %in% CXOuts,"CX output","Other CX"),
             category.from=ifelse(type.from %in% CXOuts,"CX output","Other CX"),
             category.to=factor(category.to,levels=c("Other CX","CX output")),
             category.from=factor(category.from,levels=c("Other CX","CX output")),
             layer.from=type.from %in% ty,
             layer.to=type.to %in% ty,
             motif=factor(motif,levels=c("out of CX pathways","CX parallel connection","CX canonical feedback","linked targets in CX")),
             linkCat = targetCat
      )
    
    return(fullDf %>% mutate(linkCat = factor(linkCat,
                                              levels=c("out of CX pathways",
                                                       "CX parallel connection",
                                                       "CX canonical feedback",
                                                       "linked targets in CX")),
                             roi=factor(roi,levels=c("Outside","EB","FB","PB","NO(R)","NO(L)")),
                             roiCat=ifelse(roi=="Outside","Outside","CX")))}
}

#' Plot the results of \code{getMotifsGraphDf} in a circular layout
#' @param graphDf A table of motifs as returned by \code{getMotifsGraphDf}
plotMotifs <- function(graphDf){
  focusContrib <- filter(graphDf,roi=="Outside") %>% group_by(type.to) %>% summarize(weightRelative=sum(weightRelative)) %>% ungroup()
  motifGG <- makeGraph(graphDf) 
  edgePal <- c("Plum",paletteer_d("ggthemes::Traffic")[1:3])
  names(edgePal) <- c("out of CX pathways","CX parallel connection","CX canonical feedback","linked targets in CX")
  
  motifGG <- motifGG %N>% 
    mutate(
      customSuper = ifelse(layer,"Source",supertype2))
  
  customPal <- c("Source"="grey50",supertype2Palette()$pal)
  
  return(ggraph(motifGG,layout="focus",focus=(motifGG %N>% as_tibble())$layer)+
           graphlayouts::draw_circle(max.circle = 1,col="grey90") +
           geom_edge_fan2(aes(color=motif,
                              width=weightRelative),
                          end_cap = circle(3, "mm"),linejoin = "mitre",linemitre=3,
                          arrow =arrow(length = unit(1, 'mm'),type = "closed"))+
           geom_node_point(aes(color=customSuper),size=4)+
           scale_color_manual(values=customPal) +
           guides(color="none") +
           theme_paper_map() + scale_edge_color_manual(name="motif",values=edgePal,breaks=c("out of CX pathways","CX parallel connection","CX canonical feedback","linked targets in CX"),drop=FALSE)+
           scale_edge_width(name="relative weight/pathway weight",range=c(0.2,3),limits=c(0,1),)+
           geom_node_text(aes(label=name),nudge_y = 0.15) + 
           coord_fixed(clip="off"))
}

