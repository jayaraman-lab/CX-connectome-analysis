library(graphlayouts)
library(neuprintrExtra)

## Palettes
customROIPalette <-  function(){
  roiH <- getRoiTree()
  roiOutLabels <- selectRoiSet(roiH,default_level = 0)
  rPal <- roisPalette(rois = roiOutLabels)
  rPal["CX"] <- "white"
  rPal
}

p36 <- paletteer_d("Polychrome::palette36")
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

## MESHES AND SYNAPSES

displayAnatomies <- function(neurons=NULL,synapses=NULL,ROIs,saveName,neuronPalette=NULL,synapsePalette=NULL,synapseCluster="customContributor",
                             roiRef=selectRoiSet(getRoiTree(),exceptions=list("LAL(R)"=4,"CRE(R)"=4)),
                             roiPal=customROIPalette(),...){
  
  nopen3d()
  par3d(windowRect = c(30, 30, 1530, 1530))
  for (r in ROIs){
    locMesh <- neuprint_ROI_mesh(r)
    superROI <- roiRef$level0[match(r,roiRef$level2)]
    plot3d(locMesh,color=roiPal[superROI],alpha=ifelse(superROI=="CX",0.05,0.1),xlab="",ylab="",zlab="",box=FALSE,axes=FALSE,add=T)
  }
  par3d(scale=c(1,1,1))
  if (!is.null(neurons)){
    allMeshes <- lapply(neurons,function(bodyids){
      lapply(bodyids,getNeuronMesh,cloud=TRUE)
    })
    
    for (ty in names(allMeshes)){
      for (i in seq(length(allMeshes[[ty]]))){
        shade3d(allMeshes[[ty]][[i]],col=neuronPalette[ty])
      }
    }
  }
  
  if (!is.null(synapses)){
    for (cl in unique(synapses[[synapseCluster]])){
      plot3d(filter(synapses,!!(sym(synapseCluster))==cl)[,c("x","y","z")],col=synapsePalette[cl],add=T,...)
    }
  }
  
  nview3d("ventral")
  rgl.snapshot(paste0(outputsFolder,saveName,"-front.png"))
  nview3d("right",extramat=rotationMatrix(-pi/2, 1, 0, 0))
  rgl.snapshot(paste0(outputsFolder,saveName,"-side.png"))
  rgl.close()
}

## GRAPHS and SUBGRAPHS

plotSubgraph <- function(contributors,influenceThreshold=0.005,conns=mainFFConns,targets=mainFFTargetsS,graph=CX_outTblG,...){
  mainT <- filter(conns, Path_weight>influenceThreshold & mainContributor %in% contributors & type.to %in% (graph %N>% as_tibble())$type) 
  sourceSG <- as_tbl_graph(induced_subgraph(graph,c(contributors,mainT$type.to))) %>% 
    activate(nodes) %>% mutate(databaseType=type) %>% 
    mutate(supertype = targets$supertype[match(type,targets$type)]) %>% 
    mutate(supertype = ifelse(type %in% contributors,"Source",as.character(supertype)))
  sourceSGPlot <- ggraph(sourceSG,...) + geom_edge_fan(aes(width=weightRelative,color=.N()$supertype[from])) +
    geom_node_point(aes(color=supertype),size=3) + geom_node_text(aes(label=type),repel = T,size=2) + 
    scale_color_manual(name="Supertype",values=customGraphPalette) + scale_edge_color_manual(values=customGraphPalette)+ theme_paper_map() + scale_edge_width(name="Relative weight",limits=c(0,1)) +guides(edge_color="none")
  sourceSGPlot
}

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

## GRAPHS: OUT to CX MOTIFS

getMotifsGraphDf <- function(outDf,CXDf,ty,wr="relativeWeight4"){
  CXOuts <- unique(outDf$type.from)
  outDf <- filter(outDf,type.from %in% ty & !is.na(!!sym(wr))) %>% rename(weightRelative=!!sym(wr)) %>% mutate(motif="Out of CX pathways")
  if(nrow(outDf)>0){
    CXDf <- filter(CXDf, (type.from %in% unique(c(outDf$type.from,outDf$type.to))) & 
                     (type.to %in% unique(c(outDf$type.from,outDf$type.to)))) 
    if(nrow(CXDf)==0){CXDf <- mutate(CXDf,motif=character())}else{
      CXDf <- CXDf%>%
        mutate(motif=case_when(
          paste0(type.from,type.to) %in% paste0(outDf$type.from,outDf$type.to) ~ "CX Parallel connection",
          paste0(type.from,type.to) %in% paste0(outDf$type.to,outDf$type.from) & 
            !paste0(type.from,type.to) %in% paste0(outDf$type.from,outDf$type.to) ~ "CX Canonical feedback",
          (sapply(1:length(type.to),function(i){
            any(outDf[outDf$type.to == type.from[i],]$type.from %in% outDf[outDf$type.to == type.to[i],]$type.from)}) &
             !paste0(type.from,type.to) %in% paste0(outDf$type.from,outDf$type.to) &
             !paste0(type.from,type.to) %in% paste0(outDf$type.from,outDf$type.to)) ~"Linked targets in CX"
        )
        )}
    
    fullDf <- rbind(select(outDf,any_of(names(CXDf))),select(CXDf,any_of(names(outDf)))) %>% mutate(targetCat=ifelse(type.to %in% CXOuts,"CX output links","Other CX links"),
                                                                                                    sourceCat=ifelse(type.from %in% CXOuts,"CX output links","Other CX links"),
                                                                                                    category.to=ifelse(type.to %in% CXOuts,"CX output","Other CX"),
                                                                                                    category.from=ifelse(type.from %in% CXOuts,"CX output","Other CX"),
                                                                                                    category.to=factor(category.to,levels=c("Other CX","CX output")),
                                                                                                    category.from=factor(category.from,levels=c("Other CX","CX output")),
                                                                                                    layer.from=type.from %in% ty,
                                                                                                    layer.to=type.to %in% ty,
                                                                                                    motif=factor(motif,levels=c("Out of CX pathways","CX Parallel connection","CX Canonical feedback","Linked targets in CX")),
                                                                                                    linkCat = targetCat
    )
    fullDf$linkCat[fullDf$motif=="Reciprocal"] = fullDf$sourceCat[fullDf$motif=="Reciprocal"]
    
    return(fullDf %>% mutate(linkCat = factor(linkCat,levels=c("Out of CX pathways","CX Parallel connection","CX Canonical feedback","Linked targets in CX")),
                             roi=factor(roi,levels=c("Outside","EB","FB","PB","NO(R)","NO(L)")),
                             roiCat=ifelse(roi=="Outside","Outside","CX")))}
}

plotMotifs <- function(graphDf){
  focusContrib <- filter(graphDf,roi=="Outside") %>% group_by(type.to) %>% summarize(weightRelative=sum(weightRelative)) %>% ungroup()
  motifGG <- makeGraph(graphDf) 
  edgePal <- c("Plum",paletteer_d("ggthemes::Traffic")[1:3])
  names(edgePal) <- c("Out of CX pathways","CX Parallel connection","CX Canonical feedback","Linked targets in CX")
  
  #if(length(E(motifGG))>1){
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
           theme_paper_map() + scale_edge_color_manual(name="Motif",values=edgePal,breaks=c("Out of CX pathways","CX Parallel connection","CX Canonical feedback","Linked targets in CX"),drop=FALSE)+
           scale_edge_width(name="Relative weight/pathway weight",range=c(0.2,3),limits=c(0,1),)+
           geom_node_text(aes(label=name),size=1.5,nudge_y = 0.06,hjust="inward") + coord_fixed())
  #}
}

