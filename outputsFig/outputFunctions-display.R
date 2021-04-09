library(graphlayouts)
library(tidygraph)
library(ggraph)
library(igraph)
library(neuprintrExtra)
library(neuprintr)
library(nat)
library(paletteer)

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

##
standardGraph <- function(gr,pal,colP="customSupertype",...){ggraph(gr,...) + 
  geom_edge_fan(aes(color=!!sym(paste0(colP,".from")),width=weightRelative),end_cap = circle(3, "mm"),linejoin = "mitre",linemitre=3,
                arrow =arrow(length = unit(1, 'mm'),type = "closed")) + geom_node_point(aes(color=!!sym(colP)),size=4)+
  geom_node_text(aes(label=type))+theme_paper_map()+scale_color_manual(values=pal,name="Supertype")+ 
  scale_edge_color_manual(values=pal) + 
    scale_edge_width(range=c(0.2,3),limits=c(0.001,NA),name="Relative weight") + 
  guides(edge_color="none")}

## Columns/glomeruli matrices + summaries for columnar neurons
plotGlomMat <- function(bag,type,targetFilt=mainFFTargets,grouping=c("glomerulus","column")){
  grouping <- match.arg(grouping)
  outRaw <- bag$outputs_raw %>% 
    filter(type.to %in% targetFilt$type.to & !type.to %in%  type.from) %>%
    filter(type.from %in% type) %>%
    arrange(type.to) %>% 
    mutate(glomerulus=factor(gsub("_","",str_extract(name.from,"_[L|R][1-9]_")),levels=c(paste0("L",9:1),paste0("R",1:9))),
           column=factor(gsub("_","",str_extract(name.from,"_C[1-9]")),levels=paste0("C",9:1))) %>% 
    group_by(glomerulus,type.from) %>% mutate(n_glom_from=length(unique(from))) %>%
    group_by(column,type.from) %>% mutate(n_col_from=length(unique(from))) %>% ungroup()
  
  
  groupingName <- ifelse(grouping=="glomerulus","glom","col")
  outPerGroup <- group_by(outRaw,!!(sym(grouping)),type.from,type.to,supertype3.to,supertype2.to,to,name.to,roi) %>% 
    summarize(weightRelative=sum(weightRelative),
              name.from=paste0((!!(sym(grouping)))[1]," (n=",(!!sym(paste0("n_",groupingName,"_from")))[1],")"),
              from=(!!(sym(grouping)))[1]) %>% ungroup() %>% 
    arrange(type.to,!!(sym(grouping)))
  
  if (grouping=="glomerulus"){
    orderIn <-  c(paste0("L",9:1),paste0("R",1:9))
  }else{
    orderIn <- paste0("C",9:1)
  }
  
  connMat <- plotConnectivity(outPerGroup,slctROI=outPerGroup$roi[1],grouping="neuron",xaxis="outputs",replacementLabels = "name",orderIn = orderIn,legendName="Relative weight",theme=theme_paper_grid(),facetOutputs="supertype2.to",facetInputs=facetInputs) + 
    xlab("post synaptic neuron") + ylab(paste(type,grouping)) +theme(strip.text.x=element_blank())
  
  
  groupCounts <- outPerGroup %>% group_by(name.from,type.from,!!sym(grouping)) %>% summarize(totalW=sum(weightRelative)) %>% ungroup() %>% 
    mutate(!!grouping := factor(!!(sym(grouping)),levels=orderIn)) %>% arrange(!!(sym(grouping))) %>%
    mutate(name.from = factor(name.from,levels=unique(name.from)))
  
  countPlot <- ggplot(groupCounts,aes(x=name.from,y=totalW,group=type.from)) + geom_line() +theme_paper_hgrid() + coord_flip() + 
    theme(axis.text.x=element_text(angle=90),axis.text.y=element_blank(),axis.line.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y = element_blank()) + ylim(c(0,NA)) + ylab("total relative weight")
    
  connMat + countPlot + plot_layout(guides="collect",widths=c(1,0.1))
}

## MESHES AND SYNAPSES

displayAnatomies <- function(neurons=NULL,
                             synapses=NULL,
                             ROIs,s
                             aveName="rendering",
                             saveFolder=".",
                             neuronPalette=NULL,synapsePalette=NULL,
                             synapseCluster="customContributor",
                             roiRef=selectRoiSet(getRoiTree(),exceptions=list("LAL(R)"=4,"CRE(R)"=4)),
                             alphaRois=0.2,
                             roiPal=customROIPalette(),
                             roiShininess=1,
                             size=c(1500,1500),...){
  
  nopen3d()
  par3d(windowRect = c(30, 30, size[1]+30, size[2]+30))
  for (r in ROIs){
    locMesh <- neuprint_ROI_mesh(r)
    superROI <- roiRef$level0[match(r,roiRef$level2)]
    plot3d(locMesh,color=roiPal[superROI],alpha=alphaRois,shininess=roiShininess,xlab="",ylab="",zlab="",box=FALSE,axes=FALSE,add=T)
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
  if (!is.null(saveName)){
    rgl.snapshot(file.path(saveFolder,paste0(saveName,"-front.png")))
    nview3d("right",extramat=rotationMatrix(-pi/2, 1, 0, 0))
    rgl.snapshot(file.path(saveFolder,paste0(saveName,"-side.png")))
    rgl.close()
  }
}

## GRAPHS and SUBGRAPHS

plotSubgraph <- function(contributors,
                         influenceThreshold=0.005,
                         pal=customGraphPalette,
                         colP="supertype",
                         conns=mainFFConns,
                         targets=mainFFTargetsS,
                         graph=CX_outTblG,...){
  mainT <- filter(conns, Path_weight>influenceThreshold & type.from %in% contributors & type.to %in% (graph %N>% as_tibble())$type) 
  sourceSG <- as_tbl_graph(induced_subgraph(graph,c(contributors,mainT$type.to))) %>% 
    activate(nodes) %>% mutate(databaseType=type) %>% 
    mutate(supertype = targets$supertype[match(type,targets$type)]) %>% 
    mutate(supertype = ifelse(type %in% contributors,"Source",as.character(supertype))) %E>%
    mutate(supertype.from=.N()$supertype[match(type.from,.N()$type)])
  sourceSGPlot <- standardGraph(sourceSG,pal=pal,colP=colP,...)#ggraph(sourceSG,...) + geom_edge_fan(aes(width=weightRelative,color=.N()$supertype[from])) +
   # geom_node_point(aes(color=supertype),size=3) + geom_node_text(aes(label=type),repel = T,size=2) + 
   #scale_color_manual(name="Supertype",values=c(customSupertypePalette,customGraphPalette,"Source"="grey30")) + scale_edge_color_manual(values=c(customSupertypePalette,customGraphPalette,"Source"="grey30"))+ theme_paper_map() + scale_edge_width(name="Relative weight",limits=c(0,1)) +guides(edge_color="none")
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
  outDf <- filter(outDf,type.from %in% ty & !is.na(!!sym(wr))) %>% rename(weightRelative=!!sym(wr)) %>% mutate(motif="out of CX pathways")
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
    
    fullDf <- rbind(select(outDf,any_of(names(CXDf))),select(CXDf,any_of(names(outDf)))) %>% mutate(targetCat=ifelse(type.to %in% CXOuts,"CX output links","Other CX links"),
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
    #fullDf$linkCat[fullDf$motif=="Reciprocal"] = fullDf$sourceCat[fullDf$motif=="Reciprocal"]
    
    return(fullDf %>% mutate(linkCat = factor(linkCat,levels=c("out of CX pathways","CX parallel connection","CX canonical feedback","linked targets in CX")),
                             roi=factor(roi,levels=c("Outside","EB","FB","PB","NO(R)","NO(L)")),
                             roiCat=ifelse(roi=="Outside","Outside","CX")))}
}

plotMotifs <- function(graphDf){
  focusContrib <- filter(graphDf,roi=="Outside") %>% group_by(type.to) %>% summarize(weightRelative=sum(weightRelative)) %>% ungroup()
  motifGG <- makeGraph(graphDf) 
  edgePal <- c("Plum",paletteer_d("ggthemes::Traffic")[1:3])
  names(edgePal) <- c("out of CX pathways","CX parallel connection","CX canonical feedback","linked targets in CX")
  
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
           theme_paper_map() + scale_edge_color_manual(name="motif",values=edgePal,breaks=c("out of CX pathways","CX parallel connection","CX canonical feedback","linked targets in CX"),drop=FALSE)+
           scale_edge_width(name="relative weight/pathway weight",range=c(0.2,3),limits=c(0,1),)+
           geom_node_text(aes(label=name),nudge_y = 0.15) + 
           coord_fixed(clip="off"))
  #}
}

#geom_node_text(aes(label=name),size=1.5,nudge_y = 0.06,hjust="inward") + 
