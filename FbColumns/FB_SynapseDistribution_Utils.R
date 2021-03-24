##### A collection of functions for loading, processing, and plotting FB synapse distributions ##################



###############################################################################################################################
################# Functions for loading raw synapse locations #################################################################

Get_FBlayer_ColumnarSyns <- function(BodyIDs, ROI, Layer) {
  #' Function for loading FB synapses by layer, with associated partner type and other meta information.
  
  SynLocs =  neuprint_get_synapses(BodyIDs, roi = ROI)
  
  SynLocs =  mutate(SynLocs, type=neuprint_get_meta(bodyid)$type, partner_type=neuprint_get_meta(partner)$type,
                    x=as.numeric(x),y=as.numeric(y),z=as.numeric(z))
  
  SynLocs$prepost[SynLocs$prepost==0]="Output"
  SynLocs$prepost[SynLocs$prepost==1]="Input"
  
  SynLocs$Layer=Layer
  return(SynLocs)
}



###############################################################################################################################
################# Functions getting mean arbor locations for hDelta neurons and C0 vDelta neurons #############################


Kmeans_Synapses <- function(TempSynapseData){
  
  # Perform K-means clusters (n=2 clusters to get pre/post arbors)
  TempSynLocs=TempSynapseData[c("X","Y")]
  k2 <- kmeans(TempSynLocs, centers = 2)
  TempClusters=k2$cluster
  Centers=k2$centers
  
  # assign clusters to left or right half of FB
  Cluster_LR=data.frame(Cluster=c(1,2),LR=NA)
  if (Centers[1,1] > Centers[2,1]){
    Cluster_LR$LR=c("L","R")
  } else {
    Cluster_LR$LR=c("R","L")
  }
  TempSynapseData$LR[TempClusters==1]=Cluster_LR$LR[1]
  TempSynapseData$LR[TempClusters==2]=Cluster_LR$LR[2]
  
  return(TempSynapseData)
}


ComputeMeanArbor <- function(TempSynapseData, Kmeans_layer){
  
  # Calculate median position of L/R arbors and assign as input or output arbor (or neithers, in case of vDelta)
  TempColumnPositions= TempSynapseData %>% group_by(LR) %>% summarise(X=median(X),Y=median(Y), Z=median(Z), prepost=mean(prepost), numsyns=n())
  TempColumnPositions$bodyid=TempSynapseData$bodyid[1]
  TempColumnPositions$type=TempSynapseData$type[1]
  TempColumnPositions$name=TempSynapseData$name[1]
  TempColumnPositions$Layer=Kmeans_layer
  if (startsWith(TempSynapseData$type[1],"vDelta")){
    TempColumnPositions$prepost="Neither"
  } else if (startsWith(TempSynapseData$type[1],"hDelta")) {
    InputInd=which(TempColumnPositions$prepost == min(TempColumnPositions$prepost))
    if (length(InputInd)==1){
      TempColumnPositions$prepost[InputInd] = "Input"
      if (InputInd == 1){
        TempColumnPositions$prepost[2] = "Output"
      } else {
        TempColumnPositions$prepost[1] = "Output"
      }
    } else { # if length of prepost min is 2, then both arbors are inputs
      TempColumnPositions$prepost="Input"}
  }
  
  return(TempColumnPositions)
  
}   



Plot_Kmeans <- function(TempSynapseData, TempColumnPositions, OUTLINE, CurrentLayer, PlotDirTemp){
  
  P1=ggplot() + geom_point(data=subset(TempSynapseData,LR=="L"), aes(x=X, y=Z), colour="midnightblue" ,  size=1, alpha = 0.1) + 
    geom_point(data=subset(TempSynapseData,LR=="R"), aes(x=X, y=Z), colour="black" ,  size=1, alpha = 0.05) +
    geom_point(data=(TempColumnPositions), aes(x=X, y=Z, color=prepost),  size=5, alpha = 1) +
    coord_fixed(ratio = 1) +  geom_path(data=OUTLINE, aes(x=c1, y=c2), size = 1) + 
    ggtitle(paste( TempSynapseData$type[1], " ", as.character(TempSynapseData$bodyid[1]), " ", CurrentLayer, sep="" )) +  theme_bw() +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  
  
  ggsave(paste(PlotDirTemp, TempSynapseData$type[1], " ", as.character(TempSynapseData$bodyid[1]), " ", CurrentLayer, ".png", sep=""),
         plot = P1, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 150, limitsize = TRUE)
  
}



GetColorPalette <- function(TempSynapseData, TempColNum){
  
  # get color palette
  if (TempColNum == 12 & startsWith(TempSynapseData$type[1],'hDelta')){
    Colors=color(c("#CD4F39", "#7D26CD", "#63B8FF", "#FF1493", "#9A749A", "#00868B","#CD4F39", "#7D26CD", "#63B8FF", "#FF1493", "#9A749A", "#00868B"))
  } else if (TempColNum == 8 & startsWith(TempSynapseData$type[1],'hDelta')){
    Colors=c("#CD4F39", "#7D26CD", "#63B8FF", "#FF1493", "#CD4F39", "#7D26CD", "#63B8FF", "#FF1493")
  } else if (TempColNum == 6 & startsWith(TempSynapseData$type[1],'hDelta')){
    Colors=c("#CD4F39","#7D26CD","#63B8FF","#CD4F39","#7D26CD","#63B8FF")
  } else if (startsWith(TempSynapseData$type[1],'vDelta')){
    Colors=color(c("Gray50","#CD4F39", "#7D26CD", "#63B8FF", "#FF1493", "#CD96CD", "#00868B", "#ADFF2F", "#0000FF", "#FFB90F"))
  } else if (startsWith(TempSynapseData$type[1],'PF') | startsWith(TempSynapseData$type[1],'FS') | 
             startsWith(TempSynapseData$type[1],'FC') | startsWith(TempSynapseData$type[1],'FR')){
    Colors=color(c("#CD4F39", "#7D26CD", "#63B8FF", "#FF1493", "#CD96CD", "#00868B", "#ADFF2F", "#0000FF", "#FFB90F"))
  }
  
  return(Colors)
}



GetColorFactor <- function(TempSynapseData, TempColNum){
  
  # get color palette
  if (TempColNum == 12 & startsWith(TempSynapseData$type[1],'hDelta')){
    TempSynapseData$FBcol=factor(TempSynapseData$FBcol, levels=c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12") )
  } else if (TempColNum == 8 & startsWith(TempSynapseData$type[1],'hDelta')){
    TempSynapseData$FBcol=factor(TempSynapseData$FBcol, levels=c("C1","C2","C3","C4","C5","C6","C7","C8") )
  } else if (TempColNum == 6 & startsWith(TempSynapseData$type[1],'hDelta')){
    TempSynapseData$FBcol=factor(TempSynapseData$FBcol, levels=c("C1","C2","C3","C4","C5","C6") )
  } else if (startsWith(TempSynapseData$type[1],'vDelta')){
    TempSynapseData$FBcol=factor(TempSynapseData$FBcol, levels=c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9") )
  } else if (startsWith(TempSynapseData$type[1],'PF') | startsWith(TempSynapseData$type[1],'FS') | 
             startsWith(TempSynapseData$type[1],'FC') | startsWith(TempSynapseData$type[1],'FR')){
    TempSynapseData$FBcol=factor(TempSynapseData$FBcol, levels=c("C1","C2","C3","C4","C5","C6","C7","C8","C9") )
  }
  
  return(TempSynapseData)
}



PlotColLocs <- function(TempSynapseData, TempTypeSynsAll,Colors, OUTLINE, CurrentLayer, FigDirPDF, FigDirPNG){
  
  # Conversion vector from pixels to microns
  Convert=8/1000 # 8 nm/pixel, divided by 1000 nm per um.
  
  P3=ggplot() + geom_point(data=TempSynapseData, aes(x=X*Convert, y=-Z*Convert, colour=FBcol, shape=prepost),  size=4, alpha = 0.75) + 
    coord_fixed(ratio = 1) +  geom_path(data=OUTLINE, aes(x=c1*Convert, y=-c2*Convert), size = 0.5) + 
    ggtitle(paste(TempSynapseData$type[1], " ", CurrentLayer , "   Neurons = ", length(unique(TempSynapseData$bodyid)) ,
                  " of ", as.character(length(unique(TempTypeSynsAll$bodyid))) ,
                  "   Total Syns = ", sum(TempSynapseData$numsyns),  sep=""))  +
    scale_color_manual(values=Colors, drop=FALSE) + xlim(-9000*Convert,9000*Convert) + ylim(-6000*Convert, 6000*Convert) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  ggsave(paste(FigDirPDF, "ColumnLocs_", TempSynapseData$type[1], "_", as.character(CurrentLayer), ".pdf", sep=""),
         plot = P3, device='pdf', scale = 1, width = 10, height = 5, units ="in", dpi = 500, limitsize = TRUE)
  
  ggsave(paste(FigDirPNG, "ColumnLocs_", TempSynapseData$type[1], "_", as.character(CurrentLayer), ".png", sep=""),
         plot = P3, device='png', scale = 1, width = 10, height = 5, units ="in", dpi = 200, limitsize = TRUE)
  
  
}



###############################################################################################################################
################# Functions for plotting mapping between PB glom and FB cols ##################################################



GetColorMap <- function(Grouping){
  #' Function for getting color map of FB columns or PB glomeruli
  #' 
  if (Grouping == "FBcol_Old"){
    col_vector=brewer.pal(8, "Paired")
    col_vector[9]=col_vector[1]
    col_vector=col_vector[c(1,4,7,2,5,8,3,6,9)]
  } else if (Grouping == "FBcol"){
    col_vector=color(c("#CD4F39", "#7D26CD", "#63B8FF", "#FF1493", "#CD96CD", "#00868B", "#ADFF2F", "#0000FF", "#FFB90F"))
  } else if (Grouping == "FBcol_Function"){
    col_vector=color(c("#FFF688", "#FFA010", "#FF3610", "#FBB5DE", "#7D24E7", "#5C89C7", "#C7FFEF", "#81FB35", "#FFF688"))
  } else if (Grouping == "PBglom") {
    L_col_vector=brewer.pal(8, "Paired")
    L_col_vector[9]=L_col_vector[1]
    L_col_vector=L_col_vector[c(1,4,7,2,5,8,3,6,9)]
    R_col_vector=rev(L_col_vector)
    col_vector=c(L_col_vector, R_col_vector)
  }
  return(col_vector)
}



Plot_PBglom_FBcol_Mapping <- function(PFX_Distribution, DIR, NAME){
  
  # Get number of neurons going from each glom to each column
  PFX_GlomToColumn=PFX_Distribution %>% group_by(bodyid, type, PBglom, FBcol) %>% summarise(Num_layers = n(), FBcol_Con= mean(FBcol_Con))
  PFX_GlomToColumn=PFX_GlomToColumn %>% group_by(type, PBglom, FBcol) %>% summarise(Num_Neurons = n(), FBcol_Con= mean(FBcol_Con))
  
  
  PB_FB_Types=unique(PFX_GlomToColumn$type)
  for (nnn in 1:length(PB_FB_Types)){
    
    
    # Data for mapping
    PB_FB_Mapping=subset(PFX_GlomToColumn, type==PB_FB_Types[nnn] )
    
    
    # Set up node names and positions
    nodes = data.frame(name = c("L9","L8","L7","L6","L5","L4","L3","L2","L1",
                                "R1","R2","R3","R4","R5","R6","R7","R8","R9",
                                "C1","C2","C3","C4","C5","C6","C7","C8","C9") )
    nodes$x <- c( seq(from= 9 ,to = 1), seq(from= -1 ,to = -9),  seq(from= -4 ,to = 4))
    nodes$y <- c( rep(1, 18), rep(0, 9))
    
    
    # Assign node from/to pairs (numerical)
    PB_FB_Mapping$from = match(PB_FB_Mapping$PBglom, nodes$name)  
    PB_FB_Mapping$to   = match(PB_FB_Mapping$FBcol,  nodes$name)
    
    
    # build graph object
    graph = tbl_graph(nodes, PB_FB_Mapping)
    
    
    # Get column vector
    col_vector=GetColorMap("FBcol_Function")
    PB_FB_Mapping$Color = order(PB_FB_Mapping$FBcol)
    
    
    P1<-ggraph(graph,layout="manual",x=nodes$x,y=nodes$y) +
      geom_edge_diagonal(aes(width=Num_Neurons,color=FBcol),alpha=0.5,strength=0.5) +
      geom_node_point(size=5)  + scale_edge_color_manual(values=col_vector)  +
      geom_node_text(aes(label=name),angle=40,size=4, nudge_y = c(rep(0.06,18),rep(-0.06,9)) ) +
      theme_classic() + 
      theme(legend.text=element_text(size=6),legend.title=element_text(size=6),
            axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),axis.title.y=element_blank()) + 
      ggtitle(PB_FB_Mapping$type[1])
    
    ggsave(paste(DIR, "Graph_", PB_FB_Mapping$type[1], "_", NAME, ".png",sep=""), plot = P1, device='png', scale = 1, width = 8, height = 4.5, units ="in", dpi = 500, limitsize = TRUE) 
    ggsave(paste(DIR, "Graph_", PB_FB_Mapping$type[1], "_", NAME, ".pdf",sep=""), plot = P1, device='pdf', scale = 1, width = 8, height = 4.5, units ="in", dpi = 500, limitsize = TRUE) 
    
    
    PB_FB_Mapping$COL=as.numeric(sapply(PB_FB_Mapping$FBcol, substring, 2, 2))

    
    R_PB=subset(PB_FB_Mapping, startsWith(as.character(PBglom), "R"))
    L_PB=subset(PB_FB_Mapping, startsWith(as.character(PBglom), "L"))
    L_PB$PBglom=factor(L_PB$PBglom, levels = sort(unique(L_PB$PBglom), decreasing = TRUE))
    
    P2 <- ggplot() + geom_point(data=R_PB, aes(x=FBcol_Con, y=PBglom, size=1)) + 
      geom_point(data=R_PB, aes(x=COL, y=PBglom, size=1), color="blue") +
      scale_x_continuous(breaks=seq(from=1,to=9), limits=c(0,10)) 
    
    P3 <- ggplot() + geom_point(data=L_PB, aes(x=FBcol_Con, y=PBglom, size=1)) + 
      geom_point(data=L_PB, aes(x=COL, y=PBglom, size=1), color="blue") +
      scale_x_reverse(breaks=seq(from=1,to=9), limits=c(10,0)) 
      
    
    P12=ggarrange(P2,P3,ncol=1,nrow=2)
    P123=ggarrange(P1, P12, ncol=2, nrow=1)
    
    
    ggsave(paste(DIR, "Mapping_", PB_FB_Mapping$type[1], "_", NAME, ".png",sep=""), plot = P123, device='png', scale = 1, width = 12, height = 5, units ="in", dpi = 500, limitsize = TRUE) 
    ggsave(paste(DIR, "Mapping_", PB_FB_Mapping$type[1], "_", NAME, ".pdf",sep=""), plot = P123, device='pdf', scale = 1, width = 12, height = 5, units ="in", dpi = 500, limitsize = TRUE) 
    
    
  }
  return(PFX_GlomToColumn)
}


