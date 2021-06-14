##### A collection of functions analyzing FB synapse distributions, connectivity, etc.  ######################################


###############################################################################################################################
################# Function for loading raw synapse locations ##################################################################


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
################# Function for parsing columnar neuron names ##################################################################


Assign_FBcol_PBglom <- function(DF, name_field, side_field, PBglom_field, FBcol_field){
  #' Function for parsing neuron names to extract PB glomeruli and FB column info.
  
  # Assign side
  DF[[side_field]]=NA
  DF[[side_field]][grepl("_R",DF[[name_field]])]<-"R"
  DF[[side_field]][grepl("_L",DF[[name_field]])]<-"L"
  
  # Assign PB glomeruli
  DF[[PBglom_field]]=NA
  RightColumns=str_extract(DF[[name_field]], "_R(\\d+)")
  LeftColumns=str_extract(DF[[name_field]], "_L(\\d+)")
  DF[[PBglom_field]][which(!is.na(RightColumns))]=RightColumns[which(!is.na(RightColumns))]
  DF[[PBglom_field]][which(!is.na(LeftColumns))]=LeftColumns[which(!is.na(LeftColumns))]
  DF[[PBglom_field]]=sapply(DF[[PBglom_field]], substring, 2, 3)
  
  # Assign FB columns
  DF[[FBcol_field]]=NA
  DF[[FBcol_field]]=str_extract(DF[[name_field]], "_C(\\d+)")
  DF[[FBcol_field]]=sapply(DF[[FBcol_field]], substring, 2, 4)
  
  return(DF)
}


###############################################################################################################################
################# Functions getting median arbor locations for hDelta neurons and C0 vDelta neurons #############################


Kmeans_Synapses <- function(TempSynapseData){
  #' Function for performing K-means clustering to get the two arbor positions for hDelta and C0 vDelta neurons
  
  # Perform K-means clustering (n=2 clusters to get pre/post arbors)
  TempSynLocs=TempSynapseData[c("X","Y")]
  k2 <- kmeans(TempSynLocs, centers = 2)
  TempClusters=k2$cluster
  Centers=k2$centers
  
  # Assign clusters to left or right half of FB
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
  #' Function for calculating the median position of the left and right arbors for hDelta and C0 vDelta neurons

  # Calculate median column position
  TempColumnPositions= TempSynapseData %>% group_by(LR) %>% summarise(X=median(X),Y=median(Y), Z=median(Z), prepost=mean(prepost), numsyns=n())
  TempColumnPositions$bodyid=TempSynapseData$bodyid[1]
  TempColumnPositions$type=TempSynapseData$type[1]
  TempColumnPositions$name=TempSynapseData$name[1]
  TempColumnPositions$Layer=Kmeans_layer
  
  # Assign arbors as input or output (or neithers, in case of vDelta)
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
  #' Function for plotting the K-means clustered synapse locations to make sure arbors are being separated properly
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
  #' Function that gets columnar color palette according to each neuron type's columnar number
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
  #' Function that sets factor levels on FB columns according to neuron type
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


PlotColLocs <- function(TempSynapseData, TempTypeSynsAll,Colors, OUTLINE, CurrentLayer, FigDir){
  #' Function for plotting the median synapse locations of all neurons of a given neuron type
  
  # Conversion factor from pixels to microns
  Convert=8/1000 # 8 nm/pixel, divided by 1000 nm per um.
  
  # plot each neurons median synapse location
  P3=ggplot() + geom_point(data=TempSynapseData, aes(x=X*Convert, y=-Z*Convert, colour=FBcol, shape=prepost),  size=4, alpha = 0.75) + 
    coord_fixed(ratio = 1) +  geom_path(data=OUTLINE, aes(x=c1*Convert, y=-c2*Convert), size = 0.5) + 
    ggtitle(paste(TempSynapseData$type[1], " ", CurrentLayer , "   Neurons = ", length(unique(TempSynapseData$bodyid)) ,
                  " of ", as.character(length(unique(TempTypeSynsAll$bodyid))) ,
                  "   Total Syns = ", sum(TempSynapseData$numsyns),  sep=""))  +
    scale_color_manual(values=Colors, drop=FALSE) + xlim(-9000*Convert,9000*Convert) + ylim(-6000*Convert, 6000*Convert) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  ggsave(paste(FigDir, "ColumnLocs_", TempSynapseData$type[1], "_", as.character(CurrentLayer), ".png", sep=""),
         plot = P3, device='png', scale = 1, width = 10, height = 5, units ="in", dpi = 200, limitsize = TRUE)
}


###############################################################################################################################
################# Function for plotting mapping between PB glom and FB cols ###################################################


Plot_PBglom_FBcol_Mapping <- function(PFX_Distribution, DIR){
  #' Function that makes a graph showing the mapping from PB glomeruli to FB columns for PB-FB-XX neuron types
  
  # Get rid of neurons without assigned PB glomeruli (mostly neurons with "irreg" in their name and odd wiring patterns)
  PFX_Distribution=subset(PFX_Distribution, !is.na(PBglom))
  PFX_Distribution$PBglom=factor(PFX_Distribution$PBglom, levels=c("L9","L8","L7","L6","L5","L4","L3","L2","L1",
                                                                   "R1","R2","R3","R4","R5","R6","R7","R8","R9"))
  
  # Get number of neurons going from each PB glomerulus to each FB column
  PFX_GlomToColumn=PFX_Distribution %>% group_by(type, PBglom, FBcol) %>% summarise(Num_Neurons = n()) 
  
  # Loop over all PB-FB-XX neuron types and plot their PB-FB mapping
  PB_FB_Types=unique(PFX_GlomToColumn$type)
  for (nnn in 1:length(PB_FB_Types)){
    
    # Data for this neuron type
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
    
    # Build graph object
    graph = tbl_graph(nodes, PB_FB_Mapping)
    
    # Get color palette 
    col_vector=color(rev(c("#FFF688", "#FFA010", "#FF3610", "#FBB5DE", "#7D24E7", "#5C89C7", "#C7FFEF", "#81FB35", "#FFF688",
                       "#FFF688", "#FFA010", "#FF3610", "#FBB5DE", "#7D24E7", "#5C89C7", "#C7FFEF", "#81FB35", "#FFF688")))
    
    # Plot PB-FB mapping
    P1<-ggraph(graph,layout="manual",x=nodes$x,y=nodes$y) +
      geom_edge_diagonal(aes(width=Num_Neurons,color=PBglom),alpha=0.5,strength=0.5) +
      geom_node_point(size=5)  + scale_edge_color_manual(values=col_vector, drop=FALSE)  +
      geom_node_text(aes(label=name),angle=40,size=4, nudge_y = c(rep(0.06,18),rep(-0.06,9)) ) +
      theme_classic() + 
      theme(legend.text=element_text(size=6),legend.title=element_text(size=6),
            axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),axis.title.y=element_blank()) + 
      ggtitle(PB_FB_Mapping$type[1])
    ggsave(paste(DIR, "Graph_", PB_FB_Mapping$type[1], ".png",sep=""), plot = P1, device='png', scale = 1, width = 10, height = 5, units ="in", dpi = 500, limitsize = TRUE) 
  }
}


###############################################################################################################################
################# Functions for performing PCA on column-to-column connection matrix ##########################################



motifPCA <- function(PCA_In, Norm, Pre_Post){
  #' Function for doing PCA on Vectors
  
  # Binarize or normalize data
  if (Norm=="Binary"){
    PCA_In[PCA_In>0.001]=1 #0.001
  } 
  
  # Calculate covariance matrix across dimensions 
  cov = data.frame(cov(PCA_In))
  
  # Calculate eigenvectors and values (columns contain PCs, ranked from most variance accounted for to least)
  covEigen = eigen(cov) 
  
  # Project data into PC space
  NewData=PCA_In %*% covEigen$vectors
  
  # Make dataframe
  PCA_Out=cbind(Pre_Post,as.data.frame(NewData))
  
  # Compute variance explained by first two PCs
  PC_Cov=cov(NewData)
  PC_Var=diag(PC_Cov)
  PC12_Var=sum(PC_Var[c(1,2)])/sum(PC_Var)
  print(paste(as.character(PC12_Var*100), " % of varianced explained by first two PCs",sep=""))
  
  OUT=list(PCA_Out, covEigen)
  return(OUT)
}


GetPC <- function(covEigen,PCNUM){
  #' Function that converts a PC back into matrix form
  PC=covEigen$vectors[,PCNUM]
  PC=matrix(PC,nrow=10,ncol=10)
  rownames(PC)=colnames(PC)=paste("C",0:9,sep="")
  PC=melt(PC)
  return(PC)
}


Plot_ColCol_Matix <- function(PFLn_Input_Network){
  #' Function for plotting column-to-column connectivity matrix
  P1 <- ggplot(PFLn_Input_Network) +
    scale_fill_gradient2(low="thistle", mid="blueviolet", high="black", 
                         midpoint =0.5*max(PFLn_Input_Network$weightRelative), 
                         limits=c(0,max(PFLn_Input_Network$weightRelative)*1)) +
    geom_tile(aes(to_name,from_name,fill=weightRelative)) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(), 
          strip.placement = "outside", strip.background = element_rect(fill=NA, colour="grey50")) +
    facet_grid(reorder(type.from, desc(type.from)) ~ type.to, space="free", scales="free",switch="both")
  return(P1)
}


Plot_PC_Matrix <- function(PCx, PC_Var, PC_Ind){
  #' A function for plotting a principal component in matrix form
  PX=ggplot(PCx) + theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="blue", high="red", limits=c(-0.5,0.5), oob=squish) + coord_fixed() +
    geom_tile(aes(Var2,Var1,fill=value)) + ggtitle(paste("PC", as.character(PC_Ind), " ", round(PC_Var$PropVar[PC_Ind]*1000)/1000)) + theme(legend.position="none") +
    xlab("post column") + ylab("pre column")
  return(PX)
}

