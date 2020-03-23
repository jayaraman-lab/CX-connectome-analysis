################### A collection of functions for loading, processing, and plotting FB synapse distributions ##################



###############################################################################################################################
################## Function for loading FB synapses by layer ##################################################################



Get_FBlayer_ColumnarSyns <- function(BodyIDs, ROI, Layer) {
  
  
  SynLocs =  neuprint_get_synapses(BodyIDs, roi = ROI)
  
  SynLocs =  mutate(SynLocs, type=neuprint_get_meta(bodyid)$type, partner_type=neuprint_get_meta(partner)$type,
                    x=as.numeric(x),y=as.numeric(y),z=as.numeric(z))
  
  SynLocs$prepost[SynLocs$prepost==0]="Output"
  SynLocs$prepost[SynLocs$prepost==1]="Input"
  
  SynLocs$Layer=Layer
  return(SynLocs)
}



###############################################################################################################################
################# Function for plotting synapse distribution for all FB layers ################################################



PlotSyns_byLayer <- function(Plot_Syns, Type, PlotName, Grouping, DIR, TOPLOT_TF, TOPLOTLEGEND_TF){
  
  P1 = Plot_Layer(Plot_Syns, "L1", Outline_L1_XZ, paste("L1 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P2 = Plot_Layer(Plot_Syns, "L2", Outline_L2_XZ, paste("L2 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P3 = Plot_Layer(Plot_Syns, "L3", Outline_L3_XZ, paste("L3 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P4 = Plot_Layer(Plot_Syns, "L4", Outline_L4_XZ, paste("L4 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P5 = Plot_Layer(Plot_Syns, "L5", Outline_L5_XZ, paste("L5 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P6 = Plot_Layer(Plot_Syns, "L6", Outline_L6_XZ, paste("L6 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P7 = Plot_Layer(Plot_Syns, "L7", Outline_L7_XZ, paste("L7 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P8 = Plot_Layer(Plot_Syns, "L8", Outline_L8_XZ, paste("L8 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P9 = Plot_Layer(Plot_Syns, "L9", Outline_L9_XZ, paste("L9 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  
  P12<-ggarrange(plotlist=list(P9,P8,P7,P6,P5,P4,P3,P2,P1), ncol =3, nrow = 3)
  ggsave(paste(DIR, PlotName,".png",sep=""),
         plot = P12, device='png', scale = 1, width = 12, height = 10, units ="in", dpi = 500, limitsize = TRUE)
  
}



###############################################################################################################################
################# Function for plotting synapse distribution of each individual FB layer ######################################



Plot_Layer <- function(Plot_Syns, Layer_To_Plot, Outline, Title, Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF){
  
  # Get data
  PlotData=subset(Plot_Syns, Layer==Layer_To_Plot)
  
  # Get color map
  col_vector=GetColorMap(Grouping)
  
  # Get color variable
  if (Grouping == "FBcol"){
    ColorVar=PlotData$FBcol
  } else if (Grouping == "PBglom") {
    ColorVar=PlotData$PBglom
  }
  
  # Make plot
  P1 <- Plot_Synapse_Distribution(PlotData, "XZ", Outline, ColorVar, col_vector, DIR, PlotName, Layer_To_Plot, TOPLOT_TF, TOPLOTLEGEND_TF)

  return(P1)
}



###############################################################################################################################
################# Generic plotting function for synapse distributions  ################################################



Plot_Synapse_Distribution <- function(PlotData, VIEW, Outline, ColorVar, col_vector, DIR, PlotName, PlotName2, TOPLOT_TF, TOPLOTLEGEND_TF, Size=1){
  
  # Get view
  if (VIEW=="XY"){
    PlotData = PlotData[!colnames(PlotData)=="Z"]
    colnames(PlotData)[colnames(PlotData)=="Y"]="Z"
    PlotData$Z=-PlotData$Z
  }
  
  # Plot just the data for now
  P1 <- ggplot() + geom_point(data=PlotData, aes(x=X, y=Z, color=ColorVar) ,  size=Size, alpha = 1, stroke = 0, shape=16) + 
    geom_path(data=Outline, aes(x=c1, y=c2), size = 2) + xlim(c(-9000, 9000)) +  coord_fixed(ratio = 1) 

  # Save with scaled visible if TOPLOTLEGEND_TF us TRUE
  if (TOPLOTLEGEND_TF==TRUE){
    ggsave(paste(DIR, PlotName,"_",PlotName2 , "_Legend.png",sep=""), plot = P1, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE) 
  }  
  
  # Format the rest of the plot, get rid of axes
  P1 <- P1 + theme_void() + guides(fill=FALSE) + theme(legend.position = "none") +
      scale_color_manual(values=col_vector)  +  theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )
  
  # Save plot
  if (TOPLOT_TF==TRUE){
    ggsave(paste(DIR, PlotName,"_",PlotName2 , ".png",sep=""), plot = P1, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE,  bg = "transparent")
  } 
  
  return(P1)
}



###############################################################################################################################
################# Function for getting color map of FB columns or PB glomeruli ################################################



GetColorMap <- function(Grouping){
  if (Grouping == "FBcol"){
    col_vector=brewer.pal(8, "Paired")
    col_vector[9]=col_vector[1]
    col_vector=col_vector[c(1,4,7,2,5,8,3,6,9)]
  } else if (Grouping == "PBglom") {
    L_col_vector=brewer.pal(8, "Paired")
    L_col_vector[9]=L_col_vector[1]
    L_col_vector=L_col_vector[c(1,4,7,2,5,8,3,6,9)]
    R_col_vector=rev(L_col_vector)
    col_vector=c(L_col_vector, R_col_vector)
  }
  return(col_vector)
}



###############################################################################################################################
################# Function for plotting synapse distribution for each FB layer ################################################



PlotSyns_byColumn <- function(Plot_Syns, Title, PlotName, DIR, TOPLOT_TF, TOPLOTLEGEND_TF){
  
  # Get color vector
  col_vector=GetColorMap("FBcol")
  
  ### XY View
  P1 <- Plot_Synapse_Distribution(Plot_Syns, "XY", Outline_FB_XY, Plot_Syns$FBcol, col_vector, DIR, PlotName, "XY", TOPLOT_TF, TOPLOTLEGEND_TF)
  
  ### XZ View
  P1 <- Plot_Synapse_Distribution(Plot_Syns, "XZ", Outline_FB_XZ, Plot_Syns$FBcol, col_vector, DIR, PlotName, "XZ", TOPLOT_TF, TOPLOTLEGEND_TF)
  
}


###############################################################################################################################
####### Functions for plotting column positions as defined by Fx neuron synapses locations ####################################


PlotColumn_Types<- function(Plot_Syns, PlotName, DIR, PLOT){
  Types=unique(Plot_Syns$type)
  for (nnn in 1:length(Types)){
    Temp_Type=subset(Plot_Syns, type == as.character(Types[nnn]))
    PlotColumn_Aves(Temp_Type, paste(PlotName, "_", Types[nnn],sep="") , DIR, PLOT)
  }
}



PlotColumn_Aves <- function(Plot_Syns, PlotName, DIR, PLOT){
  
  # get layers
  Layers_All=sort(unique(as.character(Plot_Syns$Layer)))
  
  # Make plots
  l1=Plot_Crosses(Plot_Syns, 1, DIR, PlotName, PLOT)
  l2=Plot_Crosses(Plot_Syns, 2, DIR, PlotName, PLOT)
  l3=Plot_Crosses(Plot_Syns, 3, DIR, PlotName, PLOT)
  l4=Plot_Crosses(Plot_Syns, 4, DIR, PlotName, PLOT)
  l5=Plot_Crosses(Plot_Syns, 5, DIR, PlotName, PLOT)
  l6=Plot_Crosses(Plot_Syns, 6, DIR, PlotName, PLOT)
  l7=Plot_Crosses(Plot_Syns, 7, DIR, PlotName, PLOT)
  l8=Plot_Crosses(Plot_Syns, 8, DIR, PlotName, PLOT)
  
  
  # Save plot
  PP=ggarrange(l8,l7,l6, l5, l4, l3, l2, l1,  nrow = 3, ncol = 3)
  ggsave(paste(DIR, PlotName, "_AllLayer" , ".png",sep=""), plot = PP, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE,  bg = "transparent") 

}


Plot_Crosses <- function(Plot_Syns, LayerToPlot, DIR, PlotName, PLOT){
  
  # Get color vector
  col_vector=GetColorMap("FBcol")
  
  # Get mean positions for this layer
  Temp_Dist=subset(Plot_Syns, Layer == paste("L",as.character(LayerToPlot),sep="")) 
  colnames(Temp_Dist)[which(colnames(Temp_Dist) %in% c("X_Mean","Z_Mean"))] = c("X","Z")
  
  # Get layer outline 
  eval(parse( text= paste("Outline=Outline_L", as.character(LayerToPlot),"_XZ",sep="") ))
  
  # Plot data
  p1 <- ggplot() + geom_path(data=Outline, aes(x=c1, y=c2), size = 1) + xlim(c(-9000, 9000)) +  coord_fixed(ratio = 1) 
  
  if (length(Temp_Dist$FBcol)>0){
    P_All=list()
    for (nnn in 1:length(Temp_Dist$FBcol) ){
      Temp_Neuron_Dist=Temp_Dist[nnn,]
      CCC= col_vector[which(levels(Temp_Neuron_Dist$FBcol) == Temp_Neuron_Dist$FBcol[1])]
      
      p1 = p1 +  
        geom_line(data=data.frame(XXX = c(Temp_Neuron_Dist$Cross_XL_X, Temp_Neuron_Dist$Cross_XR_X),
                                  YYY=c(Temp_Neuron_Dist$Cross_XL_Z, Temp_Neuron_Dist$Cross_XR_Z)), aes(x=XXX, y=YYY), size = 0.5, color= CCC ) +
        geom_line(data=data.frame(XXX = c(Temp_Neuron_Dist$Cross_YU_X, Temp_Neuron_Dist$Cross_YD_X),
                                  YYY=c(Temp_Neuron_Dist$Cross_YU_Z, Temp_Neuron_Dist$Cross_YD_Z) ), aes(x=XXX, y=YYY), size = 0.5 , color= CCC) + 
        geom_point(data=Temp_Neuron_Dist, aes(x=X,y=Z),size=2, color= CCC)
    }
    
  }
  
  
  p1 = p1 + theme_void() + guides(fill=FALSE) + theme(legend.position = "none") +
    scale_color_manual(values=col_vector)  +  theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )
  
  if (PLOT==TRUE){
  ggsave(paste(DIR, PlotName, "_L", as.character(LayerToPlot) , ".png",sep=""), plot = p1, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE,  bg = "transparent")
  }
  return(p1)
}



###############################################################################################################################
### Third of three functions for assigning FB columns to PB-FB-XX and D0 neurons, here based on closest type average  #########


Assign_FB_Columns <- function(N_Distribution, FX_Column_Positions_Filt, DIR){
  
  N_Distribution$Distance2FxCol=NA
  N_Distribution$FBcol=as.character(N_Distribution$FBcol)
  Column_Names=colnames(N_Distribution)
  
  # Get unique body IDs
  Uni_BodyID=unique(N_Distribution$bodyid)
  
  # Dataframe mapping bodyid to FBcol
  Bodyid_To_FBcol=data.frame(bodyid=numeric(), FBcol=as.character(), Distance2FxCol=numeric(), FBcol_Con=numeric())
  
  
  # Loop over each neuron
  for (nnn in 1:length(Uni_BodyID)){
    
    # Get all synapses from all layers for this neuron
    Neuron_Data=subset(N_Distribution, bodyid==Uni_BodyID[nnn])
    
    # Get the synapses for this neuron from the layer with the most synapses for this neuron type
    MaxLayer=Neuron_Data$Layer[which(Neuron_Data$NumberOfSyns == max(Neuron_Data$NumberOfSyns))]
    Layer_Data=subset(Neuron_Data, Layer==MaxLayer)
    
    # Subset FX data based on this layer
    TempCol_Positions=subset(FX_Column_Positions_Filt, Layer==Layer_Data$Layer)
    
    # Compute distance to every FXs neurons synapse mean in this layer
    Temp_Distances_All= data.frame(Col= TempCol_Positions$FBcol, Distance= ((TempCol_Positions$X_Mean - Layer_Data$X_Mean)^2 +
                                                                              (TempCol_Positions$Z_Mean - Layer_Data$Z_Mean)^2 +
                                                                              (TempCol_Positions$Y_Mean - Layer_Data$Y_Mean)^2)^0.5)  # could do some normalization here
    
    # Convert column number to numeric
    Temp_Distances_All$Col=as.numeric(sapply(Temp_Distances_All$Col, substring, 2, 2))
    
    # Compute average distance by column
    Temp_Distances_Mean=Temp_Distances_All %>% group_by(Col) %>% summarise(Distance=mean(Distance)) #quantile(Distance,0.1)
    
    # Interpolate data
    #Temp_Distances_Mean_UP=as.data.frame(approx(Temp_Distances_Mean$Col, y = Temp_Distances_Mean$Distance, xout=seq(1,9,0.1), method="linear"))
    Temp_Distances_Mean_UP=as.data.frame(spline(x = Temp_Distances_Mean$Col, y = Temp_Distances_Mean$Distance,  n = 91, method = "fmm", xmin = 0.6, xmax = 9.4))
    colnames(Temp_Distances_Mean_UP)=c("Col","Distance")
    
    # Take min of function
    Min_Distance=min(Temp_Distances_Mean_UP$Distance)
    Min_IND= which(Temp_Distances_Mean_UP$Distance == Min_Distance)
    Min_C=Temp_Distances_Mean_UP$Col[Min_IND]
    Min_Col=round(Min_C)
    
    # Take min distance as closest column
    Temp_df=data.frame(bodyid=Uni_BodyID[nnn], FBcol = paste("C",as.character(Min_Col),sep=""), Distance2FxCol= Min_Distance, FBcol_Con=Min_C)
    Bodyid_To_FBcol=rbind(Bodyid_To_FBcol, Temp_df)
    
    
    if (Neuron_Data$type[1] %in% c("PFGs","PFL1","PFL2","PFL3")){
      P1<-ggplot() + geom_point(data=Temp_Distances_All, aes(x=Col, y=Distance), position = position_jitter(width = 0.1), size=0.2) +
        geom_point(data=Temp_Distances_Mean, aes(x=Col, y=Distance),color="red") +
        ggtitle(paste(Neuron_Data$type[1], " ", Neuron_Data$PBglom[1], " ", Neuron_Data$bodyid[1] )) +
        geom_line(data=Temp_Distances_Mean_UP, aes(x=Col,y=Distance)) + scale_x_continuous(breaks=seq(1,9,1)) + 
        geom_point(data=data.frame(XXX = Min_C, YYY = Min_Distance), aes(x=XXX,y=YYY), color="blue", size=2) + 
        geom_point(data=data.frame(XXX = Min_Col, YYY = Min_Distance), aes(x=XXX,y=YYY), color="green", size=2)
      
      
      OUTLINE = eval(parse( text= paste("Outline=Outline_", as.character(MaxLayer),"_XZ",sep="") ))
      col_vector=GetColorMap("FBcol")
      
      P2<-ggplot() + geom_point(data=Layer_Data, aes(x=X_Mean, y=Z_Mean),  size=5, alpha = 1, stroke = 0, shape=16) +
        geom_point(data=TempCol_Positions, aes(x=X_Mean, y=Z_Mean, color=FBcol),  size=2, alpha = 1, stroke = 0, shape=16) + 
        geom_path(data=OUTLINE, aes(x=c1, y=c2), size = 1) + xlim(c(-9000, 9000)) +  coord_fixed(ratio = 1) + theme_void() + guides(fill=FALSE) + theme(legend.position = "none") +
        scale_color_manual(values=col_vector)  +  theme(
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
        )  +
        ggtitle( paste("Layer ", as.character(MaxLayer),  "   PBglom ", as.character(Neuron_Data$PBglom[1]) ,
                       "   FBcol ", as.character(Temp_df$FBcol[1]), "   ", as.character(Neuron_Data$type[1]), "   ",  as.character(Neuron_Data$bodyid[1]),sep="") )
      
      PP=ggarrange(P1, P2, ncol=2, nrow=1)
      
      ggsave(paste(DIR, Neuron_Data$type[1], "_", Neuron_Data$PBglom[1], "_" ,Neuron_Data$bodyid[1] ,".png",sep=""),
             plot = PP, device='png', scale = 1, width = 10, height = 5, units ="in", dpi = 300, limitsize = TRUE)
      
    }
    
  }
  
  N_Distribution=N_Distribution[ , !colnames(N_Distribution) %in% c("FBcol","Distance2FxCol")]
  N_Distribution_New=merge(N_Distribution, Bodyid_To_FBcol, by="bodyid")
  N_Distribution_New=N_Distribution_New[c(Column_Names,"FBcol_Con")]
  
  N_Distribution_New$FBcol=factor(N_Distribution_New$FBcol, levels = c("C1","C2","C3","C4","C5","C6","C7","C8","C9") )
  
  return(N_Distribution_New)
}

###############################################################################################################################
### Function for plotting mapping between PB glom and FB cols #################################################################



Plot_PBglom_FBcol_Mapping <- function(PFX_GlomToColumn, DIR, NAME){
  
  
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
    col_vector=GetColorMap("FBcol")
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
    
    
    
  }
}


###############################################################################################################################
### Function for getting the width of each neuron's synapse cloud in a direction tangent to the FB layer ######################


Get_SynLayerDistribution <- function(All_Neurons, Thresh, PlotDir){
  
  
  Syn_Distribution <- data.frame(bodyid=as.numeric(), type=as.character(), Layer=as.character(), 
                                 Side=as.character(), PBglom=as.character(), FBcol=as.character(),
                                 X_Mean=as.numeric(), Y_Mean=as.numeric(), Z_Mean=as.numeric(), 
                                 X_STD_Orth=as.numeric(), Z_STD_Orth=as.numeric(),
                                 X_HW_Orth=as.numeric(), Z_HW_Orth=as.numeric(),
                                 Cross_XL_X=as.numeric(),  Cross_XL_Z=as.numeric(),
                                 Cross_XR_X=as.numeric(),  Cross_XR_Z=as.numeric(),
                                 Cross_YU_X=as.numeric(),  Cross_YU_Z=as.numeric(),
                                 Cross_YD_X=as.numeric(),  Cross_YD_Z=as.numeric(),
                                 Box_UL_X=as.numeric(),    Box_UL_Z=as.numeric(),
                                 Box_UR_X=as.numeric(),    Box_UR_Z=as.numeric(),
                                 Box_LL_X=as.numeric(),    Box_LL_Z=as.numeric(),
                                 Box_LR_X=as.numeric(),    Box_LR_Z=as.numeric(),
                                 LayerLength=as.numeric(), NumberOfSyns=as.numeric())
  
  
  Neuron_Types=unique(All_Neurons$type)
  Layers_All=sort(unique(All_Neurons$Layer))
  
  # Loop over neuron types
  for (ttt in 1:length(Neuron_Types)){
    SubData=subset(All_Neurons, type == Neuron_Types[ttt])
    
    
    # Loop of layers
    for (lll in 1:8){
      LayerData=subset(SubData, Layer==Layers_All[lll])
      eval(parse( text= paste("MID=Mid_L", as.character(lll),sep="") ))
      eval(parse( text= paste("OUTLINE=Outline_L", as.character(lll),"_XZ",sep="") ))
      
      
      # Loop over individual neurons
      NeuronIDs=unique(LayerData$bodyid)
      for (nnn in 1:length(NeuronIDs)){
        NeuronData=subset(LayerData, bodyid==NeuronIDs[nnn])
        
        # Rename Z as Y or it gets confusing when looking at plots
        Mean_YY=mean(NeuronData$Y)
        NeuronData=NeuronData[colnames(NeuronData) != "Y"]
        colnames(NeuronData)[colnames(NeuronData) == "Z"] <-"Y"
        
        if (length(NeuronData$X)>Thresh){
          
        ############# Project data into new coordinate frame whose main axis is locally parallel to FB layer ############# 
          
          # Get X range
          X_RangeMin = quantile(NeuronData$X, probs = 0.1) 
          X_RangeMax = quantile(NeuronData$X, probs = 0.9) 
          
          # Get slope of mid range 
          Mid_subset=subset(MID, X>X_RangeMin & X<X_RangeMax)
          
          # If no values in range, use two closest values
          if (length(Mid_subset$X)<2){
            Center=mean(c(X_RangeMin,X_RangeMax))
            INDS= order(abs(MID$X - Center))[c(1,2)]
            Mid_subset=MID[INDS,]
          }
          
          # Compute slope
          TempSlope=median( diff(Mid_subset$Y)/diff(Mid_subset$X) )
          
          # Define new axis (vectors should be length 1 and orthogonal to keep shape of point cloud)
          # Also, keep Ax_1 pointing towards right and Ax_2 pointing up, so that axes are consistent across clouds
          Ax_1=c(1,TempSlope)
          Ax_1=Ax_1/sum(Ax_1^2)^0.5

          Ax_2=c(1, -Ax_1[1]/Ax_1[2])
          Ax_2=Ax_2/sum(Ax_2^2)^0.5
          if (Ax_2[2]<0){
            Ax_2[1]=-Ax_2[1]
            Ax_2[2]=-Ax_2[2]} 

          
          # Make sure vectors are orthogonal
          if (Ax_1 %*% Ax_2 > 10^-15){print(paste("Vectors not orthogonal ttt=", as.character(ttt), "   lll = ", as.character(lll), "   nnn = ", as.character(nnn)))}
          
          # Project onto new x and y axes
          NeuronData$X_New= as.matrix(NeuronData[c("X","Y")]) %*% Ax_1  # This
          NeuronData$Y_New= as.matrix(NeuronData[c("X","Y")]) %*% Ax_2  # This
          

          # Define vectors to get from new space back to old space
          Ax1_NewToOld= t(Ax_1 * as.matrix(c(1,-1)))
          Ax2_NewToOld= t(Ax_2 * as.matrix(c(-1,1)))
          
          
        ############# Compute values we want to save for this neuron ############# 
          
          # Mean in original space 
          TempMean=as.data.frame( t(colMeans( data.matrix(NeuronData[c("X","Y")])) )) # This
          
          # Mean and STD along new axes
          TempMean_New=as.data.frame(t(colMeans(NeuronData[c("X_New","Y_New")])))
          TempSTD_New=as.data.frame(t(sapply(NeuronData[c("X_New","Y_New")], sd)))
          
          # Ranges along new axes
          TempRange_New=data.frame(X = abs(quantile(NeuronData$X_New, probs = 0.02) - quantile(NeuronData$X_New, probs = 0.98))/2,
                                   Y = abs(quantile(NeuronData$Y_New, probs = 0.02) - quantile(NeuronData$Y_New, probs = 0.98))/2)
          
          
          # Compute Cross points (4 points making up a cross around synapse mean position)
          Cross_XL= data.frame(X = as.numeric(as.matrix((TempMean_New -  c(TempRange_New$X,0))) %*% t(Ax1_NewToOld)) , 
                               Y = as.numeric(as.matrix((TempMean_New -  c(TempRange_New$X,0))) %*% t(Ax2_NewToOld)) )
          Cross_XR= data.frame(X = as.numeric(as.matrix((TempMean_New +  c(TempRange_New$X,0))) %*% t(Ax1_NewToOld)) , 
                               Y = as.numeric(as.matrix((TempMean_New +  c(TempRange_New$X,0))) %*% t(Ax2_NewToOld)) )
          Cross_YU= data.frame(X = as.numeric(as.matrix((TempMean_New +  c(0,TempRange_New$Y))) %*% t(Ax1_NewToOld)) , 
                               Y = as.numeric(as.matrix((TempMean_New +  c(0,TempRange_New$Y))) %*% t(Ax2_NewToOld)) )
          Cross_YD= data.frame(X = as.numeric(as.matrix((TempMean_New -  c(0,TempRange_New$Y))) %*% t(Ax1_NewToOld)) , 
                               Y = as.numeric(as.matrix((TempMean_New -  c(0,TempRange_New$Y))) %*% t(Ax2_NewToOld)) )
          
          # Compute box points (4 points making up a box around synapse mean position)
          Box_UL= data.frame(X = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(-1,1) )) %*% t(Ax1_NewToOld)) , 
                             Y = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(-1,1) )) %*% t(Ax2_NewToOld)) )
          Box_UR= data.frame(X = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(1,1) )) %*% t(Ax1_NewToOld)) , 
                             Y = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(1,1) )) %*% t(Ax2_NewToOld)) )
          Box_LL= data.frame(X = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(-1,-1) )) %*% t(Ax1_NewToOld)) , 
                             Y = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(-1,-1) )) %*% t(Ax2_NewToOld)) )
          Box_LR= data.frame(X = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(1,-1) )) %*% t(Ax1_NewToOld)) , 
                             Y = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(1,-1) )) %*% t(Ax2_NewToOld)) )
          
          
          
        ############# Plot some random examples ############# 
          
          # Debugging plot
          PLOT=0
          if ( mod(ttt*lll*nnn,50) ==0 | abs(TempMean$X)<150){
            p1<-ggplot() +  geom_path(data=OUTLINE, aes(x=c1, y=c2), size = 1, color="red") +
              geom_point(data=NeuronData, aes(x=X,y=Y)) + 
              geom_path(data=MID, aes(x=X, y=Y), size = 1, color="orange") + 
              geom_path(data=Mid_subset, aes(x=X, y=Y), size = 1, color="green") + 
              geom_point(data=TempMean, aes(x=X,y=Y,size=50),color="red") + 
              geom_point(data=Cross_XL, aes(x=X,y=Y,size=50),color="cornflowerblue") + 
              geom_point(data=Cross_XR, aes(x=X,y=Y,size=50),color="blueviolet") + 
              geom_point(data=Cross_YU, aes(x=X,y=Y,size=50),color="darkred") + 
              geom_point(data=Cross_YD, aes(x=X,y=Y,size=50),color="gold") + 
              geom_point(data=Box_UL, aes(x=X,y=Y,size=50),color="darkslateblue")  + 
              geom_point(data=Box_UR, aes(x=X,y=Y,size=50),color="darkorange") + 
              geom_point(data=Box_LL, aes(x=X,y=Y,size=50),color="deeppink") +
              geom_point(data=Box_LR, aes(x=X,y=Y,size=50),color="darkslategray") + coord_fixed(ratio = 1) +
              ggtitle(paste(NeuronData$type[1],"_Layer_",as.character(lll),sep=""))
             
            ggsave(paste(PlotDir, "Cloud_", NeuronData$type[1],"_Layer_",as.character(lll),"_", as.character(ttt*lll*nnn),".png",sep=""),
                   plot = p1, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 300, limitsize = TRUE)
            
            #print(paste("Plotting ttt=", as.character(ttt), "   lll = ", as.character(lll), "   nnn = ", as.character(nnn)))
            
            # p2<-ggplot() + geom_point(data=NeuronData, aes(x=X_New,y=Y_New))
            
            #ggarrange(p1,p2, nrow = 1, ncol = 2)
            
            remove(p1)
          
            }
          
          
        ############# Add data to dataframe ############# 
          
          TempDF <- data.frame(bodyid=NeuronData$bodyid[1], type=NeuronData$type[1], Layer=NeuronData$Layer[1], 
                               Side=NeuronData$Side[1], PBglom=NeuronData$PBglom[1], FBcol=NeuronData$FBcol[1],
                                         X_Mean=TempMean$X, Y_Mean=Mean_YY, Z_Mean=TempMean$Y, 
                                         X_STD_Orth=TempSTD_New$X_New, Z_STD_Orth=TempSTD_New$Y_New,
                                         X_HW_Orth=TempRange_New$X, Z_HW_Orth=TempRange_New$Y,
                                         Cross_XL_X=Cross_XL$X,    Cross_XL_Z=Cross_XL$Y,
                                         Cross_XR_X=Cross_XR$X,    Cross_XR_Z=Cross_XR$Y,
                                         Cross_YU_X=Cross_YU$X,    Cross_YU_Z=Cross_YU$Y,
                                         Cross_YD_X=Cross_YD$X,    Cross_YD_Z=Cross_YD$Y,
                                         Box_UL_X=Box_UL$X,        Box_UL_Z=Box_UL$Y,
                                         Box_UR_X=Box_UR$X,        Box_UR_Z=Box_UR$Y,
                                         Box_LL_X=Box_LL$X,        Box_LL_Z=Box_LL$Y,
                                         Box_LR_X=Box_LR$X,        Box_LR_Z=Box_LR$Y,
                                         LayerLength = sum((diff(MID$X)^2 + diff(MID$Y)^2)^0.5),
                                         NumberOfSyns = length(NeuronData$X) )

          Syn_Distribution=rbind(Syn_Distribution,TempDF)
          
        }
      }
    }
  }
  return(Syn_Distribution)
}

