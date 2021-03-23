##### A collection of functions for loading, processing, and plotting FB synapse distributions ##################


Get_FBlayer_ColumnarSyns <- function(BodyIDs, ROI, Layer) {
  #' Function for loading FB synapses by layer, with associated partner type and other meta information.
  #' @param BodyIDs A list of bodyIDs
  #' @param ROI A string for an ROI (e.g. "FBl1")
  #' @param Layer A string for an FB layer (e.g. ""FBl1")
  
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
  
  P1 = Plot_Layer(Plot_Syns, "L1", Outline_L1_XZ, paste("L1 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P2 = Plot_Layer(Plot_Syns, "L2", Outline_L2_XZ, paste("L2 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P3 = Plot_Layer(Plot_Syns, "L3", Outline_L3_XZ, paste("L3 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P4 = Plot_Layer(Plot_Syns, "L4", Outline_L4_XZ, paste("L4 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P5 = Plot_Layer(Plot_Syns, "L5", Outline_L5_XZ, paste("L5 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P6 = Plot_Layer(Plot_Syns, "L6", Outline_L6_XZ, paste("L6 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P7 = Plot_Layer(Plot_Syns, "L7", Outline_L7_XZ, paste("L7 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P8 = Plot_Layer(Plot_Syns, "L8", Outline_L8_XZ, paste("L8 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P9 = Plot_Layer(Plot_Syns, "L9", Outline_L9_XZ, paste("L9 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  
  P12<-ggarrange(plotlist=list(P9,P8,P7,P6,P5,P4,P3,P2,P1), ncol =3, nrow = 3)
  ggsave(paste(DIR, PlotName,"_XZ.pdf",sep=""),
         plot = P12, device='pdf', scale = 1, width = 10, height = 12, units ="in", dpi = 500, limitsize = TRUE)
  ggsave(paste(DIR, PlotName,"_XZ.png",sep=""),
         plot = P12, device='png', scale = 1, width = 10, height = 12, units ="in", dpi = 500, limitsize = TRUE)
  
  P1b = Plot_Layer(Plot_Syns, "L1", Outline_L1_XY, paste("L1 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P2b = Plot_Layer(Plot_Syns, "L2", Outline_L2_XY, paste("L2 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P3b = Plot_Layer(Plot_Syns, "L3", Outline_L3_XY, paste("L3 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P4b = Plot_Layer(Plot_Syns, "L4", Outline_L4_XY, paste("L4 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P5b = Plot_Layer(Plot_Syns, "L5", Outline_L5_XY, paste("L5 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P6b = Plot_Layer(Plot_Syns, "L6", Outline_L6_XY, paste("L6 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P7b = Plot_Layer(Plot_Syns, "L7", Outline_L7_XY, paste("L7 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P8b = Plot_Layer(Plot_Syns, "L8", Outline_L8_XY, paste("L8 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P9b = Plot_Layer(Plot_Syns, "L9", Outline_L9_XY, paste("L9 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  
  P12b<-ggarrange(plotlist=list(P9b,P8b,P7b,P6b,P5b,P4b,P3b,P2b,P1b), ncol =3, nrow = 3)
  ggsave(paste(DIR, PlotName,"_XY.pdf",sep=""),
         plot = P12b, device='pdf', scale = 1, width = 10, height = 12, units ="in", dpi = 500, limitsize = TRUE,  useDingbats=FALSE)
  ggsave(paste(DIR, PlotName,"_XY.pdf",sep=""),
         plot = P12b, device='pdf', scale = 1, width = 10, height = 12, units ="in", dpi = 500, limitsize = TRUE,  useDingbats=FALSE)
}



###############################################################################################################################
################# Function for plotting synapse distribution of each individual FB layer ######################################



Plot_Layer <- function(Plot_Syns, Layer_To_Plot, Outline, Title, Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, View){
  
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
  P1 <- Plot_Synapse_Distribution(PlotData, View, Outline, ColorVar, col_vector, DIR, PlotName, Layer_To_Plot, TOPLOT_TF, TOPLOTLEGEND_TF)

  return(P1)
}



###############################################################################################################################
################# Generic plotting function for synapse distributions  ################################################



Plot_Synapse_Distribution <- function(PlotData, VIEW, Outline, ColorVar, col_vector, DIR, PlotName, PlotName2, TOPLOT_TF, TOPLOTLEGEND_TF, Size=2){
  
  # Get view
  if (VIEW=="XY"){
    PlotData = PlotData[!colnames(PlotData)=="Z"]
    colnames(PlotData)[colnames(PlotData)=="Y"]="Z"
    PlotData$Z=-PlotData$Z
  }
  
  # Plot just the data for now
  P1 <- ggplot() + geom_point(data=PlotData, aes(x=X, y=Z, color=ColorVar, shape=Side) ,  size=Size, alpha = 1, stroke = 0) + 
    geom_path(data=Outline, aes(x=c1, y=c2), size = 1) + xlim(c(-9000, 9000)) +  coord_fixed(ratio = 1) 

  # Save with scaled visible if TOPLOTLEGEND_TF us TRUE
  if (TOPLOTLEGEND_TF==TRUE){
    ggsave(paste(DIR, PlotName,"_",PlotName2 , "_Legend.pdf",sep=""), plot = P1, device='pdf', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE) 
  }  
  
  # Format the rest of the plot, get rid of axes
  P1 <- P1 + theme_void() + guides(fill=FALSE) + theme(legend.position = "none") +
      scale_color_manual(values=col_vector, drop=FALSE)  +  theme(
      panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent", color = NA) # get rid of legend panel bg
    )
  
  # Save plot
  if (TOPLOT_TF==TRUE){
    ggsave(paste(DIR, PlotName,"_",PlotName2, VIEW, ".pdf",sep=""), plot = P1, device='pdf', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE,  bg = "transparent")
  } 
  
  return(P1)
}



###############################################################################################################################
################# Function for getting color map of FB columns or PB glomeruli ################################################



GetColorMap <- function(Grouping){
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


PlotColumn_Types<- function(Plot_Syns, PlotName, DIR, PLOT, ID){
  Types=unique(Plot_Syns$type)
  for (nnn in 1:length(Types)){
    Temp_Type=subset(Plot_Syns, type == as.character(Types[nnn]))
    T_ID=ID[Plot_Syns$type == as.character(Types[nnn]) ]
    PlotColumn_Aves(Temp_Type, paste(PlotName, "_", Types[nnn],sep="") , DIR, PLOT, T_ID)
  }
}



PlotColumn_Aves <- function(Plot_Syns, PlotName, DIR, PLOT, T_ID){
  
  
  # get layers
  Layers_All=sort(unique(as.character(Plot_Syns$Layer)))
  
  # Make plots
  l1=Plot_Crosses(Plot_Syns, 1, DIR, PlotName, PLOT, T_ID)
  l2=Plot_Crosses(Plot_Syns, 2, DIR, PlotName, PLOT, T_ID)
  l3=Plot_Crosses(Plot_Syns, 3, DIR, PlotName, PLOT, T_ID)
  l4=Plot_Crosses(Plot_Syns, 4, DIR, PlotName, PLOT, T_ID)
  l5=Plot_Crosses(Plot_Syns, 5, DIR, PlotName, PLOT, T_ID)
  l6=Plot_Crosses(Plot_Syns, 6, DIR, PlotName, PLOT, T_ID)
  l7=Plot_Crosses(Plot_Syns, 7, DIR, PlotName, PLOT, T_ID)
  l8=Plot_Crosses(Plot_Syns, 8, DIR, PlotName, PLOT, T_ID)
  
  
  # Save plot
  PP=ggarrange(l8,l7,l6, l5, l4, l3, l2, l1,  nrow = 3, ncol = 3)
  ggsave(paste(DIR, PlotName, "_AllLayer" , ".png",sep=""), plot = PP, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE,  bg = "transparent") 

}


Plot_Crosses <- function(Plot_Syns, LayerToPlot, DIR, PlotName, PLOT, ID){
  
  # Get color vector
  col_vector=GetColorMap("FBcol")
  
  # Get mean positions for this layer
  Temp_Dist=subset(Plot_Syns, Layer == paste("L",as.character(LayerToPlot),sep="")) 
  colnames(Temp_Dist)[which(colnames(Temp_Dist) %in% c("X_Mean","Z_Mean"))] = c("X","Z")
  Temp_ID=ID[Plot_Syns$Layer == paste("L",as.character(LayerToPlot),sep="")]
  
  # Get layer outline 
  eval(parse( text= paste("Outline=Outline_L", as.character(LayerToPlot),"_XZ",sep="") ))
  
  # Get cross lines
  Horizontal=data.frame(X=c(Temp_Dist[,"Cross_XL_X"], Temp_Dist[,"Cross_XR_X"]), Y=c(Temp_Dist[,"Cross_XL_Z"], Temp_Dist[,"Cross_XR_Z"]),
                        ID=c(Temp_ID, Temp_ID), COLOR=c(as.character(Temp_Dist[,"FBcol"]), as.character(Temp_Dist[,"FBcol"])))
  Vertical=data.frame(X=c(Temp_Dist[,"Cross_YU_X"], Temp_Dist[,"Cross_YD_X"]), Y=c(Temp_Dist[,"Cross_YU_Z"], Temp_Dist[,"Cross_YD_Z"]),
                        ID=c(Temp_ID, Temp_ID), COLOR=c(as.character(Temp_Dist[,"FBcol"]), as.character(Temp_Dist[,"FBcol"])))
  

  # Plot data
  p1 <- ggplot() + geom_path(data=Outline, aes(x=c1, y=c2), size = 1) + xlim(c(-9000, 9000)) +  coord_fixed(ratio = 1) +
    geom_line(data=Horizontal, aes(x=X, y=Y, group=ID, color=COLOR)) + 
    geom_line(data=Vertical, aes(x=X, y=Y, group=ID, color=COLOR)) +
    geom_point(data=Temp_Dist, aes(x=X,y=Z, color= FBcol),size=2 ) +
    theme_void() + guides(fill=FALSE) + theme(legend.position = "none") +
    scale_color_manual(values=col_vector, drop=FALSE)  +  theme(
      panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent", color = NA), # get rid of legend panel bg
      panel.border = element_rect( fill = NULL,
                                   colour = NULL,
                                   size = NULL,
                                   linetype = NULL,
                                   color = NULL,
                                   inherit.blank = TRUE)
    ) 
                     
  

  
  if (PLOT==TRUE){
  ggsave(paste(DIR, PlotName, "_L", as.character(LayerToPlot) , ".png",sep=""), plot = p1, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE,  bg = "transparent")
  }
  return(p1)
}



###############################################################################################################################
### Function for assigning FB columns to D0 neurons, here based on input column ###############################################

D0_OutputPlot_Connectivity <- function(D0_FB_Outputs){
  
  OutputsPlot_Sub <- ggplot(D0_FB_Outputs) + theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="thistle", mid="blueviolet", high="black",
                         midpoint =0.5*max(D0_FB_Outputs$weightRelative),
                         limits=c(0,max(D0_FB_Outputs$weightRelative))) + 
    geom_tile(aes(FBcol.to, InputName,  fill=weightRelative)) + 
    facet_grid(reorder(type.from, desc(type.from)) ~ type.to, space="free", scales="free",switch="both") +
    theme(axis.text.y = element_blank()) + xlab("FB Columnar Type (POST)") + ylab("Delta0 Type (PRE)")
  
  return(OutputsPlot_Sub)
}


D0_InputPlot_Connectivity <- function(D0_FB_Inputs){
  
  InputsPlot_Sub <- ggplot(D0_FB_Inputs) + theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="thistle", mid="blueviolet", high="black",
                         midpoint =0.5*max(D0_FB_Inputs$weightRelative),
                         limits=c(0,max(D0_FB_Inputs$weightRelative))) + 
    geom_tile(aes(OutputName, FBcol.from, fill=weightRelative)) + 
    facet_grid(reorder(type.from, desc(type.from)) ~ type.to, space="free", scales="free",switch="both") +
    theme(axis.text.x = element_blank()) + ylab("FB Columnar Type (PRE)") + xlab("Delta0 Type (POST)")
  
  return(InputsPlot_Sub)
}




###############################################################################################################################
### Function for plotting mapping between PB glom and FB cols #################################################################



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


